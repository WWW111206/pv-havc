% =========================================================================
% SCRIPT: run_validation.m (Robust + Bounded + Warm-start)
% - 修正 get_all_outputs 的返回顺序接收错误
% - 用 lsqnonlin 带边界与并行有限差分（较 fsolve 更稳）
% - 首步失败时做粗网格预搜索 warm start
% =========================================================================
clear; clc; close all;
addpath('../models');

% 并行池（可选：threads 模式启动快）
if isempty(gcp('nocreate'))
    try
        parpool('local');
    catch
        % 忽略
    end
end

% --- 1. 初始化 ---
params = load_parameters();
load('../data/input_data.mat', 'time_sim_minutes', 'T_amb_sim'); % 只加载环境温度
load('../data/validation_data.mat'); % 加载实验数据

% --- 2. 运行 2000 / 3500 RPM 验证工况 ---
fprintf('--- 开始运行 2000 RPM 验证工况 ---\n');
[power_sim_2000, temp_sim_2000] = run_single_validation(2000, T_amb_sim, params);

fprintf('\n--- 开始运行 3500 RPM 验证工况 ---\n');
[power_sim_3500, temp_sim_3500] = run_single_validation(3500, T_amb_sim, params);

% --- 3. 结果对比与绘图 (复现图9) ---
time_vector_val = 0:1:(200); % 验证时长为200分钟

figure('Name', 'Model Validation Results (Figure 9 Replication)');

% 压缩机功率
subplot(2, 1, 1);
plot(time_exp_2000, power_exp_2000, 'ro', 'DisplayName', '2000rpm Experiment'); hold on;
plot(time_vector_val, power_sim_2000, 'r-', 'LineWidth', 1.5, 'DisplayName', '2000rpm Simulation');
plot(time_exp_3500, power_exp_3500, 'bs', 'DisplayName', '3500rpm Experiment');
plot(time_vector_val, power_sim_3500, 'b-', 'LineWidth', 1.5, 'DisplayName', '3500rpm Simulation');
grid on; xlabel('Time (min)'); ylabel('Compressor input power (W)');
title('Compressor Input Power Comparison'); legend('Location','best');

% 水温
subplot(2, 1, 2);
plot(time_exp_2000, temp_exp_2000, 'ro', 'DisplayName', '2000rpm Experiment'); hold on;
plot(time_vector_val, temp_sim_2000, 'r-', 'LineWidth', 1.5, 'DisplayName', '2000rpm Simulation');
plot(time_exp_3500, temp_exp_3500, 'bs', 'DisplayName', '3500rpm Experiment');
plot(time_vector_val, temp_sim_3500, 'b-', 'LineWidth', 1.5, 'DisplayName', '3500rpm Simulation');
grid on; xlabel('Time (min)'); ylabel('Water temperature (°C)');
title('Water Temperature Comparison'); legend('Location','best');

% --- 4. 计算 MAE ---
mae_power_2000 = mean(abs(interp1(time_vector_val, power_sim_2000, time_exp_2000) - power_exp_2000));
mae_temp_2000  = mean(abs(interp1(time_vector_val, temp_sim_2000, time_exp_2000) - temp_exp_2000));
mae_power_3500 = mean(abs(interp1(time_vector_val, power_sim_3500, time_exp_3500) - power_exp_3500));
mae_temp_3500  = mean(abs(interp1(time_vector_val, temp_sim_3500, time_exp_3500) - temp_exp_3500));

fprintf('\n--- MAE Results ---\n');
fprintf('2000 RPM Power MAE: %.2f W\n', mae_power_2000);
fprintf('2000 RPM Temp  MAE: %.2f °C\n', mae_temp_2000);
fprintf('3500 RPM Power MAE: %.2f W\n', mae_power_3500);
fprintf('3500 RPM Temp  MAE: %.2f °C\n', mae_temp_3500);


% --- 辅助函数：用于执行单次验证仿真 ---
function [power_results, temp_results] = run_single_validation(rpm, T_amb_profile, params)
    num_steps = 201; % 0..200 min
    power_results = zeros(1, num_steps);
    temp_results  = zeros(1, num_steps);
    temp_results(1) = 25.0; % 初始水温

    x0 = [12e5; 2.5e5]; % 初始 [p_con, p_eva] Pa
    lb = [ 2e5; 0.5e5];
    ub = [40e5; 20e5];

    opts = optimoptions('lsqnonlin', ...
        'Display','off', ...
        'UseParallel', true, ...
        'FiniteDifferenceType','forward', ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-8, ...
        'MaxIterations', 150);

    for i = 1:num_steps
        fprintf('  - Time: %d min\n', i-1);

        % 边界条件
        bcs.T_amb = T_amb_profile(i) + 273.15;      % K
        bcs.T_air_chamber = temp_results(i) + 273.15; % 假设箱内空气≈水温

        % 恒速控制
        control_input.mode  = 'CSP';
        control_input.omega = rpm * 2*pi/60;

        fun = @(x) system_residuals(x, bcs, params, control_input);
        [x_sol, ~, residual, exitflag] = lsqnonlin(fun, x0, lb, ub, opts);

        if exitflag <= 0 && i == 1
            % 首步失败：做一次粗网格 warm start，再求解一次
            x0_try = coarse_warm_start(bcs, params, control_input, lb, ub);
            [x_sol, ~, residual, exitflag] = lsqnonlin(fun, x0_try, lb, ub, opts);
        end

        if exitflag > 0
            x0 = x_sol; % 更新下一步初值
            % 取输出（修正：按正确顺序接收返回值）
            [comp_out, cond_out, cap_out, eva_out] = get_all_outputs(x_sol, bcs, params, control_input); %#ok<ASGLU>
            % 功率
            power_results(i) = comp_out.P_in_elec;
            % 一阶水温更新
            if i < num_steps
                dt = 60;             % s
                mass_water = 4.5;    % kg (9 * 0.5 kg)
                cp_water   = 4186;   % J/kg.K
                Q_cooling  = max(abs(eva_out.Q_total), 0);
                delta_T    = Q_cooling * dt / (mass_water * cp_water);
                temp_results(i+1) = temp_results(i) - delta_T;
            end
        else
            fprintf('    Solver failed at step %d (||F||=%.3e). Using previous values.\n', i, norm(residual));
            if i > 1
                power_results(i) = power_results(i-1);
                if i < num_steps, temp_results(i+1) = temp_results(i); end
            else
                power_results(i) = NaN;
                if i < num_steps, temp_results(i+1) = temp_results(i); end
            end
        end
    end
end

function x0 = coarse_warm_start(bcs, params, control_input, lb, ub)
    % 在物理合理范围做一个粗网格搜索，找 ||F|| 最小的点作为 warm start
    p_con_vec = linspace(8e5, 18e5, 6);   % 8..18 bar
    p_eva_vec = linspace(1.5e5, 4e5, 6);  % 1.5..4 bar
    bestJ = inf; x0 = [min(ub(1),max(lb(1),12e5)); min(ub(2),max(lb(2),2.5e5))];
    fun = @(x) system_residuals(x, bcs, params, control_input);

    for pc = p_con_vec
        for pe = p_eva_vec
            if pc <= pe, continue; end
            r = fun([pc; pe]);
            J = sum(r.^2);
            if isfinite(J) && J < bestJ
                bestJ = J; x0 = [pc; pe];
            end
        end
    end
end

function [comp_out, cond_out, cap_out, eva_out] = get_all_outputs(x, bcs, params, control_input)
    p_con = x(1); p_eva = x(2);
    h_comp_in = get_fluid_props('H','P',p_eva,'Q',1,'R134a');
    comp_out = model_compressor(control_input, p_eva, p_con, h_comp_in, params);
    cond_out = model_heat_exchanger('condenser', comp_out.mdot_ref, p_con, comp_out.h_ref_out, bcs.T_amb, params);
    cap_out  = model_capillary(p_con, cond_out.h_out_ref, p_eva, params);
    eva_out  = model_heat_exchanger('evaporator', cap_out.mdot_ref, p_eva, cap_out.h_out, bcs.T_air_chamber, params);
end

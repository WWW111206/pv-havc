% =========================================================================
% SCRIPT: run_validation.m (Corrected Version)
% DESCRIPTION:
% This script runs the simulation under constant compressor speed conditions
% (2000 rpm and 3500 rpm) to validate the model against the experimental
% data from Figure 9 of the source paper.
% =========================================================================
clear; clc; close all;
addpath('../models');

% --- 1. 初始化 ---
params = load_parameters();
load('../data/input_data.mat', 'time_sim_minutes', 'T_amb_sim'); % 只加载环境温度
load('../data/validation_data.mat'); % 加载实验数据

% --- 2. 运行 2000 RPM 验证工况 ---
fprintf('--- 开始运行 2000 RPM 验证工况 ---\n');
[power_sim_2000, temp_sim_2000] = run_single_validation(2000, T_amb_sim, params);

% --- 3. 运行 3500 RPM 验证工况 ---
fprintf('\n--- 开始运行 3500 RPM 验证工况 ---\n');
[power_sim_3500, temp_sim_3500] = run_single_validation(3500, T_amb_sim, params);

% --- 4. 结果对比与绘图 (复现图9) ---
time_vector_val = 0:1:(200); % 验证时长为200分钟

figure('Name', 'Model Validation Results (Figure 9 Replication)');

% 绘制压缩机功率
subplot(2, 1, 1);
plot(time_exp_2000, power_exp_2000, 'ro', 'DisplayName', '2000rpm Experiment');
hold on;
plot(time_vector_val, power_sim_2000, 'r-', 'LineWidth', 1.5, 'DisplayName', '2000rpm Simulation');
plot(time_exp_3500, power_exp_3500, 'bs', 'DisplayName', '3500rpm Experiment');
plot(time_vector_val, power_sim_3500, 'b-', 'LineWidth', 1.5, 'DisplayName', '3500rpm Simulation');
grid on;
xlabel('Time (min)');
ylabel('Compressor input power (W)');
title('Compressor Input Power Comparison');
legend('Location', 'best');

% 绘制水温
subplot(2, 1, 2);
plot(time_exp_2000, temp_exp_2000, 'ro', 'DisplayName', '2000rpm Experiment');
hold on;
plot(time_vector_val, temp_sim_2000, 'r-', 'LineWidth', 1.5, 'DisplayName', '2000rpm Simulation');
plot(time_exp_3500, temp_exp_3500, 'bs', 'DisplayName', '3500rpm Experiment');
plot(time_vector_val, temp_sim_3500, 'b-', 'LineWidth', 1.5, 'DisplayName', '3500rpm Simulation');
grid on;
xlabel('Time (min)');
ylabel('Water temperature (°C)');
title('Water Temperature Comparison');
legend('Location', 'best');

% --- 5. 计算平均绝对误差 (MAE) ---
mae_power_2000 = mean(abs(interp1(time_vector_val, power_sim_2000, time_exp_2000) - power_exp_2000));
mae_temp_2000 = mean(abs(interp1(time_vector_val, temp_sim_2000, time_exp_2000) - temp_exp_2000));
mae_power_3500 = mean(abs(interp1(time_vector_val, power_sim_3500, time_exp_3500) - power_exp_3500));
mae_temp_3500 = mean(abs(interp1(time_vector_val, temp_sim_3500, time_exp_3500) - temp_exp_3500));

fprintf('\n--- MAE Results ---\n');
fprintf('2000 RPM Power MAE: %.2f W (Paper: 1.12 W)\n', mae_power_2000);
fprintf('2000 RPM Temp MAE:  %.2f °C (Paper: 0.35 °C)\n', mae_temp_2000);
fprintf('3500 RPM Power MAE: %.2f W (Paper: 1.40 W)\n', mae_power_3500);
fprintf('3500 RPM Temp MAE:  %.2f °C (Paper: 0.28 °C)\n', mae_temp_3500);


% --- 辅助函数：用于执行单次验证仿真 ---
function [power_results, temp_results] = run_single_validation(rpm, T_amb_profile, params)
    num_steps = 201; % 0 to 200 minutes
    
    % 初始化结果存储
    power_results = zeros(1, num_steps);
    temp_results = zeros(1, num_steps);
    temp_results(1) = 25.0; % 初始水温
    
    x0 = [12e5; 2.5e5]; % 初始压力猜测值 [p_con, p_eva] in Pa

    for i = 1:num_steps
        fprintf('  - Time: %d min\n', i-1);
        
        % 边界条件
        bcs.T_amb = T_amb_profile(i) + 273.15; % K
        bcs.T_air_chamber = temp_results(i) + 273.15; % 假设箱内空气温度等于水温
        
        % 控制输入：恒定转速
        control_input.mode = 'CSP';
        control_input.omega = rpm * 2 * pi / 60; % rad/s
        
        % 求解
        residual_fun = @(x) system_residuals(x, bcs, params, control_input);
        options = optimoptions('fsolve', 'Display', 'none');
        [x_sol, ~, exitflag] = fsolve(residual_fun, x0, options);
        
        if exitflag > 0
            x0 = x_sol; % 更新下一次迭代的初始猜测值
            
            % 重新计算一次以获取输出
            [~, comp_out, ~, eva_out] = get_all_outputs(x_sol, bcs, params, control_input);
            
            % 存储功率
            power_results(i) = comp_out.P_in_elec;
            
            % 简单的一阶模型更新水温
            if i < num_steps
                dt = 60; % 60 seconds
                mass_water = 4.5; % 9 * 0.5 kg
                cp_water = 4186; % J/kg.K
                Q_cooling = abs(eva_out.Q_total);
                delta_T = Q_cooling * dt / (mass_water * cp_water);
                temp_results(i+1) = temp_results(i) - delta_T;
            end
        else
            fprintf('    Solver failed at step %d. Using previous values.\n', i);
            % <-- 修正: 增加判断，防止在第一个时间步(i=1)出错
            if i > 1
                power_results(i) = power_results(i-1);
                if i < num_steps, temp_results(i+1) = temp_results(i); end
            else % 如果第一个时间步就失败了
                power_results(i) = NaN; % 标记为无效值
                if i < num_steps, temp_results(i+1) = temp_results(i); end
            end
        end
    end
end

function [comp_out, cond_out, cap_out, eva_out] = get_all_outputs(x, bcs, params, control_input)
    p_con = x(1); p_eva = x(2);
    h_comp_in = get_fluid_props('H', 'P', p_eva, 'Q', 1, 'R134a');
    comp_out = model_compressor(control_input, p_eva, p_con, h_comp_in, params);
    cond_out = model_heat_exchanger('condenser', comp_out.mdot_ref, p_con, comp_out.h_ref_out, bcs.T_amb, params);
    cap_out = model_capillary(p_con, cond_out.h_out_ref, p_eva, params);
    eva_out = model_heat_exchanger('evaporator', cap_out.mdot_ref, p_eva, cap_out.h_out, bcs.T_air_chamber, params);
end
clear; clc; close all;
addpath('../models');

% 1) 参数与输入数据
params = load_parameters();
S = load('../data/input_data.mat'); % T_amb_sim(°C), Ira_sim(W/m2), time_sim_minutes(min)
time_vector = S.time_sim_minutes * 60; % s
num_steps   = numel(time_vector);

% 2) 结果存储与初值
p_con_results = zeros(1, num_steps);
p_eva_results = zeros(1, num_steps);
x0 = [10e5; 2e5];

for i = 1:num_steps
    fprintf('正在计算时间步 %d / %d\n', i, num_steps);

    bcs.T_amb = S.T_amb_sim(i) + 273.15; % K
    bcs.Ira   = S.Ira_sim(i);            % W/m^2

    % 简化：箱内空气温度取环境温度（若有更精细的箱体/水箱模型可替换）
    bcs.T_air_chamber = bcs.T_amb;

    % PV 温度（若设置了 NOCT，可替换为 NOCT 模型）
    if isfield(params.pv, 'T_const')
        T_PV = params.pv.T_const;
    else
        T_PV = bcs.T_amb;
    end

    % MPPT 目标功率
    pv_out = model_pv_panel(bcs.Ira, T_PV, params);
    eta    = params.electronics.eta_mppt;
    control_input.mode = 'MPPT';
    control_input.P_in = eta * pv_out.P_mpp;

    residual_fun = @(x) system_residuals(x, bcs, params, control_input);
    options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-6);

    [x_sol, ~, exitflag] = fsolve(residual_fun, x0, options);

    if exitflag > 0
        p_con_results(i) = x_sol(1);
        p_eva_results(i) = x_sol(2);
        x0 = x_sol;
    else
        fprintf('时间步 %d 求解失败，exitflag = %d\n', i, exitflag);
        if i > 1
            p_con_results(i) = p_con_results(i-1);
            p_eva_results(i) = p_eva_results(i-1);
        end
    end
end

% 3) 简单可视化
figure;
plot(time_vector/3600, p_con_results/1e5, 'b-', 'LineWidth',1.2); hold on;
plot(time_vector/3600, p_eva_results/1e5, 'r-', 'LineWidth',1.2);
xlabel('时间 (小时)'); ylabel('压力 (bar)');
legend('冷凝压力','蒸发压力'); grid on; title('日内运行压力');

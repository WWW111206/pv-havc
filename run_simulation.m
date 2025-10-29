clear; clc; close all;
addpath('../models'); % 将模型文件夹添加到路径

% 1. 初始化
params = load_parameters();
load('../data/input_data.mat'); % 加载T_amb_sim, Ira_sim, time_sim_minutes

% 2. 设置仿真参数
time_vector = time_sim_minutes * 60; % 转换为秒
dt = 60; % 时间步长
num_steps = length(time_vector);

% 3. 初始化结果存储数组和初始猜测值
p_con_results = zeros(1, num_steps);
p_eva_results = zeros(1, num_steps);
x0 = [10e5; 2e5]; % 初始压力猜测值 [p_con, p_eva] in Pa

% 4. 时间循环
for i = 1:num_steps
    fprintf('正在计算时间步 %d / %d\n', i, num_steps);

    % 获取当前边界条件
    bcs.T_amb = T_amb_sim(i) + 273.15; % K
    bcs.Ira = Ira_sim(i); % W/m^2

    % 设置控制输入 (此处以MPPT为例)
    pv_out = model_pv_panel(bcs.Ira, params.pv.T_const, params);
    control_input.mode = 'MPPT';
    control_input.P_in = 0.98 * pv_out.P_mpp;

    % 定义匿名函数句柄
    residual_fun = @(x) system_residuals(x, bcs, params, control_input);

    % 设置求解器选项
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-6);

    % 求解
    [x_sol, ~, exitflag] = fsolve(residual_fun, x0, options);

    if exitflag > 0
        % 存储结果
        p_con_results(i) = x_sol(1);
        p_eva_results(i) = x_sol(2);
        % 更新下一次迭代的初始猜测值
        x0 = x_sol;
    else
        fprintf('时间步 %d 求解失败，exitflag = %d\n', i, exitflag);
        % 保持上一步结果
        if i > 1
            p_con_results(i) = p_con_results(i-1);
            p_eva_results(i) = p_eva_results(i-1);
        end
    end
end

% 5. 结果可视化
figure;
plot(time_vector/3600, p_con_results/1e5);
xlabel('时间 (小时)');
ylabel('冷凝压力 (bar)');
title('冷凝压力随时间变化');
grid on;
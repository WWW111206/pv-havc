% =========================================================================
% SCRIPT: prepare_input_data.m (Final Corrected and Populated Version)
% DATE: 2025-10-29
% =========================================================================

clear; clc; close all;

disp('开始处理和生成仿真输入数据...');

% --- 1. 定义从图表中手动提取的原始数据点 ---
% 注意：下面给出的是示例/占位的数字化点（7:00 到 17:00）。
% 请用你从论文图10中实际提取的点替换这些数组，或从 CSV 文件读取。

% X轴: 时间 (小时制，从 7 表示 7:00 到 17 表示 17:00)
time_hours_ira_raw = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];

% Y轴: 太阳辐照度 (W/m^2)，示例值（早晚为 0，中午峰值）
irradiance_raw = [0, 120, 300, 520, 700, 900, 820, 600, 380, 120, 0];

% 校验数组一致性
if numel(time_hours_ira_raw) ~= numel(irradiance_raw)
    error('时间数据点与辐照度数据点数量不匹配，请检查 time_hours_ira_raw 与 irradiance_raw.');
end

% --- 2. 创建统一的仿真时间轴 ---
% 仿真时间范围与图10一致，从早上 7:00 到下午 17:00
total_duration_minutes = (time_hours_ira_raw(end) - time_hours_ira_raw(1)) * 60;
time_sim_minutes = 0:1:total_duration_minutes; % 仿真时间向量 (单位: 分钟)

% --- 3. 数据转换与插值 ---
% a) 处理太阳辐照度数据
% 将小时制的时间转换为相对于仿真开始时间（7:00）的分钟数
time_minutes_ira_rel = (time_hours_ira_raw - time_hours_ira_raw(1)) * 60;

% 使用线性插值方法，根据仿真时间轴生成平滑的太阳辐照度曲线
Ira_sim = interp1(time_minutes_ira_rel, irradiance_raw, time_sim_minutes, 'linear', 'extrap');

% 清理数据：确保辐照度非负
Ira_sim(Ira_sim < 0) = 0;

% b) 处理环境温度数据
% 使用论文中给出的典型值及图8的形状来生成一个示例环境温度曲线
mid_point_sim = total_duration_minutes / 2;
% 对应 7:00, ~9:30, 12:00, ~14:30, 17:00（单位：分钟相对）
time_amb_full_day = [0, mid_point_sim/2, mid_point_sim, mid_point_sim*1.5, total_duration_minutes];
temp_amb_full_day = [25.0, 27.0, 29.0, 30.0, 26.0]; % (°C)，可替换为更精确的点

% 对假设的温度点进行平滑插值
T_amb_sim = interp1(time_amb_full_day, temp_amb_full_day, time_sim_minutes, 'pchip');

% --- 4. 可视化检查 ---
figure('Name', 'Interpolated Simulation Inputs', 'NumberTitle', 'off');

% 绘制太阳辐照度
subplot(2, 1, 1);
plot(time_hours_ira_raw, irradiance_raw, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '原始数据点');
hold on;
plot(time_sim_minutes/60 + time_hours_ira_raw(1), Ira_sim, 'b-', 'LineWidth', 1.5, 'DisplayName', '插值后曲线');
title('太阳辐照度 (源自图10)');
xlabel('一天中的时间 (小时)');
ylabel('辐照度 (W/m^2)');
grid on;
legend('Location','northwest');
xlim([time_hours_ira_raw(1), time_hours_ira_raw(end)]);

% 绘制环境温度
subplot(2, 1, 2);
plot(time_amb_full_day/60 + time_hours_ira_raw(1), temp_amb_full_day, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '关键温度点');
hold on;
plot(time_sim_minutes/60 + time_hours_ira_raw(1), T_amb_sim, 'b-', 'LineWidth', 1.5, 'DisplayName', '插值后曲线');
title('环境温度 (基于图8和论文描述扩展)');
xlabel('一天中的时间 (小时)');
ylabel('温度 (°C)');
grid on;
legend('Location','northwest');
xlim([time_hours_ira_raw(1), time_hours_ira_raw(end)]);

% --- 5. 保存结果 ---
% 将处理好的数据保存到 ../data 文件夹下的 .mat 文件中
data_dir = fullfile('..', 'data');
if ~exist(data_dir, 'dir')
   mkdir(data_dir);
end
output_filename = fullfile(data_dir, 'input_data.mat');

save(output_filename, 'time_sim_minutes', 'T_amb_sim', 'Ira_sim');

fprintf('\n数据处理完成！\n');
fprintf('仿真输入数据已成功保存至: %s\n', output_filename);
% =========================================================================
% SCRIPT: prepare_validation_data.m
% DESCRIPTION:
% This script digitizes and saves the experimental data from Figure 9 of
% the source paper for model validation purposes.
% =========================================================================
clear; clc; close all;

disp('正在处理和保存用于验证的实验数据...');

% --- 数据来源: 论文图9 ---
% 注意：下面给出的 time_exp_2000 和 time_exp_3500 为示例时间点（单位：min）。
% 请用你实际数字化得到的时间向量替换这些示例数据。

% 2000 rpm 实验数据
% 示例时间点（分钟），长度应与 power_exp_2000 / temp_exp_2000 一致
time_exp_2000 = [0, 5, 10, 15, 20, 25, 30, 35, 40]; % min
power_exp_2000 = [30.5, 30.0, 29.5, 29.0, 28.8, 28.6, 28.4, 28.2, 28.0]; % W
temp_exp_2000 = [25.0, 19.0, 14.5, 11.0, 8.5, 6.5, 5.0, 4.0, 3.2]; % degC

% 3500 rpm 实验数据
% 示例时间点（分钟），长度应与 power_exp_3500 / temp_exp_3500 一致
time_exp_3500 = [0, 5, 10, 15, 20, 25, 30, 35, 40]; % min
power_exp_3500 = [70.5, 69.5, 68.0, 66.5, 65.0, 63.8, 62.8, 62.0, 61.2]; % W
temp_exp_3500 = [25.0, 15.5, 9.0, 5.0, 2.5, 0.8, -0.5, -1.5, -2.3]; % degC

% --- 校验数组一致性 ---
if ~(numel(time_exp_2000) == numel(power_exp_2000) && numel(time_exp_2000) == numel(temp_exp_2000))
    error('长度不匹配：time_exp_2000, power_exp_2000, temp_exp_2000 的长度必须相等。');
end
if ~(numel(time_exp_3500) == numel(power_exp_3500) && numel(time_exp_3500) == numel(temp_exp_3500))
    error('长度不匹配：time_exp_3500, power_exp_3500, temp_exp_3500 的长度必须相等。');
end

% --- 可选：将行向量转为列向量（更常用的时间序列形式） ---
time_exp_2000 = time_exp_2000(:);
power_exp_2000 = power_exp_2000(:);
temp_exp_2000  = temp_exp_2000(:);

time_exp_3500 = time_exp_3500(:);
power_exp_3500 = power_exp_3500(:);
temp_exp_3500  = temp_exp_3500(:);

% --- 保存到 ../data 文件夹 ---
data_dir = fullfile('..', 'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
output_filename = fullfile(data_dir, 'validation_data.mat');

save(output_filename, 'time_exp_2000', 'power_exp_2000', 'temp_exp_2000', ...
                      'time_exp_3500', 'power_exp_3500', 'temp_exp_3500');

fprintf('验证数据已成功保存至: %s\n', output_filename);

% --- 简单可视化（可选） ---
figure('Name','Validation Data','NumberTitle','off');

subplot(2,2,1);
plot(time_exp_2000, power_exp_2000, 'o-');
title('Power @ 2000 rpm');
xlabel('Time (min)'); ylabel('Power (W)');
grid on;

subplot(2,2,2);
plot(time_exp_2000, temp_exp_2000, 'o-');
title('Temperature @ 2000 rpm');
xlabel('Time (min)'); ylabel('Temp (°C)');
grid on;

subplot(2,2,3);
plot(time_exp_3500, power_exp_3500, 'o-');
title('Power @ 3500 rpm');
xlabel('Time (min)'); ylabel('Power (W)');
grid on;

subplot(2,2,4);
plot(time_exp_3500, temp_exp_3500, 'o-');
title('Temperature @ 3500 rpm');
xlabel('Time (min)'); ylabel('Temp (°C)');
grid on;

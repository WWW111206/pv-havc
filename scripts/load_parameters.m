function params = load_parameters()
    % 系统参数（SI）

    % --- 压缩机 (BD35F) ---
    params.compressor.R  = 4.4;       % Ohm
    params.compressor.Ke = 0.12;      % V*s/rad
    params.compressor.La = 3.2e-4;    % H
    params.compressor.piston_D = 15.0e-3; % m
    params.compressor.piston_S = 11.5e-3; % m

    % --- 冷凝器 ---
    params.condenser.L    = 4.0;      % m
    params.condenser.D_in = 5.0e-3;   % m
    params.condenser.D_out= 6.0e-3;   % m

    % --- 毛细管 ---
    params.capillary.L    = 1.8;      % m
    params.capillary.D_in = 1.0e-3;   % m

    % --- 蒸发器 ---
    params.evaporator.L    = 10.0;    % m
    params.evaporator.D_in = 7.0e-3;  % m
    params.evaporator.D_out= 8.0e-3;  % m
    params.evaporator.superheat_K = 5.0; % 目标过热度(可调)

    % --- 光伏板 @STC ---
    params.pv.U_oc0 = 21.6;  % V
    params.pv.I_sc0 = 7.72;  % A
    params.pv.U_m0  = 17.5;  % V
    params.pv.I_m0  = 6.85;  % A
    % PV 温度假设（若无 NOCT 模型）
    params.pv.T_const = 45 + 273.15; % K，可按需替换为 NOCT 估算

    % --- 电力电子 ---
    params.electronics.eta_mppt = 0.98;

    disp('系统参数加载完成。');
end

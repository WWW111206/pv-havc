function params = load_parameters()
    % 该函数返回一个包含所有系统参数的结构体
    % 所有单位均为国际单位制 (SI)

    % --- 压缩机 (BD35F) - 来自论文表1 ---
    params.compressor.R = 4.4; % 线路电阻 (Ohm)
    params.compressor.Ke = 0.12; % 反电动势系数 (Wb or V*s/rad)
    params.compressor.La = 3.2e-4; % 电感 (H)
    params.compressor.piston_D = 15.0e-3; % 活塞直径 (m)
    params.compressor.piston_S = 11.5e-3; % <-- 修正: 添加缺失的活塞冲程参数 (m)

    % --- 冷凝器 ---
    params.condenser.L = 4.0; % 长度 (m)
    params.condenser.D_in = 5.0e-3; % 内径 (m)
    params.condenser.D_out = 6.0e-3; % 外径 (m)

    % --- 毛细管 ---
    params.capillary.L = 1.8; % 长度 (m)
    params.capillary.D_in = 1.0e-3; % 内径 (m)

    % --- 蒸发器 ---
    params.evaporator.L = 10.0; % 长度 (m)
    params.evaporator.D_in = 7.0e-3; % 内径 (m)
    params.evaporator.D_out = 8.0e-3; % 外径 (m)

    % --- 光伏板 (Polysilicon PV panel @STC) ---
    params.pv.U_oc0 = 21.6; % 开路电压 (V)
    params.pv.I_sc0 = 7.72; % 短路电流 (A)
    params.pv.U_m0 = 17.5; % 最大功率点电压 (V)
    params.pv.I_m0 = 6.85; % 最大功率点电流 (A)
    
    disp('系统参数加载完成。');
end
function pv_out = model_pv_panel(Ira_PV, T_PV, params)
% 简化 PV 板模型（论文公式(3)-(8)风格）
% 保护：
% - 低辐照度直接输出 0（避免 log 的负值域）
% - 对数项加入下限保护

    Isc0 = params.pv.I_sc0; Uoc0 = params.pv.U_oc0;
    Im0  = params.pv.I_m0;  Um0  = params.pv.U_m0;

    Ira0 = 1000;
    T0   = 25 + 273.15;

    a = 0.0025; b = 0.5; c = 0.00288;

    if Ira_PV <= 1 % 几乎无光照
        pv_out.P_mpp = 0;
        pv_out.U_m   = 0;
        pv_out.I_m   = 0;
        return;
    end

    Grel = Ira_PV / Ira0;
    dT   = T_PV - T0;

    Isc = Isc0 * Grel * (1 + a * dT);

    % 对数项保护
    arg = 1 + b * (Grel - 1);
    arg = max(arg, 1e-6);

    Uoc = Uoc0 * (1 - c * dT) * log(arg);
    Um  = Um0  * (1 - c * dT) * log(arg);
    Im  = Im0  * Grel * (1 + a * dT);

    Um  = max(Um, 0);
    Im  = max(Im, 0);

    pv_out.P_mpp = Um * Im;
    pv_out.U_m   = Um;
    pv_out.I_m   = Im;
end

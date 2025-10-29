function pv_out = model_pv_panel(Ira_PV, T_PV, params)
    % 实现论文方程(5)-(8)和(3)-(4)
    % 输入: Ira_PV (W/m^2), T_PV (K)
    % 输出: 包含P_mpp, U_m, I_m的结构体

    % 从params结构体中获取标准测试条件下的参数
    Isc0 = params.pv.I_sc0; Uoc0 = params.pv.U_oc0;
    Im0 = params.pv.I_m0; Um0 = params.pv.U_m0;
    Ira_PV0 = 1000; T_PV0 = 25 + 273.15;

    % 典型值 a, b, c
    a = 0.0025; b = 0.5; c = 0.00288;

    % 方程 (5)-(8)
    Isc = Isc0 * (Ira_PV / Ira_PV0) * (1 + a * (T_PV - T_PV0));
    Uoc = Uoc0 * (1 - c * (T_PV - T_PV0)) * log(1 + b * (Ira_PV / Ira_PV0 - 1));
    Im = Im0 * (Ira_PV / Ira_PV0) * (1 + a * (T_PV - T_PV0));
    Um = Um0 * (1 - c * (T_PV - T_PV0)) * log(1 + b * (Ira_PV / Ira_PV0 - 1));

    % 输出
    pv_out.P_mpp = Um * Im;
    pv_out.U_m = Um;
    pv_out.I_m = Im;
end

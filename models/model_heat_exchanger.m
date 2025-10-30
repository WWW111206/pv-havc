function hx_out = model_heat_exchanger(type, mdot_ref, p_ref, h_in_ref, T_air_in, params)
% 1D 稳态有限体积换热器模型（冷凝器/蒸发器）
% 修复：
% - 以制冷剂侧换热量 dQ_ref = alpha_ref*A_in*(T_ref - T_wall) 为准，修正符号
% - i=1 时的异常分支不再访问 T_wall(i-1)
% - 属性名与 CoolProp 键兼容（通过 get_fluid_props 映射）

fluid = 'R134a';
N = 50;

% 几何
if strcmpi(type,'condenser')
    L = params.condenser.L;  D_in = params.condenser.D_in;  D_out = params.condenser.D_out;
else
    L = params.evaporator.L; D_in = params.evaporator.D_in; D_out = params.evaporator.D_out;
end

dL   = L / N;
A_in = pi * D_in * dL;
A_out= pi * D_out * dL;

% 边界：无流量时直接返回
if mdot_ref < 1e-9
    hx_out.h_out_ref  = h_in_ref;
    hx_out.Q_total    = 0;
    hx_out.T_wall_avg = T_air_in;
    return;
end

% 初始化
h_ref  = zeros(N+1,1); h_ref(1) = h_in_ref;
T_wall = zeros(N,1);

for i = 1:N
    h_out_guess = h_ref(i);

    for iter = 1:30
        h_avg = 0.5*(h_ref(i) + h_out_guess);

        try
            T_ref_avg = get_fluid_props('T','P',p_ref,'H',h_avg,fluid);
            x_avg     = get_fluid_props('Q','P',p_ref,'H',h_avg,fluid); %#ok<NASGU> (可按需使用)
        catch ME
            warning('CoolProp失败 @segment %d: %s。沿用上一段结果。', i, ME.message);
            h_ref(i+1) = h_ref(i);
            T_wall(i)  = (i>1) * T_wall(i-1) + (i==1) * T_air_in;
            break;
        end

        % 计算换热系数
        if strcmpi(type,'condenser')
            alpha_ref = calc_alpha_condenser_ref(p_ref, h_avg, mdot_ref, D_in, params);
            alpha_air = calc_alpha_condenser_air(T_air_in, T_ref_avg, D_out, params);
        else
            T_wall_guess = 0.5*(T_ref_avg + T_air_in);
            q_guess      = calc_alpha_evaporator_air(T_air_in, T_wall_guess, L, params) ...
                         * (T_wall_guess - T_air_in);
            alpha_ref = calc_alpha_evaporator_ref(p_ref, h_avg, mdot_ref, D_in, q_guess, params);
            alpha_air = calc_alpha_evaporator_air(T_air_in, T_wall_guess, L, params);
        end

        % 墙面能量平衡：alpha_ref*A_in*(T_ref-T_wall) = alpha_air*A_out*(T_wall-T_air)
        T_wall_new = (alpha_ref*A_in*T_ref_avg + alpha_air*A_out*T_air_in) ...
                   / (alpha_ref*A_in + alpha_air*A_out);

        % 使用制冷剂侧换热量（正值表示制冷剂向外放热）
        dQ_ref = alpha_ref * A_in * (T_ref_avg - T_wall_new);

        % 焓更新（冷凝器：dQ_ref>0 → 焓降低；蒸发器：dQ_ref<0 → 焓升高）
        h_out_new = h_ref(i) - dQ_ref / mdot_ref;

        if abs(h_out_new - h_out_guess) < 0.5 % J/kg 收敛阈值
            h_out_guess = h_out_new;
            T_wall(i)   = T_wall_new;
            break;
        end
        h_out_guess = h_out_new;
    end

    h_ref(i+1) = h_out_guess;
    if T_wall(i) == 0
        T_wall(i) = (i>1) * T_wall(i-1) + (i==1) * T_air_in;
    end
end

hx_out.h_out_ref  = h_ref(N+1);
hx_out.Q_total    = mdot_ref * (h_ref(N+1) - h_ref(1));
hx_out.T_wall_avg = mean(T_wall);

end

% --------------------- 辅助相关关联式 ---------------------

function alpha = calc_alpha_condenser_ref(p, h, mdot, D, params)
    fluid = 'R134a';
    x = get_fluid_props('Q','P',p,'H',h,fluid);

    if x>=0 && x<=1
        p_crit = get_fluid_props('Pcrit','',0,'',0,fluid);
        p_r = p / p_crit;

        k_L  = get_fluid_props('CONDUCTIVITY','P',p,'Q',0,fluid);
        mu_L = get_fluid_props('VISCOSITY','P',p,'Q',0,fluid);
        Pr_L = get_fluid_props('PRANDTL','P',p,'Q',0,fluid);

        G    = mdot / (pi*(D/2)^2);
        Re_L = G * D / mu_L;

        Nu_i = 0.023 * Re_L^0.8 * Pr_L^0.4;
        alpha_i = Nu_i * k_L / D;

        alpha = alpha_i * ((1-x)^0.8 + (3.8*x^0.76*(1-x)^0.04)/(p_r^0.38));
    else
        T  = get_fluid_props('T','P',p,'H',h,fluid);
        mu = get_fluid_props('VISCOSITY','T',T,'P',p,fluid);
        k  = get_fluid_props('CONDUCTIVITY','T',T,'P',p,fluid);
        Pr = get_fluid_props('PRANDTL','T',T,'P',p,fluid);

        G  = mdot / (pi*(D/2)^2);
        Re = G * D / mu;

        if Re < 2300, Nu = 3.66; else, Nu = 0.023*Re^0.8*Pr^0.3; end
        alpha = Nu * k / D;
    end
end

function alpha = calc_alpha_evaporator_ref(p, h, mdot, D, q, params)
    fluid = 'R134a';
    x = get_fluid_props('Q','P',p,'H',h,fluid);

    if x>0 && x<1
        rho_L = get_fluid_props('D','P',p,'Q',0,fluid);
        rho_G = get_fluid_props('D','P',p,'Q',1,fluid);
        mu_L  = get_fluid_props('VISCOSITY','P',p,'Q',0,fluid);
        mu_G  = get_fluid_props('VISCOSITY','P',p,'Q',1,fluid);
        h_L   = get_fluid_props('H','P',p,'Q',0,fluid);
        h_G   = get_fluid_props('H','P',p,'Q',1,fluid);
        Lat   = h_G - h_L;

        Bo  = q / max(mdot*Lat, 1e-6);
        Xtt = (x/(1-x))^0.9 * (rho_G/rho_L)^0.5 * (mu_L/mu_G)^0.1;

        a1  = 1 + 24000*Bo^1.16 + 1.37*Xtt^-0.86;
        G   = mdot / (pi*(D/2)^2);
        ReL = G * D / mu_L;
        a2  = (1 + 1.15e6 * a1^2 * ReL^1.17)^-1;

        M   = 102.03; % g/mol for R134a
        alpha = a1 * a2 * (p/1e5)^0.24 * M^-0.5 * max(q,1e-6)^0.7;
    else
        T  = get_fluid_props('T','P',p,'H',h,fluid);
        mu = get_fluid_props('VISCOSITY','T',T,'P',p,fluid);
        k  = get_fluid_props('CONDUCTIVITY','T',T,'P',p,fluid);
        Pr = get_fluid_props('PRANDTL','T',T,'P',p,fluid);

        G  = mdot / (pi*(D/2)^2);
        Re = G * D / mu;

        if Re < 2300, Nu = 4.36; else, Nu = 0.023*Re^0.8*Pr^0.4; end
        alpha = Nu * k / D;
    end
end

function alpha = calc_alpha_condenser_air(T_air, T_wall, D_con_o, params)
    % 强迫对流（参数取典型值）
    k_air  = 0.0263;
    mu_air = 1.85e-5;
    rho_air= 1.16;
    d1 = 1.5e-3; d2 = 25e-3; N_rows = 1; V_air = 1.5;
    Re = rho_air * V_air * D_con_o / mu_air;

    alpha = 0.982 * Re^0.424 * (d1/D_con_o)^-0.0887 ...
          * (N_rows*d2/D_con_o)^-0.1590 * (k_air/D_con_o);
end

function alpha = calc_alpha_evaporator_air(T_air, T_wall, L_eva, params)
    % 自然对流关联式
    g = 9.81;
    T_f = 0.5*(T_air + T_wall);
    beta = 1/max(T_f, 200); % 防止数值异常
    k_air  = 0.0263;
    mu_air = 1.85e-5;
    rho_air= 1.16;
    Pr = 0.71;
    nu = mu_air / rho_air;

    Gr = (g * beta * abs(T_wall - T_air) * L_eva^3) / (nu^2);
    term1 = 0.75 * Pr^0.5;
    term2 = (0.609 + 1.221*Pr^0.5 + 1.238*Pr)^0.25;

    alpha = (4/3) * (k_air / L_eva) * (max(Gr,1e-6)/4)^0.25 * (term1/term2);
end

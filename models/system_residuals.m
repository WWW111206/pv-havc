function F = system_residuals(x, bcs, params, control_input)
% 闭环残差：质量守恒 + 蒸发器出口过热度目标
% 未知量：x = [p_con ; p_eva] (Pa)

    p_con = x(1);
    p_eva = x(2);

    % 基本边界检查
    if (p_con <= p_eva) || (p_eva <= 1e4) || (p_con > 40e5)
        F = [1e6; 1e6];
        return;
    end

    fluid = 'R134a';

    % 若未给箱内空气温度，默认等于环境
    if ~isfield(bcs, 'T_air_chamber') || isempty(bcs.T_air_chamber)
        bcs.T_air_chamber = bcs.T_amb;
    end

    try
        % 1) 压缩机吸气（初值取饱和蒸汽）
        h_comp_in = get_fluid_props('H','P',p_eva,'Q',1,fluid);

        % 2) 压缩机
        comp_out   = model_compressor(control_input, p_eva, p_con, h_comp_in, params);
        mdot_comp  = comp_out.mdot_ref;
        h_comp_out = comp_out.h_ref_out;

        % 3) 冷凝器
        cond_out   = model_heat_exchanger('condenser', mdot_comp, p_con, h_comp_out, bcs.T_amb, params);
        h_cond_out = cond_out.h_out_ref;

        % 4) 毛细管
        cap_out   = model_capillary(p_con, h_cond_out, p_eva, params);
        mdot_cap  = cap_out.mdot_ref;

        % 5) 蒸发器（用毛细管出口焓作为入口）
        eva_out   = model_heat_exchanger('evaporator', mdot_cap, p_eva, cap_out.h_out, bcs.T_air_chamber, params);
        h_eva_out = eva_out.h_out_ref;

    catch ME
        warning('组件模型失败 @ (p_con=%.2fbar, p_eva=%.2fbar): %s', p_con/1e5, p_eva/1e5, ME.message);
        F = [1e6; 1e6];
        return;
    end

    % 残差1：质量守恒
    F1 = mdot_comp - mdot_cap;
    if abs(mdot_comp) > 1e-8
        F1 = F1 / mdot_comp; % 归一化
    end

    % 残差2：蒸发器出口过热度（目标 superheat_K）
    try
        T_sat_eva = get_fluid_props('T','P',p_eva,'Q',1,fluid);
        T_target  = T_sat_eva + params.evaporator.superheat_K;
        h_target  = get_fluid_props('H','P',p_eva,'T',T_target,fluid);
    catch
        h_target = h_comp_in; % 兜底
    end
    F2 = (h_eva_out - h_target) / max(abs(h_target), 1);

    F = [F1; F2];
end

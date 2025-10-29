function props = get_fluid_props(query, param1_char, param1_val, param2_char, param2_val, fluid)
% 封装CoolProp的PropsSI函数
% 示例: h = get_fluid_props('H', 'P', 101325, 'T', 300, 'R134a');
    try
        % 调用Python中的CoolProp库
        props = py.CoolProp.CoolProp.PropsSI(query, param1_char, param1_val, param2_char, param2_val, fluid);
    catch ME
        error('CoolProp调用失败: %s', ME.message);
    end
end

function val = get_fluid_props(query, p1_char, p1_val, p2_char, p2_val, fluid)
% 封装 CoolProp 的 PropsSI，支持常数性质与短名称映射，并返回 double
% 用法示例：
%   h = get_fluid_props('H', 'P', 101325, 'T', 300, 'R134a');
%   Pcrit = get_fluid_props('Pcrit', '', 0, '', 0, 'R134a'); % 常数性质
%
% 注意：本函数依赖 MATLAB 的 Python 接口，请先配置 Python 与 CoolProp:
%   pyenv;  以及 pip install CoolProp

    % 映射简写到 CoolProp 标准键
    mapKey = @(k) local_map_key(k);

    q = mapKey(query);
    k1 = mapKey(p1_char);
    k2 = mapKey(p2_char);

    try
        % 情形一：请求常数性质（如 Pcrit、Tcrit 等），允许 p1_char/p2_char 为空
        if (isempty(k1) || strcmp(k1,'')) && (isempty(k2) || strcmp(k2,''))
            % Python 接口支持 PropsSI('Pcrit', 'R134a')
            out = py.CoolProp.CoolProp.PropsSI(q, fluid);
        else
            % 常规模式：6 参数
            out = py.CoolProp.CoolProp.PropsSI(q, k1, double(p1_val), k2, double(p2_val), fluid);
        end
        % 转 double
        try
            val = double(out);
        catch
            % 兜底：尝试将 Python 数值转为 float 再转 double
            val = double(py.float(out));
        end
    catch ME
        error('CoolProp调用失败(%s): %s', q, ME.message);
    end
end

function key = local_map_key(keyIn)
    if nargin==0 || isempty(keyIn)
        key = '';
        return;
    end
    k = upper(string(keyIn));

    % 输出/输入属性键映射（常用简写到标准键）
    switch k
        % 基本性质
        case {'T','P','H','S','D','Q'}
            key = char(k);

        % 常数性质
        case {'PCRIT','TCRIT','RHOLIMIT','OMEGA','MM','ACF','P_TRIPLE','T_TRIPLE'}
            key = char(k);

        % 比热
        case {'C','CP','CPMASS'}
            key = 'Cpmass';
        case {'O','CVMASS'}
            key = 'Cvmass';

        % 黏度、导热、普朗特
        case {'V','MU','VISCOSITY'}
            key = 'VISCOSITY';
        case {'L','K','CONDUCTIVITY'}
            key = 'CONDUCTIVITY';
        case {'PR','PRANDTL'}
            key = 'PRANDTL';

        % 其他别名（可按需扩展）
        otherwise
            key = char(keyIn); % 直接透传，CoolProp 会校验
    end
end

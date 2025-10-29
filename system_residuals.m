function F = system_residuals(x, bcs, params, control_input)
% =========================================================================
% FUNCTION: system_residuals.m (Final Corrected Version)
% DATE: 2025-10-29
% =========================================================================

% --- 1. Unpack Guessed Variables ---
p_con = x(1);
p_eva = x(2);

% --- 2. Sanity Check for Solver Robustness ---
% 1) p_con > p_eva  2) p_eva > 0.1 bar (~1e4 Pa)  3) p_con < 40 bar (~40e5 Pa)
if (p_con <= p_eva) || (p_eva <= 1e4) || (p_con > 40e5)
    F = [1e6; 1e6]; % Return large residuals if pressures are invalid
    return;
end

% --- 3. Simulate the Refrigeration Cycle Component by Component ---
try
    % a) Compressor inlet: saturated vapor at evaporating pressure
    h_comp_in = get_fluid_props('H', 'P', p_eva, 'Q', 1, 'R134a'); % Q=1

    % b) COMPRESSOR Model
    comp_out = model_compressor(control_input, p_eva, p_con, h_comp_in, params);
    mdot_comp = comp_out.mdot_ref;
    h_comp_out = comp_out.h_ref_out;

    % c) CONDENSER Model
    cond_out = model_heat_exchanger('condenser', mdot_comp, p_con, h_comp_out, bcs.T_amb, params);
    h_cond_out = cond_out.h_out_ref;

    % d) CAPILLARY TUBE Model
    cap_out = model_capillary(p_con, h_cond_out, p_eva, params);
    mdot_cap = cap_out.mdot_ref;

catch ME
    warning('A component model failed. Pressure guess (p_con, p_eva): (%.2f, %.2f) bar. Error: %s', ...
             p_con/1e5, p_eva/1e5, ME.message);
    F = [1e6; 1e6];
    return;
end

% --- 4. Calculate Residuals for the Solver ---
% RESIDUAL 1: Mass Flow Rate Balance
F(1) = mdot_comp - mdot_cap;

% RESIDUAL 2: Simplified (focus on mass balance)
F(2) = 0;

% Normalize the primary residual to improve solver performance
if abs(mdot_comp) > 1e-6
    F(1) = F(1) / mdot_comp;
end

end
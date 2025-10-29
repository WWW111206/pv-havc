function cap_out = model_capillary(p_con, h_in, p_eva, params) %#ok<INUSD>
% =========================================================================
% FUNCTION: model_capillary.m
% AUTHOR:
% DATE: 2025-10-29
%
% DESCRIPTION:
% This function models the refrigerant flow through the capillary tube.
% It calculates the mass flow rate based on the empirical correlation from
% the ASHRAE handbook (Eq. 19 in the source paper) and assumes an
% isenthalpic (constant enthalpy) expansion process (Eq. 20).
%
% NOTE ON UNUSED PARAMETER 'p_eva':
% The MATLAB editor will flag 'p_eva' as an unused input. This is expected
% because the specific empirical correlation (Eq. 19) does not depend on
% the outlet pressure. However, p_eva is a primary variable in the overall
% system solver, so it is passed to this function to maintain a consistent
% interface across all component models. The %#ok<INUSD> tag is used to
% suppress this expected warning.
%
% INPUTS:
%   p_con      - Condensing pressure at capillary inlet (Pa).
%   h_in       - Specific enthalpy of refrigerant at inlet (J/kg).
%   p_eva      - Evaporating pressure at capillary outlet (Pa) [Unused].
%   params     - A struct containing all system and geometric parameters.
%
% OUTPUTS:
%   cap_out    - A struct containing the capillary's performance:
%   .mdot_ref - Refrigerant mass flow rate (kg/s).
%   .h_out    - Specific enthalpy of refrigerant at outlet (J/kg).
% =========================================================================

% --- 1. Unpack Parameters ---
D_cap = params.capillary.D_in; % Capillary inside diameter (m) 
L_cap = params.capillary.L;    % Capillary length (m) 
fluid = 'R134a';

% --- 2. Calculate Degree of Subcooling (ΔT) ---
% Subcooling is the difference between the saturation temperature at the
% condenser pressure and the actual temperature of the refrigerant entering
% the capillary.

try
    % Get saturation temperature corresponding to the condensing pressure
    T_sat_at_pcon = get_fluid_props('T', 'P', p_con, 'Q', 0, fluid); % Q=0 for saturated liquid
    
    % Get the actual temperature of the refrigerant at the inlet
    T_in_actual = get_fluid_props('T', 'P', p_con, 'H', h_in, fluid);
catch ME
    error('CoolProp failed in capillary model: %s', ME.message);
end

% Calculate the degree of subcooling
delta_T = T_sat_at_pcon - T_in_actual;

% Ensure subcooling is not negative. If it is, the fluid is not subcooled,
% and the correlation might not be valid. We set it to a small positive
% number to avoid math errors with non-integer exponents.
if delta_T <= 0
    delta_T = 1e-6;
end

% --- 3. Calculate Mass Flow Rate (Eq. 19) ---
% The formula requires pressure in bar, so we convert from Pa.
p_con_bar = p_con * 1e-5;

% ASHRAE empirical correlation from the paper (Eq. 19) 
mdot_ref = 0.006369 * (p_con_bar^0.916) * (delta_T^0.632) * (D_cap^1.832) * (L_cap^-0.416);

% --- 4. Apply Isenthalpic Expansion (Eq. 20) ---
% The expansion process in a capillary tube is assumed to be adiabatic
% with no work done, so the enthalpy remains constant.
h_out = h_in;

% --- 5. Populate Output Structure ---
cap_out.mdot_ref = mdot_ref;
cap_out.h_out = h_out;

end
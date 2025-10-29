function hx_out = model_heat_exchanger(type, mdot_ref, p_ref, h_in_ref, T_air_in, params)
% =========================================================================
% FUNCTION: model_heat_exchanger.m
% AUTHOR:
% DATE: 2025-10-29
%
% DESCRIPTION:
% This function models a heat exchanger (condenser or evaporator) using a
% 1D steady-state finite volume method. It discretizes the heat exchanger
% into a number of control volumes and solves the energy balance for each
% volume iteratively along the length of the tube.
%
% The function is generic and can be configured to model either the
% condenser or the evaporator by specifying the 'type' argument. It uses
% the specific heat transfer correlations from the source paper.
%
% METHOD:
% 1. Discretize the heat exchanger tube into 'N' segments.
% 2. March from the inlet (segment 1) to the outlet (segment N).
% 3. For each segment, iteratively solve for the wall temperature and
%    refrigerant outlet enthalpy by balancing the heat transfer from the
%    refrigerant to the wall, and from the wall to the ambient air.
% 4. Use specific correlations for refrigerant-side and air-side heat
%    transfer coefficients based on the heat exchanger type and flow regime.
%
% INPUTS:
%   type       - String: 'condenser' or 'evaporator'.
%   mdot_ref   - Refrigerant mass flow rate (kg/s).
%   p_ref      - Refrigerant pressure (Pa), assumed constant.
%   h_in_ref   - Specific enthalpy of refrigerant at inlet (J/kg).
%   T_air_in   - Ambient air temperature (K).
%   params     - A struct containing all system and geometric parameters.
%
% OUTPUTS:
%   hx_out     - A struct containing the heat exchanger's performance:
%   .h_out_ref - Specific enthalpy of refrigerant at outlet (J/kg).
%   .Q_total   - Total heat transfer rate (W).
%   .T_wall_avg- Average wall temperature (K).
% =========================================================================

% --- 1. Initialization and Parameter Unpacking ---
fluid = 'R134a';
N = 50; % Number of control volumes for discretization

% Unpack parameters based on the heat exchanger type
if strcmpi(type, 'condenser')
    L = params.condenser.L;         % Total length (m) 
    D_in = params.condenser.D_in;   % Inner diameter (m) 
    D_out = params.condenser.D_out; % Outer diameter (m) 
else % 'evaporator'
    L = params.evaporator.L;        % Total length (m) 
    D_in = params.evaporator.D_in;  % Inner diameter (m) 
    D_out = params.evaporator.D_out;% Outer diameter (m) 
end

dL = L / N; % Length of one control volume (m)
A_in = pi * D_in * dL; % Inner surface area of one segment (m^2)
A_out = pi * D_out * dL; % Outer surface area of one segment (m^2)

% Handle zero mass flow case to prevent errors
if mdot_ref < 1e-9
    hx_out.h_out_ref = h_in_ref;
    hx_out.Q_total = 0;
    hx_out.T_wall_avg = T_air_in;
    return;
end

% Initialize state vectors
h_ref = zeros(N+1, 1); % Enthalpy at the boundaries of each segment
T_wall = zeros(N, 1);  % Average wall temperature of each segment

h_ref(1) = h_in_ref;   % Set inlet condition

% --- 2. Spatial Marching Loop (Finite Volume Method) ---
for i = 1:N
    % Use outlet enthalpy of previous segment as initial guess for current
    h_out_guess = h_ref(i);
    
    % Inner loop to solve for the state of the current segment (i)
    for iter = 1:20 % Max 20 iterations for convergence
        
        % Average enthalpy in the current control volume
        h_avg = (h_ref(i) + h_out_guess) / 2;
        
        % --- Get Refrigerant Properties ---
        try
            T_ref_avg = get_fluid_props('T', 'P', p_ref, 'H', h_avg, fluid);
            x_avg = get_fluid_props('Q', 'P', p_ref, 'H', h_avg, fluid);
        catch ME
            warning('CoolProp failed in heat exchanger model at segment %d. Using previous segment properties. Error: %s', i, ME.message);
            h_ref(i+1) = h_ref(i);
            T_wall(i) = T_wall(i-1);
            continue; % Skip to next segment
        end
        
        % --- Calculate Heat Transfer Coefficients (alpha) ---
        if strcmpi(type, 'condenser')
            alpha_ref = calc_alpha_condenser_ref(p_ref, h_avg, mdot_ref, D_in, params);
            alpha_air = calc_alpha_condenser_air(T_air_in, T_ref_avg, D_out, params);
        else % 'evaporator'
            % For evaporator, heat flux 'q' is needed for alpha_ref, which depends on T_wall.
            % We first estimate T_wall to calculate q, then iterate.
            T_wall_guess = (T_ref_avg + T_air_in) / 2;
            q_guess = calc_alpha_evaporator_air(T_air_in, T_wall_guess, L, params) * (T_wall_guess - T_air_in);
            alpha_ref = calc_alpha_evaporator_ref(p_ref, h_avg, mdot_ref, D_in, q_guess, params);
            alpha_air = calc_alpha_evaporator_air(T_air_in, T_wall_guess, L, params);
        end
        
        % --- Solve for Wall Temperature and Heat Transfer ---
        % From energy balance at the wall: alpha_ref*A_in*(T_ref-T_wall) = alpha_air*A_out*(T_wall-T_air)
        % This assumes negligible wall resistance, consistent with the paper's formulation.
        T_wall_new = (alpha_ref * A_in * T_ref_avg + alpha_air * A_out * T_air_in) / (alpha_ref * A_in + alpha_air * A_out);
        
        % Total heat transfer in this segment (from air to wall)
        dQ = alpha_air * A_out * (T_air_in - T_wall_new);
        
        % Update refrigerant outlet enthalpy
        % For condenser, dQ is negative (heat rejected). For evaporator, dQ is positive (heat absorbed).
        h_out_new = h_ref(i) - dQ / mdot_ref;
        
        % Check for convergence of the inner loop
        if abs(h_out_new - h_out_guess) < 1 % Tolerance of 1 J/kg
            break;
        end
        
        h_out_guess = h_out_new; % Update guess for next iteration
    end
    
    % Store the converged values for the current segment
    h_ref(i+1) = h_out_new;
    T_wall(i) = T_wall_new;
end

% --- 3. Populate Output Structure ---
hx_out.h_out_ref = h_ref(N+1);
hx_out.Q_total = mdot_ref * (h_ref(N+1) - h_ref(1));
hx_out.T_wall_avg = mean(T_wall);

end

% =========================================================================
%                        --- HELPER FUNCTIONS ---
% =========================================================================

function alpha = calc_alpha_condenser_ref(p, h, mdot, D, params)
    % Implements refrigerant-side heat transfer for condenser
    % Uses Shah correlation for two-phase and Dittus-Boelter for single-phase.
    
    fluid = 'R134a';
    x = get_fluid_props('Q', 'P', p, 'H', h, fluid);
    
    if x >= 0 && x <= 1 % Two-phase region (Eq. 25)
        p_crit = get_fluid_props('Pcrit', '', 0, '', 0, fluid);
        p_r = p / p_crit; % Reduced pressure
        
        k_L = get_fluid_props('L', 'P', p, 'Q', 0, fluid); % Thermal conductivity of liquid
        mu_L = get_fluid_props('V', 'P', p, 'Q', 0, fluid); % Viscosity of liquid
        Pr_L = get_fluid_props('Prandtl', 'P', p, 'Q', 0, fluid); % Prandtl number of liquid
        
        G = mdot / (pi*(D/2)^2); % Mass flux
        Re_L = G * D / mu_L; % Reynolds number if all flow is liquid
        
        % Shah correlation (1979) as cited in paper (Eq. 25)
        alpha_i = 0.023 * Re_L^0.8 * Pr_L^0.4 * (k_L / D);
        alpha = alpha_i * ((1-x)^0.8 + (3.8 * x^0.76 * (1-x)^0.04) / (p_r^0.38));
        
    else % Single-phase region (liquid or vapor)
        % Using standard Dittus-Boelter correlation (a reasonable assumption)
        T = get_fluid_props('T', 'P', p, 'H', h, fluid);
        mu = get_fluid_props('V', 'T', T, 'P', p, fluid);
        k = get_fluid_props('L', 'T', T, 'P', p, fluid);
        Pr = get_fluid_props('Prandtl', 'T', T, 'P', p, fluid);
        
        G = mdot / (pi*(D/2)^2);
        Re = G * D / mu;
        
        if Re < 2300, Nu = 3.66; else, Nu = 0.023 * Re^0.8 * Pr^0.3; end % n=0.3 for cooling
        alpha = Nu * k / D;
    end
end

function alpha = calc_alpha_evaporator_ref(p, h, mdot, D, q, params)
    % Implements refrigerant-side heat transfer for evaporator
    % Uses Wang and Touber correlation (Eq. 26-30)
    
    fluid = 'R134a';
    x = get_fluid_props('Q', 'P', p, 'H', h, fluid);
    
    if x > 0 && x < 1 % Two-phase region
        % Get liquid and vapor properties
        rho_L = get_fluid_props('D', 'P', p, 'Q', 0, fluid);
        rho_G = get_fluid_props('D', 'P', p, 'Q', 1, fluid);
        mu_L = get_fluid_props('V', 'P', p, 'Q', 0, fluid);
        mu_G = get_fluid_props('V', 'P', p, 'Q', 1, fluid);
        h_L = get_fluid_props('H', 'P', p, 'Q', 0, fluid);
        h_G = get_fluid_props('H', 'P', p, 'Q', 1, fluid);
        Lat_ref = h_G - h_L; % Latent heat
        
        % Boiling number (Eq. 28)
        Bo = q / (mdot * Lat_ref);
        
        % Lockhart-Martinelli parameter (Eq. 29)
        Xtt = (x / (1-x))^0.9 * (rho_G / rho_L)^0.5 * (mu_L / mu_G)^0.1;
        
        % Coefficients a1 and a2 (Eq. 27, 30)
        a1 = 1 + 24000 * Bo^1.16 + 1.37 * Xtt^-0.86;
        G = mdot / (pi*(D/2)^2);
        Re_L = G * D / mu_L;
        a2 = (1 + 1.15e6 * a1^2 * Re_L^1.17)^-1;
        
        % Heat transfer coefficient (Eq. 26)
        M = 102.03; % Molar mass of R134a (g/mol)
        alpha = a1 * a2 * (p/1e5)^0.24 * M^-0.5 * q^0.7;
        
    else % Single-phase region (liquid or vapor)
        % Using standard Dittus-Boelter correlation
        T = get_fluid_props('T', 'P', p, 'H', h, fluid);
        mu = get_fluid_props('V', 'T', T, 'P', p, fluid);
        k = get_fluid_props('L', 'T', T, 'P', p, fluid);
        Pr = get_fluid_props('Prandtl', 'T', T, 'P', p, fluid);
        
        G = mdot / (pi*(D/2)^2);
        Re = G * D / mu;
        
        if Re < 2300, Nu = 4.36; else, Nu = 0.023 * Re^0.8 * Pr^0.4; end % n=0.4 for heating
        alpha = Nu * k / D;
    end
end

function alpha = calc_alpha_condenser_air(T_air, T_wall, D_con_o, params)
    % Implements air-side heat transfer for condenser (forced convection)
    % Eq. 31
    
    % Assumed air properties at film temperature
    T_film = (T_air + T_wall) / 2;
    k_air = 0.0263; % W/m.K at approx 300K
    mu_air = 1.85e-5; % Pa.s
    rho_air = 1.16; % kg/m^3
    
    % Assumed parameters from typical fin-tube geometry, as they are not in Table 1
    d1 = 1.5e-3; % fin spacing (m)
    d2 = 25e-3;  % pipe spacing (m)
    N_rows = 1;  % number of pipe rows
    V_air = 1.5; % Assumed air velocity (m/s) for fan
    
    Re_air = rho_air * V_air * D_con_o / mu_air;
    
    alpha = 0.982 * Re_air^0.424 * (d1/D_con_o)^-0.0887 * (N_rows*d2/D_con_o)^-0.1590 * (k_air/D_con_o);
end

function alpha = calc_alpha_evaporator_air(T_air, T_wall, L_eva, params)
    % Implements air-side heat transfer for evaporator (natural convection)
    % Eq. 32-33
    
    g = 9.81; % m/s^2
    
    % Air properties at film temperature
    T_film = (T_air + T_wall) / 2;
    beta_air = 1 / T_film; % Thermal expansion coefficient for ideal gas (1/K)
    k_air = 0.0263; % W/m.K at approx 300K
    mu_air = 1.85e-5; % Pa.s
    rho_air = 1.16; % kg/m^3
    Pr_air = 0.71;
    nu_air = mu_air / rho_air; % kinematic viscosity
    
    % Grashof number (Eq. 33)
    Gr_eva = (g * beta_air * abs(T_wall - T_air) * L_eva^3) / (nu_air^2);
    
    % Heat transfer coefficient (Eq. 32)
    term1 = 0.75 * Pr_air^0.5;
    term2 = (0.609 + 1.221*Pr_air^0.5 + 1.238*Pr_air)^0.25;
    
    alpha = (4/3) * (k_air / L_eva) * (Gr_eva / 4)^0.25 * (term1 / term2);
end

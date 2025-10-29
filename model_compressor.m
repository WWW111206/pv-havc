function comp_out = model_compressor(control_input, p_eva, p_con, h_ref_in, params)
% =========================================================================
% FUNCTION: model_compressor.m
% AUTHOR:
% DATE: 2025-10-29
%
% DESCRIPTION:
% This function models the behavior of the variable-speed DC compressor
% based on the equations (9) through (17) from the source paper. It can
% operate in two modes:
%   1. 'CSP' (Compressor Speed Prediction): Calculates the required
%      electrical power for a given speed.
%   2. 'MPPT' (Maximum Power Point Tracking): Solves for the operating
%      speed that corresponds to a given electrical input power.
%
% INPUTS:
%   control_input - A struct with control parameters:
%    .mode       - String: 'CSP' or 'MPPT'
%    .omega      - Compressor angular velocity (rad/s) for 'CSP' mode
%    .P_in       - Compressor electrical input power (W) for 'MPPT' mode
%   p_eva         - Evaporating pressure (Pa)
%   p_con         - Condensing pressure (Pa)
%   h_ref_in      - Specific enthalpy of refrigerant at compressor inlet (J/kg)
%   params        - A struct containing all system parameters
%
% OUTPUTS:
%   comp_out      - A struct containing the compressor's performance metrics:
%    .omega      - Resulting angular velocity (rad/s)
%    .P_in_elec  - Resulting electrical input power (W)
%    .mdot_ref   - Refrigerant mass flow rate (kg/s)
%    .h_ref_out  - Specific enthalpy of refrigerant at outlet (J/kg)
%    .I_motor    - Motor current (A)
%    .U_motor    - Motor voltage (V)
%    .Tor_load   - Calculated average load torque (N*m)
% =========================================================================

% --- 1. Unpack Parameters and Define Constants ---
R  = params.compressor.R;         % Line resistance (Ohm) 
Ke = params.compressor.Ke;        % Back EMF coefficient (V*s/rad) 
% For a PMDC motor, torque constant Kt (Nm/A) = Ke (Vs/rad) in SI units
Kt = Ke;
piston_D = params.compressor.piston_D; % Piston diameter (m) 
piston_S = params.compressor.piston_S; % Piston stroke (m) 
lambda_d = 0.98;                  % Leakage coefficient (dimensionless) 
fluid = 'R134a';                  % Refrigerant type

% --- 2. Calculate Inlet Properties and Compressor Geometry ---
% Get refrigerant properties at the inlet using the CoolProp wrapper
try
    rho_ref_in = get_fluid_props('D', 'P', p_eva, 'H', h_ref_in, fluid);
    cp_ref_in  = get_fluid_props('C', 'P', p_eva, 'H', h_ref_in, fluid);
    cv_ref_in  = get_fluid_props('O', 'P', p_eva, 'H', h_ref_in, fluid); % 'O' for Cvmass
catch ME
    error('CoolProp failed in compressor model: %s', ME.message);
end

% Adiabatic index (kappa)
if cv_ref_in > 0
    kappa = cp_ref_in / cv_ref_in;
else
    kappa = 1.13; % Fallback value if Cv is not available
end

% Swept volume (V_max in paper)
A_pis = pi * (piston_D/2)^2;
V_swept = A_pis * piston_S;

% --- 3. Calculate Compressor Load Torque ---
% This part calculates the average torque required to compress the refrigerant,
% which is independent of speed in this formulation.

% Volumetric efficiency coefficient due to clearance volume (Eq. 16)
lambda_v = 1 - 0.1 * ((p_con / p_eva)^(1/kappa) - 1);
if lambda_v < 0, lambda_v = 0; end % Cannot be negative

% Total volumetric efficiency (Eq. 15, assuming pressure and temp coeffs are 1)
lambda = lambda_v * lambda_d;

% Indicated work per cycle for an ideal reciprocating compressor. This is
% a standard thermodynamic formula representing the area of the P-V diagram.
work_per_cycle = (kappa / (kappa-1)) * p_eva * V_swept *...
                 ((p_con/p_eva)^((kappa-1)/kappa) - 1);

% Average load torque is the work per revolution (Eq. 10, 11, 12, 13 averaged)
Tor_load = work_per_cycle / (2*pi);

% --- 4. Solve for Operating Point Based on Control Mode ---
switch upper(control_input.mode)
    case 'CSP'
        % For CSP mode, speed is the input. We calculate the required power.
        omega = control_input.omega;
        
        % Calculate motor current from torque balance (Eq. 10, steady state)
        I_motor = Tor_load / Kt;
        
        % Calculate motor voltage from voltage equation (Eq. 9, steady state)
        U_motor = R * I_motor + Ke * omega;
        
        % Calculate total electrical power input
        P_in_elec = U_motor * I_motor;

    case 'MPPT'
        % For MPPT mode, power is the input. We must solve for the speed.
        P_in_elec = control_input.P_in;
        
        % We need to find omega where the motor's output torque at that speed
        % (given P_in_elec) equals the compressor's load torque.
        % Residual function: Tor_motor(omega) - Tor_load = 0
        residual_fun = @(w) calc_motor_torque(w, P_in_elec, R, Ke) - Tor_load;
        
        % Use fsolve to find the root (omega)
        options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-8);
        omega_guess = 2800 * 2*pi/60; % Initial guess: 2800 rpm, a reasonable mid-point
        
        [omega, ~, exitflag] = fsolve(residual_fun, omega_guess, options);
        
        % Check for solver convergence
        if exitflag <= 0
            warning('Compressor model fsolve did not converge for P_in = %.2f W. Outputting NaN.', P_in_elec);
            omega = NaN;
        end
        
        % With the solved omega, calculate the resulting U and I for output
        if isfinite(omega)
            U_motor = (Ke*omega + sqrt((Ke*omega)^2 + 4*R*P_in_elec))/2;
            I_motor = P_in_elec / U_motor;
        else
            U_motor = NaN;
            I_motor = NaN;
        end

    otherwise
        error('Unknown control mode specified for compressor model: %s', control_input.mode);
end

% --- 5. Calculate Final Outputs (Common for Both Modes) ---
if isfinite(omega) && omega > 0
    % Refrigerant mass flow rate (Eq. 14)
    mdot_ref = lambda * rho_ref_in * V_swept * (omega / (2*pi));
    
    % Outlet enthalpy from energy balance (Eq. 17)
    if mdot_ref > 1e-9 % Avoid division by zero for very low/zero flow
        h_ref_out = h_ref_in + (P_in_elec - I_motor^2 * R) / mdot_ref;
    else
        h_ref_out = h_ref_in; % No flow, no enthalpy change
    end
else
    % If solver failed or speed is zero, set outputs to zero/default
    omega = 0;
    mdot_ref = 0;
    h_ref_out = h_ref_in;
    I_motor = 0;
    U_motor = 0;
end

% --- 6. Populate Output Structure ---
comp_out.omega = omega;
comp_out.P_in_elec = P_in_elec;
comp_out.mdot_ref = mdot_ref;
comp_out.h_ref_out = h_ref_out;
comp_out.I_motor = I_motor;
comp_out.U_motor = U_motor;
comp_out.Tor_load = Tor_load;

end

% --- Helper function to calculate motor torque for fsolve in MPPT mode ---
function Tor_motor = calc_motor_torque(omega, P_in, R, Ke)
    % This function calculates the mechanical torque produced by the motor
    % for a given speed (omega) and electrical input power (P_in).

    if omega < 1e-3 % Avoid division by zero at or near standstill
        % At zero speed, all power is dissipated as heat, no mechanical torque
        Tor_motor = 0; 
        return;
    end
    
    % From P_in = U*I and U = R*I + Ke*omega, we can derive a quadratic
    % equation for U: U^2 - (Ke*omega)*U - R*P_in = 0.
    % Solving for U (we take the positive root):
    U = (Ke*omega + sqrt((Ke*omega)^2 + 4*R*P_in))/2;
    
    % Calculate current
    I = P_in / U;
    
    % Mechanical Power = Electrical Power - Resistive Losses
    P_mech = P_in - I^2*R;
    
    % Torque = Mechanical Power / Angular Velocity
    Tor_motor = P_mech / omega;
end
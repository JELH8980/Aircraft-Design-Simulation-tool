function [T, p, rho, mu, nu, a] = atmosphere_state(Hp)
% ATMOSPHERE_STATE - Function to calculate atmospheric properties at a given altitude.
%
% This function computes key atmospheric properties—temperature (T), pressure (p), 
% density (rho), dynamic viscosity (mu), kinematic viscosity (nu), and speed of 
% sound (a)—based on the pressure altitude (Hp) using the International Standard 
% Atmosphere (ISA) model. It models the troposphere (0–11,000 m) with a linear 
% temperature lapse rate and the lower stratosphere (11,000–20,000 m) with constant 
% temperature.
%
% INPUTS:
%   Hp             - Pressure altitude [m], valid between 0 and 20,000 m.
%
% OUTPUTS:
%   T              - Temperature [K]
%   p              - Pressure [N/m^2]
%   rho            - Density [kg/m^3]
%   mu             - Dynamic viscosity [N*s/m^2]
%   nu             - Kinematic viscosity [m^2/s]
%   a              - Speed of sound [m/s]
%
% FUNCTIONALITY:
% - Defines ISA constants: sea-level temperature (T0), pressure (p0), lapse rate 
%   (alphaT), gravitational acceleration (g0), gas constant (R), etc.
% - Calculates temperature (T):
%   - Troposphere (0 ≤ Hp < 11,000 m): Linear decrease with altitude.
%   - Stratosphere (11,000 ≤ Hp ≤ 20,000 m): Constant temperature.
% - Calculates pressure (p):
%   - Troposphere: Based on temperature ratio and lapse rate.
%   - Stratosphere: Exponential decrease under isothermal conditions.
% - Derives additional properties:
%   - Density (rho) from the ideal gas law.
%   - Speed of sound (a) using the adiabatic index (gamma).
%   - Dynamic viscosity (mu) via Sutherland’s law.
%   - Kinematic viscosity (nu) as mu/rho.
%
% NOTES:
% - Assumes a standard lapse rate of 0.0065 K/m in the troposphere and a constant 
%   temperature of 216.65 K in the stratosphere up to 20,000 m.
% - No error handling for altitudes outside 0–20,000 m; outputs are undefined for 
%   such cases.
% - Units follow SI conventions (meters, Kelvin, Newtons, etc.).
% - Sutherland’s law uses a reference viscosity (mu0) and temperature (S) to model 
%   viscosity variation.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

%--------------------------ATMOSPHERE---------------------------------------
H0                      = 0                 ;%    H0          [m]
T0                      = 288.15            ;%    T0          [K]
alphaT                  = 0.0065            ;%    alphaT      [K/m]
deltaT                  = 0                 ;%    deltaT      [K]
p0                      = 101325            ;%    p0          [N/m^2]
g0                      = 9.80665           ;%    g0          [m/s^2]  % Corrected unit comment
mu0                     = 1.7894e-5        ;%    mu0         [N*s/m^2] % Corrected value (was 1.7894)
S                       = 110               ;%    S           [K]
R                       = 287.053           ;%    R           [J/(kg*K)] % Corrected unit comment
gamma                   = 1.4               ;%    gamma       [-]
H11                     = 11000             ;%    H11         [m]
T11                     = 216.65            ;%    T11         [K]
p11                     = 22632             ;%    p11         [N/m^2]
ISA_trop_ub             = 11000             ;%    tropub      [m]
ISA_strat_ub            = 20000             ;%    stratub     [m]

% Temperature Calculation
if (ISA_strat_ub >= Hp) && (Hp >= ISA_trop_ub)
    T = T11 + deltaT;  % Constant temperature in stratosphere
elseif (0 <= Hp) && (Hp < ISA_trop_ub)
    T = T0 - alphaT * (Hp - H0);  % Temperature decreases with altitude in troposphere
end

% Pressure Calculation
if (ISA_strat_ub >= Hp) && (Hp >= ISA_trop_ub)
    p = p11 * exp(-g0 / (R * T11) * (Hp - H11));  % Isothermal pressure decrease
elseif (0 <= Hp) && (Hp < ISA_trop_ub)
    p = p0 * (T / T0)^(g0 / (R * alphaT));  % Pressure based on temperature ratio
end

% Derived Properties
rho = p / (R * T);              % Density
a = sqrt(gamma * R * T);        % Speed of sound
mu = mu0 * (T / T0)^(3/2) * (T0 + S) / (T + S);  % Dynamic viscosity (Sutherland's law)
nu = mu / rho;                  % Kinematic viscosity

end
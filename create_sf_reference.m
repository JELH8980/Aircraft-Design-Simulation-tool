function reference = create_sf_reference()
% CREATE_SF_REFERENCE - Function to create a reference struct for state-space simulation.
%
% This function constructs a reference struct containing initial conditions, control 
% inputs, and model parameters for a state-space simulation of an aircraft. It 
% retrieves data from the MATLAB base workspace (condition, model, signals, etc.) 
% and computes derived states such as angle of attack (alpha) and sideslip angle 
% (beta). The struct is designed to serve as a baseline reference for dynamic analysis.
%
% INPUTS:
%   None           - Relies on variables in the base workspace: condition, model, 
%                    latitude0, longitude0, signal_in_delta_a, signal_in_delta_e, 
%                    signal_in_delta_r, signal_in_psi, signal_in_theta, autopilot, tsim.
%
% OUTPUT:
%   reference      - A struct containing initial states, control inputs, inertial 
%                    properties, and metadata for simulation.
%
% FUNCTIONALITY:
% - Defines gravitational acceleration (g = 9.80665 m/s^2).
% - Retrieves flight condition (h, u, w, theta, delta_e, Tcom), model parameters 
%   (mass, inertia), and input signals from the base workspace using `evalin`.
% - Initializes states: 
%   - Velocities: u, v=0, w, wd=0.
%   - Attitudes: phi=0, theta, psi=0.
%   - Rates: P=0, Q=0, R=0.
%   - Controls: delta_e, delta_a=0, delta_r=0, Tcom.
% - Computes derived states: alpha (angle of attack), beta (sideslip angle).
% - Stores vectors:
%   - Control inputs (ctrl): [delta_e, delta_a, delta_r, Tcom]'.
%   - Position (ip): [0, 0, -h].
%   - Velocity (iv): [u, v, w].
%   - Orientation (io): [phi, theta, psi].
%   - Angular rates (iav): [P, Q, R].
%   - Full state vector (x): [u, v, w, P, Q, R, phi, theta, psi, 0, 0, -h]'.
% - Includes mass (m), inertia matrix (I), gravity (g), and geographic coordinates.
%
% NOTES:
% - Units: SI (meters, seconds, radians), except Tcom in percentage (%).
% - Assumes workspace variables are predefined and correctly formatted (e.g., 
%   condition as a struct with named fields, model.parameters with m, Ixx, etc.).
% - No error handling for missing or invalid workspace variables; `evalin` may throw 
%   exceptions if variables are undefined.
% - The inertia matrix (I) includes off-diagonal Ixz terms for cross-coupling.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

g = 9.80665;

condition = evalin('base', 'condition');
model     = evalin('base', 'model');
latitude  = evalin('base', 'latitude0');
longitude = evalin('base', 'longitude0');

reference.signal_delta_a    = evalin('base', 'signal_in_delta_a');
reference.signal_delta_e    = evalin('base', 'signal_in_delta_e');
reference.signal_delta_r    = evalin('base', 'signal_in_delta_r');
reference.signal_psi        = evalin('base', 'signal_in_psi');
reference.signal_theta      = evalin('base', 'signal_in_theta');
reference.autopilot         = evalin('base', 'autopilot');

reference.tsim              = evalin('base', 'tsim');

h = condition.("h [m]");
u = condition.("u [m/s]");
v = 0;
w = condition.("w [m/s]");
wd = 0;
phi = 0;
theta = condition.("theta [rad]");
psi = 0;
P = 0;
Q = 0;
R = 0;
delta_e = condition.("delta_e [rad]");
delta_a = 0;
delta_r = 0;
Tcom = condition.("Tcom [%]");


reference.h       = h;
reference.u       = u;
reference.v       = v;
reference.w       = w;
reference.wd      = wd;
reference.phi     = phi;
reference.theta   = theta;
reference.psi     = psi;
reference.P       = P;
reference.Q       = Q;
reference.R       = R;
reference.alpha   = atan2(w, u);
reference.beta    = asin(v/sqrt(u^2 + v^2 + w^2));
reference.delta_e = delta_e;
reference.delta_a = delta_a;
reference.delta_r = delta_r;
reference.Tcom    = Tcom;

reference.ctrl = [delta_e, delta_a, delta_r, Tcom]';


reference.ip   = [0, 0, -reference.h];

reference.iv   = [u, v, w];

reference.io   = [phi, theta, psi];

reference.iav  = [P, Q, R];

reference.x    = [u, v, w, P, Q, R, phi, theta, psi, 0, 0, -reference.h]';

reference.m    = model.parameters.m;

reference.I    = [model.parameters.Ixx,                       0, model.parameters.Ixz;
                                     0,   model.parameters.Iyy,                    0; 
                  model.parameters.Ixz,                      0, model.parameters.Izz];

reference.g    = g;

reference.latitude  = latitude;
reference.longitude = longitude;


end
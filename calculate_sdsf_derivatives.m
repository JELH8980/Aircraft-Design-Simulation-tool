function derivatives = calculate_sdsf_derivatives(model, condition)
% CALCULATE_SDSF_DERIVATIVES - Computes the aerodynamic derivatives for 
% an aircraft based on the given flight conditions and propulsion model.
%
% INPUT:
%   model    - Aircraft model containing aerodynamic and propulsion parameters.
%   condition - Flight conditions including Mach number (M), altitude (h), 
%              velocity (u), vertical speed (w), pitch angle (theta), 
%              elevator deflection (delta_e), thrust (T), and throttle 
%              command (Tcom).
%
% OUTPUT:
%   derivatives - Table of aerodynamic derivatives, including force and 
%                 moment derivatives with respect to control surface deflections 
%                 and other state variables (e.g., speed, throttle).
%
% PROCESS:
% - Flight conditions are used to calculate aerodynamic forces and moments.
% - The Vortex Lattice Method (VLM) is employed to obtain aerodynamic coefficients.
% - Derivatives with respect to control surface deflections (elevator, aileron, 
%   rudder) are computed.
% - Throttle sensitivity is evaluated by calculating the effect of changes in 
%   throttle position on the aircraft's performance.
%
% VARIABLES:
%   alpha, beta    - Angle of attack and sideslip angle.
%   state          - Struct containing flight conditions (e.g., altitude, speed, 
%                   control surface deflections).
%   results        - Struct containing computed aerodynamic coefficients.
%   lattice        - Vortex lattice setup for aerodynamic calculations.
%   Tdev           - Thrust at a perturbed Mach number.
%   CT_M, CT_delta_t - Thrust derivatives with respect to Mach number and throttle.
%
% AERODYNAMIC DERIVATIVES:
%   CX, CZ, Cm    - Force and moment coefficients (x, z forces, pitching moment).
%   Cl, Cn        - Roll and yaw moment coefficients.
%   CY            - Sideforce coefficient.
%   derivatives   - Includes derivatives like CX_u (speed), CX_a (alpha), 
%                   CX_delta_e (elevator deflection), etc.
%
% METHOD:
% - The function evaluates various sensitivity derivatives using finite differences 
%   and vortex lattice method to compute aerodynamic properties at different 
%   perturbation levels (e.g., Mach number, control deflections).

%% Acquire Straight Flight Conditions
addpath("aircraft\propulsion\")
% Create Handle for Thrust Computation
    
engine_model = str2func(strrep(model.parameters.propulsion.method, '.m', ''));  % Convert function name to handle

% Perform Loop to get Derivatives

Tmax = model.parameters.propulsion.Tmax;


M          = condition.M;
h          = condition.("h [m]");
u          = condition.("u [m/s]");
w          = condition.("w [m/s]");
theta      = condition.("theta [rad]");
delta_e    = condition.("delta_e [rad]");
T          = condition.("T [N]");
Tcom       = condition.("Tcom [%]");

[~, ~, rho, ~, ~, a] = atmosphere_state(h);

alpha   = atan2(w, u);

% Format
state.alpha     = alpha;
state.betha     = 0;
state.P         = 0;
state.Q         = 0;
state.R         = 0;
state.alphadot  = 0;
state.bethadot  = 0;

delta_a         = 0;
delta_e         = delta_e;
delta_r         = 0;


% Setting Control Surface Deflections
model = set_flaps(model, "delta_a", delta_a);
model = set_flaps(model, "delta_e", delta_e);
model = set_flaps(model, "delta_r", delta_r);

state.ALT  = h;
state.AS   = M*a;
state.rho  = rho;


if  M > 0.3
state.pgcorr = 1;                           % Correction for viscous effects (Prandtl-Gaertner correction factor)
else
state.pgcorr = 0;    
end

[delta_a_indx, ~] = find(model.parameters.delta_a.lmatrix == 1, 1, 'first');
[delta_e_indx, ~] = find(model.parameters.delta_e.lmatrix == 1, 1, 'first');
[delta_r_indx, ~] = find(model.parameters.delta_r.lmatrix == 1, 1, 'first');

delta_indices = [delta_a_indx, delta_e_indx, delta_r_indx];

[~, indices] = sort(delta_indices);

delta_a_indx = indices(delta_indices == delta_a_indx);
delta_e_indx = indices(delta_indices == delta_e_indx);
delta_r_indx = indices(delta_indices == delta_r_indx);

% Performing Vortex Lattice Calculations
[lattice, ref] = fLattice_setup(model.geo, state, 0);
results.dwcond=0;           %computation result memory structure.    
[results] = solver(results, state, model.geo, lattice, ref);
[results] = coeff_create(results, lattice, state, ref, model.geo);

A=funsteady(lattice,model.geo,ref,M);

results.CZ_a_dot=-A(1);
results.Cm_a_dot=-A(2);

%Flipping the lattice to get the sideslipe solved.
lattice.COLLOC=[lattice.COLLOC(:,1) lattice.COLLOC(:,3) lattice.COLLOC(:,2)];
lattice.N=[lattice.N(:,1) lattice.N(:,3) lattice.N(:,2)];

lattice.VORTEX2(:,:,1)=lattice.VORTEX(:,:,1);
lattice.VORTEX2(:,:,3)=lattice.VORTEX(:,:,2);
lattice.VORTEX2(:,:,2)=lattice.VORTEX(:,:,3);
lattice.VORTEX=lattice.VORTEX2;

lattice.XYZ2(:,:,1)=lattice.XYZ(:,:,1);
lattice.XYZ2(:,:,3)=lattice.XYZ(:,:,2);
lattice.XYZ2(:,:,2)=lattice.XYZ(:,:,3);
lattice.XYZ=lattice.XYZ2;

A=funsteady(lattice,model.geo,ref,M);

results.CY_b_dot=-A(1);
results.Cn_b_dot=-A(2);

factor1=results.CL/results.CZ;  %L scaling factor  
if isnan(factor1)
factor1=1;
end

factor2=results.CD/results.CZ;  %D scaling factor  
if isnan(factor2)
factor2=0;
end

results.CL_a_dot=results.CZ_a_dot*factor1;
results.CD_a_dot=results.CZ_a_dot*factor2;
results.CC_b_dot=results.CY_b_dot;

CL        = results.CL;
CD        = results.CD;

CX        = -CD*cos(alpha) + CL*sin(alpha);
CZ        = -CD*sin(alpha) - CL*cos(alpha);

Cm        = results.Cm;

CL_a      = results.CL_a;
CD_a      = results.CD_a;



CX_a      = -CD_a*cos(alpha) + CL_a*sin(alpha);
CZ_a      = -CD_a*sin(alpha) - CL_a*cos(alpha);


CD_a_dot  = results.CD_a_dot;
CL_a_dot  = results.CL_a_dot;

CZ_a_dot  = -CD_a_dot*sin(alpha) - CL_a_dot*cos(alpha);

Cm_a      = results.Cm_a;

Cm_a_dot  = results.Cm_a_dot;

CY_b      = results.CC_b;
Cl_b      = results.Cl_b;
Cn_b      = results.Cn_b;

CY_P      = results.CC_P;
Cl_P      = results.Cl_P;
Cn_P      = results.Cn_P;

CL_Q      = results.CL_Q;
CD_Q      = results.CD_Q;

CZ_Q      = -CD_Q*sin(alpha) - CL_Q*cos(alpha);

Cm_Q      = results.Cm_Q;

CY_R      = results.CC_R;
Cl_R      = results.Cl_R;
Cn_R      = results.Cn_R;


Cl_delta_a = results.Cl_d(delta_a_indx);
Cn_delta_a = results.Cn_d(delta_a_indx);

CL_delta_e = results.CL_d(delta_e_indx);
CD_delta_e = results.CD_d(delta_e_indx);

CX_delta_e = -CD_delta_e*cos(alpha) + CL_delta_e*sin(alpha);
CZ_delta_e = -CD_delta_e*sin(alpha) - CL_delta_e*cos(alpha);

Cm_delta_e = results.Cm_d(delta_e_indx);

CY_delta_r = results.CC_d(delta_r_indx);
Cl_delta_r = results.Cl_d(delta_r_indx);
Cn_delta_r = results.Cn_d(delta_r_indx);

% Speeed dependence

rel = 0.01;

Mdev = (1+rel)*M;
state.AS = Mdev * a;

% Performing Vortex Lattice Calculations
[lattice, ref] = fLattice_setup(model.geo, state, 0);   
[results] = solver(results, state, model.geo, lattice, ref);
[results] = coeff_create(results, lattice, state, ref, model.geo);

CLdev = results.CL;
CDdev = results.CD;
Cmdev = results.Cm;

CL_u        = (CLdev - CL)/(Mdev*a - M*a);
CD_u        = (CDdev - CD)/(Mdev*a - M*a);


CX_u        = -CD_u*cos(alpha) + CL_u*sin(alpha);
CZ_u        = -CD_u*sin(alpha) - CL_u*cos(alpha);

Cm_u        = (Cmdev - Cm)/(Mdev*a - M*a);

Tdev = engine_model(Tmax, Tcom, Mdev, h);

CT_M = (Tdev-T)/(Mdev - M);

% Throttle derivatives

Tcomdev = (1+rel)*Tcom;

Tdev = engine_model(Tmax, Tcomdev, M, h);

CT_delta_t = (Tdev-T)/(Tcomdev - Tcom);


derivatives = array2table([CX, CX_u, CX_a, CX_delta_e, ...
                           CZ, CZ_u, CZ_a, CZ_a_dot, CZ_Q, CZ_delta_e, ...
                           Cm, Cm_u, Cm_a, Cm_a_dot, Cm_Q, Cm_delta_e, ...
                           CY_b, CY_P, CY_R, CY_delta_r, ...
                           Cl_b, Cl_P, Cl_R, Cl_delta_a, Cl_delta_r, ...
                           Cn_b, Cn_P, Cn_R, Cn_delta_r, Cn_delta_a, ...
                           CT_M, CT_delta_t], "VariableNames", {'CX', 'CX_u', 'CX_a', 'CX_delta_e', ...
                                                                'CZ', 'CZ_u', 'CZ_a', 'CZ_a_dot', 'CZ_Q', 'CZ_delta_e', ...
                                                                'Cm', 'Cm_u', 'Cm_a', 'Cm_a_dot', 'Cm_Q', 'Cm_delta_e', ...
                                                                'CY_b', 'CY_P', 'CY_R', 'CY_delta_r', ...
                                                                'Cl_b', 'Cl_P', 'Cl_R', 'Cl_delta_a', 'Cl_delta_r', ...
                                                                'Cn_b', 'Cn_P', 'Cn_R', 'Cn_delta_r', 'Cn_delta_a', ...
                                                                'CT_M', 'CT_delta_t'});

end
function derivatives = calculate_sf_derivatives(model, condition)
% CALCULATE_SF_DERIVATIVES - Computes the aerodynamic derivatives for 
% an aircraft based on the given flight conditions using the Vortex Lattice 
% Method (VLM). The function calculates various aerodynamic coefficients 
% and their derivatives with respect to control surface deflections (aileron, 
% elevator, rudder), speed, and other flight parameters.
%
% INPUT:
%   model    - Aircraft model containing geometric, aerodynamic, and propulsion parameters.
%   condition - Flight conditions including Mach number (M), altitude (h), 
%              horizontal speed (u), vertical speed (w), and control surface deflections.
%
% OUTPUT:
%   derivatives - Table containing aerodynamic coefficients and their derivatives 
%                 (e.g., lift, drag, sideforce, moments, and control surface derivatives).
%
% VARIABLES:
%   alpha, beta     - Angle of attack and sideslip angle.
%   state           - Structure holding flight conditions (e.g., altitude, speed, control surface deflections).
%   lattice         - Vortex lattice setup for aerodynamic calculations.
%   results         - Struct containing computed aerodynamic coefficients and their derivatives.
%   derivatives     - Aerodynamic derivatives including control surface effects on forces and moments.
%
% PROCESS:
% - The function calculates the aerodynamic coefficients (CL, CD, Cm, etc.) 
%   for the current flight condition using the vortex lattice method.
% - It also computes the derivatives of these coefficients with respect to 
%   control surface deflections (delta_a, delta_e, delta_r), Mach number, 
%   and other parameters like speed and throttle.
% - The derivatives are output as a table, which provides insights into 
%   how each aerodynamic parameter changes with variations in the flight 
%   condition or control surfaces.
%
% AERODYNAMIC DERIVATIVES:
%   CD, CC, CL    - Aerodynamic coefficients for drag, sideforce, and lift.
%   Cl, Cn        - Moment coefficients for roll and yaw.
%   CD_a, CC_a    - Derivatives with respect to angle of attack.
%   CD_b, CC_b    - Derivatives with respect to sideslip.
%   CD_delta_a, CC_delta_a - Derivatives with respect to aileron deflection.
%   CD_delta_e, CC_delta_e - Derivatives with respect to elevator deflection.
%   CD_delta_r, CC_delta_r - Derivatives with respect to rudder deflection.
%   Derivatives are also computed for dynamic terms like pitch rate (Q), yaw rate (R).



M          = condition.M;
h          = condition.("h [m]");
u          = condition.("u [m/s]");
w          = condition.("w [m/s]");

[~, ~, rho, ~, ~, a] = atmosphere_state(h);

alpha     = atan(w/u);

% Format
state.alpha     = alpha;
state.betha     = 0;
state.P         = 0;
state.Q         = 0;
state.R         = 0;
state.alphadot  = 0;
state.bethadot  = 0;

delta_a         = 0;
delta_e         = condition.("delta_e [rad]");
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

CD        = results.CD;
CC        = results.CC;
CL        = results.CL;
Cl        = results.Cl;
Cm        = results.Cm;
Cn        = results.Cn;

CD_a        = results.CD_a;
CC_a        = results.CC_a;
CL_a        = results.CL_a;
Cl_a        = results.Cl_a;
Cm_a        = results.Cm_a;
Cn_a        = results.Cn_a;

CD_a_dot    = results.CL_a_dot;
CL_a_dot    = results.CL_a_dot;
Cm_a_dot    = results.Cm_a_dot;

CC_b_dot    = results.CC_b_dot;
Cn_b_dot    = results.Cn_b_dot;

CD_b        = results.CD_b;
CC_b        = results.CC_b;
CL_b        = results.CL_b;
Cl_b        = results.Cl_b;
Cm_b        = results.Cm_b;
Cn_b        = results.Cn_b;

CD_P        = results.CD_P;
CC_P        = results.CC_P;
CL_P        = results.CL_P;
Cl_P        = results.Cl_P;
Cm_P        = results.Cm_P;
Cn_P        = results.Cn_P;

CD_Q        = results.CD_Q;
CC_Q        = results.CC_Q;
CL_Q        = results.CL_Q;
Cl_Q        = results.Cl_Q;
Cm_Q        = results.Cm_Q;
Cn_Q        = results.Cn_Q;

CD_R        = results.CD_R;
CC_R        = results.CC_R;
CL_R        = results.CL_R;
Cl_R        = results.Cl_R;
Cm_R        = results.Cm_R;
Cn_R        = results.Cn_R;


CD_delta_a = results.CD_d(delta_a_indx);
CC_delta_a = results.CC_d(delta_a_indx);
CL_delta_a = results.CL_d(delta_a_indx);

CD_delta_e = results.CD_d(delta_e_indx);
CC_delta_e = results.CC_d(delta_e_indx);
CL_delta_e = results.CL_d(delta_e_indx);

CD_delta_r = results.CC_d(delta_r_indx);
CC_delta_r = results.CD_d(delta_r_indx);
CL_delta_r = results.CL_d(delta_r_indx);


Cl_delta_a = results.Cl_d(delta_a_indx);
Cm_delta_a = results.Cm_d(delta_a_indx);
Cn_delta_a = results.Cn_d(delta_a_indx);

Cl_delta_e = results.Cl_d(delta_e_indx);
Cm_delta_e = results.Cm_d(delta_e_indx);
Cn_delta_e = results.Cn_d(delta_e_indx);

Cl_delta_r = results.Cl_d(delta_r_indx);
Cm_delta_r = results.Cm_d(delta_r_indx);
Cn_delta_r = results.Cn_d(delta_r_indx);



derivatives = array2table([CD, CC, CL, Cl, Cm, Cn, ...
                           CD_a, CC_a, CL_a, Cl_a, Cm_a, Cn_a, ...
                           CD_a_dot, CL_a_dot, Cm_a_dot, ...
                           CC_b_dot, Cn_b_dot, ...
                           CD_b, CC_b, CL_b, Cl_b, Cm_b, Cn_b, ...
                           CD_P, CC_P, CL_P, Cl_P, Cm_P, Cn_P, ...
                           CD_Q, CC_Q, CL_Q, Cl_Q, Cm_Q, Cn_Q, ...
                           CD_R, CC_R, CL_R, Cl_R, Cm_R, Cn_R, ...
                           CD_delta_a, CC_delta_a, CL_delta_a, Cl_delta_a, Cm_delta_a, Cn_delta_a, ...
                           CD_delta_e, CC_delta_e, CL_delta_e, Cl_delta_e, Cm_delta_e, Cn_delta_e, ...
                           CD_delta_r, CC_delta_r, CL_delta_r, Cl_delta_r, Cm_delta_r, Cn_delta_r], ...
                          "VariableNames", {'CD', 'CC', 'CL', 'Cl', 'Cm', 'Cn', ...
                                            'CD_a', 'CC_a', 'CL_a', 'Cl_a', 'Cm_a', 'Cn_a', ...
                                            'CD_a_dot', 'CL_a_dot', 'Cm_a_dot', ...
                                            'CC_b_dot', 'Cn_b_dot', ...
                                            'CD_b', 'CC_b', 'CL_b', 'Cl_b', 'Cm_b', 'Cn_b', ...
                                            'CD_P', 'CC_P', 'CL_P', 'Cl_P', 'Cm_P', 'Cn_P', ...
                                            'CD_Q', 'CC_Q', 'CL_Q', 'Cl_Q', 'Cm_Q', 'Cn_Q', ...
                                            'CD_R', 'CC_R', 'CL_R', 'Cl_R', 'Cm_R', 'Cn_R', ...
                                            'CD_delta_a', 'CC_delta_a', 'CL_delta_a', 'Cl_delta_a', 'Cm_delta_a', 'Cn_delta_a', ...
                                            'CD_delta_e', 'CC_delta_e', 'CL_delta_e', 'Cl_delta_e', 'Cm_delta_e', 'Cn_delta_e', ...
                                            'CD_delta_r', 'CC_delta_r', 'CL_delta_r', 'Cl_delta_r', 'Cm_delta_r', 'Cn_delta_r'});

end
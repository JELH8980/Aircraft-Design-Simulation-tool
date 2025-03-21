function [Alon, Blon, Alat, Blat] = create_sf_lti(model, condition, derivatives)
% CREATE_SF_LTI - Function to create state-space linear time-invariant (LTI) models for longitudinal and lateral dynamics.
%
% This function constructs state-space matrices (A and B) for longitudinal (Alon, Blon) 
% and lateral (Alat, Blat) dynamics of an aircraft, based on a linearized model at a 
% given flight condition. It uses model parameters, atmospheric conditions, and 
% aerodynamic derivatives to compute stability and control derivatives, forming 
% the LTI system for small perturbation analysis.
%
% INPUTS:
%   model          - A struct containing aircraft parameters: mass (m), moments of 
%                    inertia (Ixx, Iyy, Izz, Ixz), reference areas (S, b, c), thrust 
%                    angle (aT), and geometric distances (hB, hR, lB, lR, hT, lT).
%   condition      - A struct with flight condition data: altitude (h [m]), Mach 
%                    number (M), initial velocity (u [m/s]), and pitch angle (theta [rad]).
%   derivatives    - A struct containing aerodynamic derivatives (e.g., CX, CZ, Cm, 
%                    CY, Cl, Cn and their partial derivatives with respect to state 
%                    and control variables).
%
% OUTPUTS:
%   Alon           - Longitudinal state matrix (4x4) [u, w, q, theta].
%   Blon           - Longitudinal control matrix (4x2) [delta_e, delta_t].
%   Alat           - Lateral state matrix (4x4) [v, p, r, phi].
%   Blat           - Lateral control matrix (4x2) [delta_r, delta_a].
%
% FUNCTIONALITY:
% - Extracts physical parameters (mass, inertia, geometry) and computes atmospheric 
%   properties (rho, a) at altitude h using `atmosphere_state`.
% - Calculates dynamic pressure (Q) and airspeed (V) from Mach number and speed of sound.
% - Defines geometric terms (a1, a2, a3, a4) for moment calculations.
% - Computes longitudinal stability derivatives (Xu, Xw, Zu, Zw, etc.) and control 
%   derivatives (Xdelta_e, Zdelta_t, etc.) for forces (X, Z) and moment (M).
% - Computes lateral stability derivatives (Yv, Yp, Yr, Lv, Lp, etc.) and control 
%   derivatives (Ldelta_a, Ndelta_r, etc.) for forces (Y) and moments (L, N).
% - Constructs state-space matrices:
%   - Longitudinal: Accounts for mass and inertia effects via Inlon matrix.
%   - Lateral: Includes cross-coupling via Ixz in Inlat matrix.
% - Applies inverse inertia matrices to finalize Alon, Blon, Alat, and Blat.
%
% NOTES:
% - Assumes a standard gravitational acceleration (g0 = 9.80665 m/s^2).
% - Longitudinal states: u (axial velocity), w (vertical velocity), q (pitch rate), 
%   theta (pitch angle). Controls: delta_e (elevator), delta_t (thrust).
% - Lateral states: v (lateral velocity), p (roll rate), r (yaw rate), phi (roll angle). 
%   Controls: delta_r (rudder), delta_a (aileron).
% - Units: SI (meters, seconds, kilograms, radians), except where derivatives use degrees 
%   (e.g., cosd, sind).
% - Depends on `atmosphere_state` for atmospheric properties.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

g0  = 9.80665;

m   = model.parameters.m;
Ix  = model.parameters.Ixx;
Iy  = model.parameters.Iyy;
Iz  = model.parameters.Izz;
Ixz = model.parameters.Ixz;

S   = model.parameters.S;
b   = model.parameters.b;
c   = model.parameters.c;

aT  = model.parameters.aT;

h = condition.("h [m]");
   
[~, ~, rho, ~, ~, a] = atmosphere_state(h);

M = condition.M;

V = M*a;

Q   = 1/2*rho*V^2;

a1 = model.parameters.hB-model.parameters.hR;

a2 = model.parameters.lR-model.parameters.lB;

a3 = cosd(aT)*(model.parameters.hB-model.parameters.hT) + sind(aT)*(model.parameters.lB-model.parameters.lT);

a4 = model.parameters.lB-model.parameters.lR;


u0           = condition.("u [m/s]");
theta0       = condition.("theta [rad]");

CT_M         = derivatives.CT_M;
CT_delta_t   = derivatives.CT_delta_t;

CX           = derivatives.CX;
CX_u         = derivatives.CX_u;
CX_a         = derivatives.CX_a;
CX_delta_e   = derivatives.CX_delta_e;


Xu           = 1/m*(Q*S/u0)*(2*CX + u0*CX_u + M/(Q*S)*cosd(aT)*CT_M);

XuA          = 1/m*(Q*S/u0)*(2*CX + u0*CX_u);

Xw           = 1/m*(Q*S/u0)*CX_a;

Xdelta_e     = 1/m*Q*S*CX_delta_e;

Xdelta_t     = 1/m*cosd(aT)*CT_delta_t;



CZ           = derivatives.CZ;
CZ_u         = derivatives.CZ_u;
CZ_a         = derivatives.CZ_a;
CZ_a_dot     = derivatives.CZ_a_dot;
CZ_Q         = derivatives.CZ_Q;
CZ_delta_e   = derivatives.CX_delta_e;


Zu           = 1/m*(Q*S/u0)*(2*CZ + u0*CZ_u - M/(Q*S)*sind(aT)*CT_M);

ZuA          = 1/m*(Q*S/u0)*(2*CZ + u0*CZ_u);

Zw           = 1/m*(Q*S/u0)*CZ_a;

Zwd          = 1/m*(Q*S/u0)*(c/(2*u0))*CZ_a_dot;

Zq           = 1/m*Q*S*(c/(2*u0))*CZ_Q;

Zdelta_e     = 1/m*Q*S*CZ_delta_e;

Zdelta_t     =-1/m*sind(aT)*CT_delta_t;



Cm           = derivatives.Cm;

Cm_u         = derivatives.Cm_u;

Cm_a         = derivatives.Cm_a;

Cm_Q         = derivatives.Cm_Q;

Cm_a_dot     = derivatives.Cm_a_dot;

Cm_delta_e   = derivatives.Cm_delta_e;


Mu           = 1/Iy*(a1*m*XuA + a2*m*ZuA + a3*M/u0*CT_M + (Q*S*c/u0)*(Cm + u0*Cm_u));

Mw           = 1/Iy*(a1*m*Xw + a2*m*Zw + (Q*S*c)/u0*Cm_a);

Mwd          = 1/Iy*(a2*m*Zwd + (Q*S*c^2)/(2*u0^2)*Cm_a_dot);

Mq           = 1/Iy*(a2*m*Zq + (Q*S*c^2)/(2*u0)*Cm_Q);

Mdelta_t     = 1/Iy*a3*CT_delta_t;

Mdelta_e     = 1/Iy*(a1*m*Xdelta_e + a2*m*Zdelta_e + Q*S*c*Cm_delta_e);
    


CY_b         = derivatives.CY_b;

CY_P         = derivatives.CY_P;

CY_R         = derivatives.CY_R;

CY_delta_r   = derivatives.CY_delta_r;


Yv           = 1/m*Q*S/u0*CY_b;

Yp           = 1/m*Q*S*b/(2*u0)*CY_P;

Yr           = 1/m*Q*S*b/(2*u0)*CY_R;

Ydelta_r     = 1/m*Q*S*CY_delta_r;



Cl_b         = derivatives.Cl_b;

Cl_P         = derivatives.Cl_P;

Cl_R         = derivatives.Cl_R;

Cl_delta_a   = derivatives.Cl_delta_a;

Cl_delta_r   = derivatives.Cl_delta_r;


Lv          = 1/Ix*Q*S*b/u0*Cl_b;

Lp          = 1/Ix*Q*S*b^2/(2*u0)*Cl_P;

Lr          = 1/Ix*Q*S*b^2/(2*u0)*Cl_R;

Ldelta_a    = 1/Ix*Q*S*b*Cl_delta_a;

Ldelta_r    = 1/Ix*Q*S*b*Cl_delta_r;



Cn_b         = derivatives.Cn_b;

Cn_P         = derivatives.Cn_P;

Cn_R         = derivatives.Cn_R;

Cn_delta_a   = derivatives.Cn_delta_a;

Cn_delta_r   = derivatives.Cn_delta_r;

Nv           = 1/Iz*(Q*S*b/u0*Cn_b + a4*m*Yv);

Np           = 1/Iz*Q*S*b^2/(2*u0)*Cn_P;

Nr           = 1/Iz*Q*S*b^2/(2*u0)*Cn_R;

Ndelta_a     = 1/Iz*Q*S*b*Cn_delta_a;

Ndelta_r     = 1/Iz*Q*S*b*Cn_delta_r;   


% Longitudinal 


Inlon        = [1,     0, 0, 0;
                0, 1-Zwd, 0, 0;
                0,  -Mwd, 1, 0;
                0,     0, 0, 1];

Alon_        = [Xu, Xw,     0, -g0*cos(theta0);
                Zu, Zw, u0+Zq, -g0*sin(theta0);
                Mu, Mw,    Mq,               0;
                 0,  0,     1,               0];

Blon_        = [Xdelta_e, Xdelta_t;
                Zdelta_e, Zdelta_t;
                Mdelta_e, Mdelta_t;
                       0,        0];

Inlon_inv = inv(Inlon);

Alon = Inlon_inv*Alon_;

Blon = Inlon_inv*Blon_;



% Lateral 
Inlat        = [  1,      0,      0,       0;
                  0,      1,-Ixz/Ix,       0
                  0,-Ixz/Iz,      1,       0;
                  0,      0,      0,       1];

Alat_        = [Yv, Yp, Yr-u0,  g0*cos(theta0);
                Lv, Lp,    Lr,               0;
                Nv, Np,    Nr,               0;
                 0,  1,     0,               0];

Blat_        = [Ydelta_r,        0;
                Ldelta_r, Ldelta_a;
                Ndelta_r, Ndelta_a;
                       0,        0];

Inlat_inv = inv(Inlat);

Alat = Inlat_inv*Alat_;

Blat = Inlat_inv*Blat_;


end
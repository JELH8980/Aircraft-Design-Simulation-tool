function FMAero = AeroFM_NL(qbar, V, alpha, alpha_dot, beta, beta_dot, P, Q, R, delta_e, delta_a, delta_r)
% AEROFM_NL - Computes nonlinear aerodynamic forces and moments in the body frame.
%
% This function calculates the aerodynamic forces (X, Y, Z) and moments (L, M, N) in the body frame for a given flight condition 
% using a nonlinear aerodynamic model. It incorporates stability and control derivatives, dynamic pressure, and dimensional 
% scaling based on reference areas and lengths. The forces and moments account for contributions from angle of attack, sideslip, 
% angular rates, control surface deflections, and their rates of change.
%
% INPUTS:
%   qbar       - Dynamic pressure (Pa), computed as 1/2 * rho * V^2.
%   V          - Total velocity (m/s).
%   alpha      - Angle of attack (rad).
%   alpha_dot  - Rate of change of angle of attack (rad/s).
%   beta       - Sideslip angle (rad).
%   beta_dot   - Rate of change of sideslip angle (rad/s).
%   P          - Roll rate (rad/s).
%   Q          - Pitch rate (rad/s).
%   R          - Yaw rate (rad/s).
%   delta_e    - Elevator deflection (rad).
%   delta_a    - Aileron deflection (rad).
%   delta_r    - Rudder deflection (rad).
%
% OUTPUTS:
%   FMAero     - A 6-element vector of aerodynamic forces and moments in the body frame: [X, Y, Z, L, M, N]', where:
%                - X, Y, Z: forces in body frame (N)
%                - L, M, N: moments in body frame (N·m)
%
% FUNCTIONALITY:
% - Initializes persistent variables (geometry and aerodynamic coefficients) on the first call from 'model' and 'interpolant' 
%   structs in the base workspace.
% - Computes wind-frame forces (D: drag, C: side force, L: lift) using interpolants for coefficients (e.g., CD, CL) and their 
%   derivatives with respect to state (alpha, beta, etc.) and control inputs (delta_e, delta_a, delta_r).
% - Applies dimensional scaling using dynamic pressure (qbar), reference area (S), wing span (b), and mean chord (c).
% - Transforms wind-frame forces to body-frame forces (X, Y, Z) using alpha and beta.
% - Computes moments (L, M, N) with additional contributions from force offsets based on reference point distances (lR, hR, lB, hB).
% - Returns the combined force and moment vector for use in 6-DOF dynamics.
%
% NOTES:
% - Persistent variables include geometric parameters (S, b, c, lR, hR, lB, hB) and aerodynamic interpolants (e.g., CD, CL_a).
% - Interpolants are assumed to be function handles or lookup tables from 'interpolant' that accept alpha and delta_e as inputs.
% - Nondimensional rates (e.g., c*Q/(2*V)) are used to scale angular rate and rate-of-change effects.
% - The function assumes a standard aerospace coordinate system with forces transformed from wind to body frame.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

    persistent initialized_flag ...
               S b c lR hR lB hB ...
               delta_a0 delta_r0 ...
               CD CC CL Cl Cm Cn ...
               CD_a CD_a_dot CD_b CD_P CD_Q CD_R CD_delta_a CD_delta_e CD_delta_r ...
               CC_a CC_b CC_b_dot CC_P CC_Q CC_R CC_delta_a CC_delta_e CC_delta_r ...
               CL_a CL_a_dot CL_b CL_P CL_Q CL_R CL_delta_a CL_delta_e CL_delta_r ...
               Cl_a Cl_b Cl_P Cl_Q Cl_R Cl_delta_a Cl_delta_e Cl_delta_r ...
               Cm_a Cm_a_dot Cm_b Cm_P Cm_Q Cm_R Cm_delta_a Cm_delta_e Cm_delta_r ...
               Cn_a Cn_b Cn_b_dot Cn_P Cn_Q Cn_R Cn_delta_a Cn_delta_e Cn_delta_r ...


               
    if isempty(initialized_flag)
        model               = evalin("base", 'model');
        S                   = model.parameters.S; 
        b                   = model.parameters.b;
        c                   = model.parameters.c;
        lR                  = model.parameters.lR;
        hR                  = model.parameters.hR;
        lB                  = model.parameters.lB;
        hB                  = model.parameters.hB;

        interpolant           = evalin("base", 'interpolant');
        % Stability and control derivatives
        CD                    = interpolant.CD;
        CD_a_dot              = interpolant.CD_a_dot;
        CD_b                  = interpolant.CD_b;
        CD_P                  = interpolant.CD_P;
        CD_Q                  = interpolant.CD_Q;
        CD_R                  = interpolant.CD_R;
        CD_delta_a            = interpolant.CD_delta_a;
        CD_delta_r            = interpolant.CD_delta_r;
        CC                    = interpolant.CC;
        CC_b                  = interpolant.CC_b;
        CC_b_dot              = interpolant.CC_b_dot;
        CC_P                  = interpolant.CC_P;
        CC_Q                  = interpolant.CC_Q;
        CC_R                  = interpolant.CC_R;
        CC_delta_a            = interpolant.CC_delta_a;
        CC_delta_r            = interpolant.CC_delta_r;
        CL                    = interpolant.CL;
        CL_a_dot              = interpolant.CL_a_dot;
        CL_b                  = interpolant.CL_b;
        CL_P                  = interpolant.CL_P;
        CL_Q                  = interpolant.CL_Q;
        CL_R                  = interpolant.CL_R;
        CL_delta_a            = interpolant.CL_delta_a;
        CL_delta_r            = interpolant.CL_delta_r;
        Cl                    = interpolant.Cl;
        Cl_b                  = interpolant.Cl_b;
        Cl_P                  = interpolant.Cl_P;
        Cl_Q                  = interpolant.Cl_Q;
        Cl_R                  = interpolant.Cl_R;
        Cl_delta_a            = interpolant.Cl_delta_a;
        Cl_delta_r            = interpolant.Cl_delta_r;
        Cm                    = interpolant.Cm;
        Cm_a_dot              = interpolant.Cm_a_dot;
        Cm_b                  = interpolant.Cm_b;
        Cm_P                  = interpolant.Cm_P;
        Cm_Q                  = interpolant.Cm_Q;
        Cm_R                  = interpolant.Cm_R;
        Cm_delta_a            = interpolant.Cm_delta_a;
        Cm_delta_r            = interpolant.Cm_delta_r;
        Cn                    = interpolant.Cn;
        Cn_a                  = interpolant.Cn_a;
        Cn_b                  = interpolant.Cn_b;
        Cn_b_dot              = interpolant.Cn_b_dot;
        Cn_P                  = interpolant.Cn_P;
        Cn_Q                  = interpolant.Cn_Q;
        Cn_R                  = interpolant.Cn_R;
        Cn_delta_a            = interpolant.Cn_delta_a;
        Cn_delta_r            = interpolant.Cn_delta_r;
        initialized_flag    = 1;
        

    end

    Dadot    = alpha_dot;

    Db       = beta;

    Dbdot    = beta_dot;

    DP       = P;

    DQ       = Q;

    DR       = R;

    Ddelta_a = delta_a;

    Ddelta_r = delta_r;

    D   = qbar * S * (CD(alpha, delta_e) + CD_Q(alpha, delta_e) * (c*DQ/(2*V)));
    C   = qbar * S * (CC_b(alpha, delta_e) * Db + CC_P(alpha, delta_e) * (b*DP/(2*V)) + CC_R(alpha, delta_e) * (b*DR/(2*V)) + CC_delta_a(alpha, delta_e) * Ddelta_a + CC_delta_r(alpha, delta_e) * Ddelta_r);
    L   = qbar * S * (CL(alpha, delta_e) + CL_b(alpha, delta_e) * Db + CL_Q(alpha, delta_e) * (c*DQ/(2*V)));
    
    
    X = -cos(alpha)*cos(beta)*D - cos(alpha)*sin(beta)*C + sin(alpha)*L;
    Y = -sin(beta)*D + cos(beta)*C;
    Z = -sin(alpha)*cos(beta)*D - sin(alpha)*sin(beta)*C - cos(alpha)*L;

    L   = qbar * S * b * (Cl_b(alpha, delta_e) * Db + Cl_P(alpha, delta_e) * (b*DP/(2*V)) + Cl_R(alpha, delta_e) * (b*DR/(2*V)) + Cl_delta_a(alpha, delta_e) * Ddelta_a + Cl_delta_r(alpha, delta_e) * Ddelta_r);
    M   = qbar * S * c * (Cm(alpha, delta_e) + Cm_a_dot(alpha, delta_e) * (c*Dadot/(2*V)) + Cm_Q(alpha, delta_e) * (c*DQ/(2*V)));
    N   = qbar * S * b * (Cn_b(alpha, delta_e) * Db + Cn_b_dot(alpha, delta_e) * (b*Dbdot/(2*V)) + Cn_P(alpha, delta_e) * (b*DP/(2*V)) + Cn_R(alpha, delta_e) * (b*DR/(2*V)) + Cn_delta_a(alpha, delta_e) * Ddelta_a + Cn_delta_r(alpha, delta_e) * Ddelta_r);

    L   = L + Y*(hR-hB);
    M   = M - X*(hR-hB) - Z*(lB-lR);
    N   = N + Y*(lB-lR);

    FMAero = [X, Y, Z, L, M, N]';
 
end
function x_dot = linear6dof_simulink(in)
    % LINEAR6DOF_SIMULINK - Simulates a linearized 6-DOF (Degrees of Freedom) dynamic model in Simulink.
%
% This function implements a linearized 6-DOF model for simulating the dynamics of a system, such as an aircraft or vehicle, 
% using longitudinal and lateral state-space models. It computes the time derivatives of the state variables based on 
% initial conditions, control inputs, and pre-defined system matrices (Alon, Blon, Alat, Blat). The model includes 
% translational and rotational kinematics to update position and orientation in a 3D space.
%
% INPUTS:
%   in         - A 16-element vector containing the current state (12 elements) and control inputs (4 elements).
%                - in(1:12) : State vector [u, v, w, p, q, r, e1, e2, e3, xE, yE, zE]
%                - in(13:16): Control inputs [delta_e, delta_a, delta_r, Tcom]
%
% OUTPUTS:
%   x_dot      - A 12-element vector of state derivatives [u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, 
%                e1_dot, e2_dot, e3_dot, xE_dot, yE_dot, zE_dot].
%
% FUNCTIONALITY:
% - Initializes persistent variables (system matrices, initial conditions) on the first call using data from the base workspace.
% - Extracts the current state and control inputs from the input vector.
% - Computes deviations from initial conditions for longitudinal and lateral dynamics.
% - Applies linearized state-space equations to calculate state derivatives for longitudinal (u, w, q, e2) 
%   and lateral (v, p, r, e1) dynamics.
% - Computes rotational kinematics for e3_dot (yaw rate) and translational kinematics for position derivatives (xE_dot, yE_dot, zE_dot).
% - Returns the full state derivative vector for use in Simulink integration.
%
% NOTES:
% - Persistent variables are used to store system matrices (Alon, Blon, Alat, Blat) and initial conditions 
%   (u0, v0, w0, etc.) between calls, initialized only once.
% - The function assumes the existence of 'Alon', 'Blon', 'Alat', 'Blat', and 'reference' structs in the base workspace.
% - Euler angles (e1, e2, e3) represent roll (phi), pitch (theta), and yaw (psi), respectively.
% - Control inputs are elevator (delta_e), aileron (delta_a), rudder (delta_r), and thrust command (Tcom).
%
% Author: Ludwig Horvath
% Date: 2/11/2025

    
    persistent initialized_flag ...
               Alon Blon Alat Blat ...
               u0 v0 w0 p0 q0 r0 e10 e20 ...
               delta_e0 delta_a0 delta_r0 Tcom0

    if isempty(initialized_flag)
        % Expecting the second argument to be a struct called 'setup'
        Alon                = evalin("base", 'Alon');
        Blon                = evalin("base", 'Blon');
        Alat                = evalin("base", 'Alat');
        Blat                = evalin("base", 'Blat');

        initialized_flag    = 1;
        
        reference           = evalin("base", 'reference');
        u0                  = reference.u;
        v0                  = reference.v;
        w0                  = reference.w;
        p0                  = reference.P;
        q0                  = reference.Q;
        r0                  = reference.R;
        e10                 = reference.phi;
        e20                 = reference.theta;
        delta_e0            = reference.delta_e;
        delta_a0            = reference.delta_a;
        delta_r0            = reference.delta_r;
        Tcom0               = reference.Tcom;

    end

    x = in(1:12);

    ctrl = in(13:16);


    delta_e = ctrl(1);
    delta_a = ctrl(2);
    delta_r = ctrl(3);
    Tcom    = ctrl(4);
    
    u = x(1);
    v = x(2);
    w = x(3);

    V = sqrt(u^2 + v^2 + w^2);
    
    p = x(4);
    q = x(5);
    r = x(6);
    
    e1 = x(7);
    e2 = x(8);
    e3 = x(9);
    
    xE = x(10);
    yE = x(11);
    zE = x(12);

    delta_xlon = [u-u0;
                  w-w0;
                  q-q0;
                  e2-e20];

    delta_ulon = [delta_e - delta_e0;
                  Tcom - Tcom0];

    xdotlon = Alon*delta_xlon + Blon*delta_ulon;

    delta_xlat = [v-v0;
                  p-p0;
                  r-r0;
                  e1-e10];

    delta_ulat = [delta_r - delta_r0;
                  delta_a - delta_a0];

    xdotlat = Alat*delta_xlat + Blat*delta_ulat;

    u_dot     = xdotlon(1);
    w_dot     = xdotlon(2);
    q_dot     = xdotlon(3);
    e2_dot    = xdotlon(4);

    v_dot     = xdotlat(1);
    p_dot     = xdotlat(2);
    r_dot     = xdotlat(3);
    e1_dot    = xdotlat(4);

    % Rotational kinematics
    e3_dot = (q*sin(e1) + r*cos(e1))/cos(e2);
    
    % Translational kinematics
    xE_dot = cos(e2)*cos(e3)*u + (sin(e1)*sin(e2)*cos(e3) - cos(e1)*sin(e3))*v + (cos(e1)*sin(e2)*cos(e3) + sin(e1)*sin(e3))*w;
    yE_dot = cos(e2)*sin(e3)*u + (sin(e1)*sin(e2)*sin(e3) + cos(e1)*cos(e3))*v + (cos(e1)*sin(e2)*sin(e3) - sin(e1)*cos(e3))*w;
    zE_dot =      (-sin(e2))*u +                             sin(e1)*cos(e2)*v +                             cos(e1)*cos(e2)*w;
    
    x_dot = [ u_dot; 
              v_dot;
              w_dot;
              p_dot;
              q_dot;
              r_dot;
             e1_dot;
             e2_dot;
             e3_dot;
             xE_dot;
             yE_dot;
             zE_dot];


end
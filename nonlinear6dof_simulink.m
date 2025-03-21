function x_dot = nonlinear6dof_simulink(in)
    % NONLINEAR6DOF_SIMULINK - Simulates a nonlinear 6-DOF (Degrees of Freedom) dynamic model in Simulink with control limits.
    %
    % This function implements a nonlinear 6-DOF model for simulating the dynamics of a system, such as an aircraft or vehicle, 
    % in Simulink. It computes the time derivatives of the state variables based on control inputs, aerodynamic forces, thrust, 
    % gravity, and inertial properties, while enforcing limits on control surface deflections and thrust. The model includes 
    % translational and rotational dynamics, as well as kinematics, to update velocity, angular rates, orientation, and position 
    % in a 3D space.
    %
    % INPUTS:
    %   in         - A 16-element vector containing the current state (12 elements) and control inputs (4 elements):
    %                - in(1:12) : State vector [u, v, w, p, q, r, e1, e2, e3, xE, yE, zE], where:
    %                  - u, v, w: velocities in body frame (m/s)
    %                  - p, q, r: angular rates (roll, pitch, yaw) in body frame (rad/s)
    %                  - e1, e2, e3: Euler angles (roll, pitch, yaw) in radians
    %                  - xE, yE, zE: positions in Earth frame (m)
    %                - in(13:16): Control inputs [delta_e, delta_a, delta_r, Tcom], where:
    %                  - delta_e: elevator deflection (rad)
    %                  - delta_a: aileron deflection (rad)
    %                  - delta_r: rudder deflection (rad)
    %                  - Tcom: thrust command (%)
    %
    % OUTPUTS:
    %   x_dot      - A 12-element vector of state derivatives: [u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, 
    %                e1_dot, e2_dot, e3_dot, xE_dot, yE_dot, zE_dot].
    %
    % FUNCTIONALITY:
    % - Initializes persistent variables (mass, inertia properties, and control limits) on the first call from a 'model' struct 
    %   in the base workspace.
    % - Applies saturation limits to control inputs (delta_e, delta_a, delta_r, Tcom) based on predefined max/min values.
    % - Computes atmospheric conditions (density, speed of sound) based on altitude (h = -zE).
    % - Calculates total velocity (V), angle of attack (alpha), sideslip angle (beta), Mach number (M), and dynamic pressure (qbar).
    % - Evaluates forces and moments from thrust, gravity, and aerodynamics using external functions (ThrustFM, GravityFM, AeroFM_NL).
    % - Computes translational dynamics (u_dot, v_dot, w_dot) using Newton’s equations with Coriolis terms.
    % - Computes rotational dynamics (p_dot, q_dot, r_dot) using Euler’s equations with inertia coupling.
    % - Computes rotational kinematics (e1_dot, e2_dot, e3_dot) for Euler angle rates.
    % - Computes translational kinematics (xE_dot, yE_dot, zE_dot) for position rates in the Earth frame.
    % - Returns the full state derivative vector for Simulink integration.
    %
    % NOTES:
    % - Persistent variables (m, Ixx, Iyy, Izz, Ixz, and control limits) are initialized once from the 'model.parameters' struct.
    % - Control limits for delta_e, delta_a, and delta_r are converted from degrees to radians during initialization.
    % - Thrust command (Tcom) is constrained between 0% and 100%.
    % - External functions (atmosphere_state, ThrustFM, GravityFM, AeroFM_NL) are assumed to be defined elsewhere.
    % - Euler angles (e1, e2, e3) represent roll (phi), pitch (theta), and yaw (psi), respectively.
    %
    % Author: Ludwig Horvath
    % Date: 2/11/2025

    persistent initialized_flag m Ixx Iyy Izz Ixz ...
                delta_e_max delta_e_min ...
                delta_a_max delta_a_min ...
                delta_r_max delta_r_min

    if isempty(initialized_flag)
        % Expecting the second argument to be a struct called 'setup'
        model               = evalin("base", 'model');
        Ixx                 = model.parameters.Ixx;
        Iyy                 = model.parameters.Iyy;
        Izz                 = model.parameters.Izz;
        Ixz                 = model.parameters.Ixz;
        m                   = model.parameters.m;
        delta_e_max         = deg2rad(model.parameters.delta_e.max);
        delta_e_min         = deg2rad(model.parameters.delta_e.min);
        delta_a_max         = deg2rad(model.parameters.delta_a.max);
        delta_a_min         = deg2rad(model.parameters.delta_a.min);
        delta_r_max         = deg2rad(model.parameters.delta_r.max);
        delta_r_min         = deg2rad(model.parameters.delta_r.min);

        initialized_flag    = 1;
        
    end

    x    = in(1:12);

    ctrl = in(13:16); 

    delta_e = ctrl(1);

    if delta_e > delta_e_max
        delta_e = delta_e_max;
    elseif delta_e < delta_e_min
        delta_e = delta_e_min;
    end

    delta_a = ctrl(2);
    
    if delta_a > delta_a_max
        delta_a = delta_a_max;
    elseif delta_a < delta_a_min
        delta_a = delta_a_min;
    end

    delta_r = ctrl(3);

    if delta_r > delta_r_max
        delta_r = delta_r_max;
    elseif delta_r < delta_r_min
        delta_r = delta_r_min;
    end

    Tcom    = ctrl(4);
    
    if Tcom > 100
        Tcom = 100;
    elseif Tcom < 0
        Tcom = 0;
    end

    u = x(1);
    v = x(2);
    w = x(3);
    
    p = x(4);
    q = x(5);
    r = x(6);
    
    e1 = x(7);
    e2 = x(8);
    e3 = x(9);
    
    xE = x(10);
    yE = x(11);
    zE = x(12);

    V = sqrt(u^2 + v^2 + w^2);
    
    alpha = atan2(w, u);
    
    beta = asin(v/V);
   
    h = -zE;
    
    [~, ~, rho, ~, ~, a] = atmosphere_state(h);
    
    M = V/a;
    
    qbar = 1/2*rho*V^2;

    FMthrust  = ThrustFM(Tcom, M, h);
    
    FMgravity = GravityFM(e1, e2);

    FMaero    = AeroFM_NL(qbar, V, alpha, 0, beta, 0, p, q, r, delta_e, delta_a, delta_r);

    Force     = FMthrust(1:3) + FMgravity(1:3) + FMaero(1:3);
    
    Moment    = FMthrust(4:6) + FMgravity(4:6) + FMaero(4:6);
    

    X = Force(1);
    Y = Force(2); 
    Z = Force(3);
    
    L = Moment(1);
    M = Moment(2); 
    N = Moment(3);

 
    
    % Translational dynamics
    u_dot = 1/m*X - q*w + r*v;
    v_dot = 1/m*Y - r*u + p*w;
    w_dot = 1/m*Z - p*v + q*u;
    
    % Rotational dynamics
    p_dot = (Izz*L + Ixz*N + (Izz - Ixz)*q*r + Iyy*q*r)/(Ixx*Izz - Ixz^2);
    q_dot = (M + Ixz*(r^2 - p^2) + (Izz - Ixx)*p*r)/Iyy;
    r_dot = (Ixz*L + Ixx*N + (Ixx - Ixz)*p*q - Iyy*p*q)/(Ixx*Izz - Ixz^2);
    
    % Rotational kinematics
    e1_dot = p + (q*sin(e1) + r*cos(e1))*tan(e2);
    e2_dot = q*cos(e1) - r*sin(e1);
    e3_dot = (q*sin(e1) + r*cos(e1))/cos(e2);
    
    % Translational kinematics
    xE_dot = cos(e2)*cos(e3)*u + (sin(e1)*sin(e2)*cos(e3) - cos(e1)*sin(e3))*v + (cos(e1)*sin(e2)*cos(e3) + sin(e1)*sin(e3))*w;
    yE_dot = cos(e2)*sin(e3)*u + (sin(e1)*sin(e2)*sin(e3) + cos(e1)*cos(e3))*v + (cos(e1)*sin(e2)*sin(e3) - sin(e1)*cos(e3))*w;
    zE_dot = (-sin(e2))*u + sin(e1)*cos(e2)*v + cos(e1)*cos(e2)*w;
        

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
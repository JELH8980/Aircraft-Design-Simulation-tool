function ctrl = controller_1(x)
% CONTROLLER_1 - This function implements a simple state feedback controller for an aircraft
% using control inputs (elevator deflection, aileron deflection, and rudder deflection) based on 
% the aircraft's state vector. The controller computes the required control surface deflections 
% to achieve desired longitudinal and lateral dynamics, with constraints applied to the deflections.
%
% INPUT:
%   x       - Aircraft state vector [u, v, w, p, q, r, phi, theta], where:
%             u: forward velocity (m/s)
%             v: side velocity (m/s)
%             w: vertical velocity (m/s)
%             p: roll rate (rad/s)
%             q: pitch rate (rad/s)
%             r: yaw rate (rad/s)
%             phi: roll angle (rad)
%             theta: pitch angle (rad)
%
% OUTPUT:
%   ctrl    - Control surface deflections [delta_e, delta_a, delta_r, 0], where:
%             delta_e: elevator deflection (rad)
%             delta_a: aileron deflection (rad)
%             delta_r: rudder deflection (rad)
%             The last entry is placeholder (set to zero).
%
% PARAMETERS:
%   The function uses persistent variables that are initialized once during the first call:
%   - Maximum and minimum control surface deflections for elevator, aileron, and rudder.
%   - Controller gains (Kp, Kq, Kr, Kyar) for longitudinal and lateral dynamics.
%   - The model parameters are retrieved from the base workspace.
%
% FUNCTIONALITY:
% - The aircraft state vector is split into longitudinal and lateral components.
% - Longitudinal control is computed using a gain matrix `Klon`.
% - Lateral control is computed using a gain matrix `Klat`.
% - The control surface deflections (elevator, aileron, rudder) are calculated from the state.
% - Constraints are applied to each control surface deflection to ensure they stay within their limits.
% - The function outputs the control surface deflections as a control command vector.
%
% NOTES:
% - The controller computes control inputs for the aircraft's stability and maneuverability, 
%   taking into account the aircraft's current states like velocity and rotational rates.
% - The control deflections are bounded by maximum and minimum limits to prevent excessive control surface deflections.


persistent delta_e_max delta_e_min ...
           delta_a_max delta_a_min ... 
           delta_r_max delta_r_min ...
           Kq Kp Kr Kyar Klat Klon ...
           persistent_flag

if isempty(persistent_flag)

    model = evalin("base", "model");
    delta_e_max = model.parameters.delta_e.max;
    delta_e_min = model.parameters.delta_e.min;
    delta_a_max = model.parameters.delta_a.max;
    delta_a_min = model.parameters.delta_a.min;
    delta_r_max = model.parameters.delta_r.max;
    delta_r_min = model.parameters.delta_r.min;

    Kp          = evalin("base", "Kp");
    Kq          = evalin("base", "Kq");
    Kr          = evalin("base", "Kr");
    Kyar        = evalin("base", "Kyar");


    Klat = [0,  0,   Kr, 0; ...
            0, Kp, Kyar, 0];
    
    Klon = [0, 0, Kq, 0; ...
            0, 0,  0, 0];
    
    persistent_flag = 1;

end

d2r = pi/180;

u = x(1);
v = x(2);
w = x(3);

p = x(4);
q = x(5);
r = x(6);

phi = x(7);
theta = x(8);

xlon = [u; w; q; theta];
xlat = [v; p; r; phi];

ctrl_lon = Klon*xlon;

delta_e = ctrl_lon(1);

ctrl_lat = Klat*xlat;

delta_r = ctrl_lat(1);
delta_a = ctrl_lat(2);


if delta_e > d2r*delta_e_max
    delta_e = delta_e_max;
elseif delta_e < d2r*delta_e_min
    delta_e = delta_e_min;
end

if delta_a > d2r*delta_a_max
    delta_a = delta_a_max;
elseif delta_a < d2r*delta_a_min
    delta_a = delta_a_min;
end

if delta_r > d2r*delta_r_max
    delta_r = delta_r_max;
elseif delta_r < d2r*delta_r_min
    delta_r = delta_r_min;
end




ctrl = [delta_e; delta_a; delta_r; 0];

end
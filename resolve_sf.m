function result = resolve_sf(model, condition)

u       = condition.("u [m/s]");
w       = condition.("w [m/s]");
theta   = condition.("theta [rad]");
h       = condition.("h [m]");
delta_e = condition.("delta_e [rad]");
Tcom    = condition.("Tcom [%]");

delta_e_max         = deg2rad(model.parameters.delta_e.max);
delta_e_min         = deg2rad(model.parameters.delta_e.min);


spec.V     = sqrt(u^2+w^2);
spec.gamma = deg2rad(condition.("gamma [deg]"));
spec.h     = condition.("h [m]");

z0 = [u, 0, w, 0, 0, 0, 0, theta, 0, 0, 0, -h, delta_e, 0, 0, Tcom]';


initialization = 0;

f = inf;

tolerance = 10e-10;
maxiter = 10;

iter = 0;
while all([f > tolerance, iter < maxiter])

if (initialization == 0)
    z_guess = z0;
else
    z_guess = zstar;

end

[zstar, f] = fminsearch(@cost_straight_flight, z_guess, ...
                        optimset('TolX', 1e-10, 'MaxFunEvals', 10000, 'MaxIter', 10000));

initialization = 1;

xstar = zstar(1:12);
ctrlstar = zstar(13:16);

x_dotstar = nonlinear6dof(xstar, ctrlstar);

iter = iter + 1;

end

h = -xstar(12);

[~, ~, ~, ~, ~, a] = atmosphere_state(h);

u = xstar(1);
v = xstar(2);
w = xstar(3);

V       = sqrt(u^2 + v^2 + w^2);
M       = V/a;

theta   = xstar(8);

alpha   = atan2(w, u);
gamma   = rad2deg(theta-alpha);

delta_e = ctrlstar(1);

Tcom    = ctrlstar(4);

addpath("aircraft\propulsion\")
method    = evalin("base", 'model.parameters.propulsion.method');
engine_model        = str2func(strrep(method, '.m', ''));  % Convert function name to handle
Tmax    = evalin("base", 'model.parameters.propulsion.Tmax');
T       = engine_model(Tmax, Tcom, M, h);


condition.M = M;
condition.("h [m]") = h;
condition.("u [m/s]") = u;
condition.("theta [rad]") = theta;
condition.("delta_e [rad]") = delta_e;
condition.("T [N]") = T;
condition.("Tcom [%]") = Tcom;
condition.("gamma [deg]") = gamma;

result.condition = condition;
result.x = xstar;
result.ctrl = ctrlstar;
result.f = f;


function f = cost_straight_flight(z)

    x    = z(1:12);
    ctrl = z(13:16);

    u = x(1);
    v = x(2);
    w = x(3);
    
    e1 = x(7);

    h = -x(12);
    
    x_dot = nonlinear6dof(x, ctrl);
    
    u_dot = x_dot(1);
    v_dot = x_dot(2);
    w_dot = x_dot(3);
    p_dot = x_dot(4);
    q_dot = x_dot(5);
    r_dot = x_dot(6);
    e1_dot = x_dot(7);
    e2_dot = x_dot(8);
    e3_dot = x_dot(9);
    
    V = sqrt(u^2+v^2+w^2);
    
    theta = x(8);
    
    alpha = atan2(w, u);
    
    gamma = theta - alpha;
    
    Q = [u_dot;
         v_dot;
         w_dot;
         p_dot;
         q_dot;
         r_dot;
         e1_dot;
         e2_dot;
         e3_dot;
         V-spec.V;
         gamma-spec.gamma;
         h-spec.h;
         v;
         e1;
         ctrl(2);
         ctrl(3)];
   
    H = eye(16);
    
    f = Q'*H*Q;

    if ctrl(1) > delta_e_max
        f = f + 1;
    elseif ctrl(1) < delta_e_min
        f = f + 1;
    end

    if ctrl(4) > 100
        f = f + 2;
    elseif ctrl(4) < 0
        f = f + 2;
    end

end

end
clear; close all; clc

model = load_model();

condition = load_sf_condition(model);

%%
result = resolve_sf(model, condition);

%%

condition = result.condition;

%%
derivatives = calculate_sf_derivatives(model, condition);

%%
derivatives_sd = calculate_sdsf_derivatives(model, condition);
[Alon, Blon, Alat, Blat] = create_sf_lti(model, condition, derivatives_sd);
%%
sf_lti_analysis(Alon, Blon, 'longitudinal');
sf_lti_analysis(Alat, Blat, 'lateral');

%%
batch = load_batch(model);

interpolant = create_interpolant(batch);

%%
[Kq, Kr, Kp, Kyar] = sf_lti_rlocus(Alon, Blon, Alat, Blat, [], [], [], []);
%%
reference = create_sf_reference();
%%

tau_p    = 0.5;
tau_q    = 0.5; 
tau_r    = 0.5;

Ksp      = 0.1;
Ksq      = 0.1;
Ksr      = 0.1; 

Kp_roll  = -1;
Kd_roll  = -10;

Kp_pitch = -1;
Ki_pitch = -0.1;

Kp_psi = 10;
Kd_psi = 10;
Ki_psi = 0.1;


%%

reference = create_sf_reference();


%%

delta_e_max         = deg2rad(model.parameters.delta_e.max);
delta_e_min         = deg2rad(model.parameters.delta_e.min);
delta_a_max         = deg2rad(model.parameters.delta_a.max);
delta_a_min         = deg2rad(model.parameters.delta_a.min);
delta_r_max         = deg2rad(model.parameters.delta_r.max);
delta_r_min         = deg2rad(model.parameters.delta_r.min);


%%

tsim = 100;
dt = 10e-4;

d2r = 180/pi;

tmax = 40;

% Define the maximum index for plotting
maxidx = ceil(tmax/dt); % Adjust this value to change the viewed region (e.g., 100, 200, etc.)


%%

ctrl = reference.ctrl;
x    = reference.x;

t = 0;

log_nl = zeros(17, ceil(tsim/dt));




try
    while t <= tsim

    ctrl = controller(x) + ctrl_ref; % Make controller
    
    xdot = nonlinear6dof(x, ctrl);
    
    x = x + dt*xdot;

    t = t + dt;

    log_nl(:,round((t+dt)/dt)) = [t; x(1:12); ctrl];
    
    end

catch ME
   disp('simulation was abrubtly cancelled')
   disp('Matlab Error Message: ')
   disp(ME)

   log_nl = log_nl(:,1:round((t+dt)/dt)-1);
end

% 
%%

log_nl = array2table(log_nl', 'VariableNames', {'t', ...
                                                'u', 'v', 'w', ...
                                                'p', 'q', 'r', ...
                                                'phi', 'theta', 'psi', ...
                                                'xE', 'yE', 'zE', ...
                                                'delta_e', 'delta_a', 'delta_r', 'Tcom'});


%%

% Figure 1: Nonlinear State-Time trajectories
figure('Name', 'Nonlinear State-Time Trajectories')

% Subplot 1: Translational velocities (u, v, w) - Nonlinear
subplot(6,1,1);
hold on
plot(log_nl.t(1:maxidx), log_nl.u(1:maxidx), 'color', 'r', 'DisplayName', 'u');
plot(log_nl.t(1:maxidx), log_nl.v(1:maxidx), 'color', 'b', 'DisplayName', 'v');
plot(log_nl.t(1:maxidx), log_nl.w(1:maxidx), 'color', 'g', 'DisplayName', 'w');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 2: Angular velocities (p, q, r) - Nonlinear
subplot(6,1,2);
hold on
plot(log_nl.t(1:maxidx), d2r*log_nl.p(1:maxidx), 'color', 'r', 'DisplayName', 'p');
plot(log_nl.t(1:maxidx), d2r*log_nl.q(1:maxidx), 'color', 'b', 'DisplayName', 'q');
plot(log_nl.t(1:maxidx), d2r*log_nl.r(1:maxidx), 'color', 'g', 'DisplayName', 'r');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 3: Euler angles (phi, theta, psi) - Nonlinear
subplot(6,1,3);
hold on
plot(log_nl.t(1:maxidx), d2r*log_nl.phi(1:maxidx), 'color', 'r', 'DisplayName', 'phi');
plot(log_nl.t(1:maxidx), d2r*log_nl.theta(1:maxidx), 'color', 'b', 'DisplayName', 'theta');
plot(log_nl.t(1:maxidx), d2r*log_nl.psi(1:maxidx), 'color', 'g', 'DisplayName', 'psi');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 4: Positions (xE, yE, zE) - Nonlinear
subplot(6,1,4);
hold on
plot(log_nl.t(1:maxidx), log_nl.xE(1:maxidx), 'color', 'r', 'DisplayName', 'x_E');
plot(log_nl.t(1:maxidx), log_nl.yE(1:maxidx), 'color', 'b', 'DisplayName', 'y_E');
plot(log_nl.t(1:maxidx), log_nl.zE(1:maxidx), 'color', 'g', 'DisplayName', 'z_E');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 5: Control Inputs (delta_e, delta_a, delta_r, Tcom)
subplot(6,1,5);
hold on
plot(log_nl.t(1:maxidx), log_nl.Tcom(1:maxidx), 'color', 'k', 'DisplayName', 'Tcom');


subplot(6,1,6);
hold on
plot(log_nl.t(1:maxidx), log_nl.delta_a(1:maxidx), 'color', 'r', 'DisplayName', 'delta_a');
plot(log_nl.t(1:maxidx), log_nl.delta_r(1:maxidx), 'color', 'g', 'DisplayName', 'delta_r');
plot(log_nl.t(1:maxidx), log_nl.delta_e(1:maxidx), 'color', 'b', 'DisplayName', 'delta_e');
hold off
legend('location', 'northwest', 'Box', 'on')

%%

ctrl_ref = reference.ctrl;
x        = reference.x;

t = 0;

log_l = zeros(17, ceil(tsim/dt));

try
    while t <= tsim

    
    ctrl = controller(x) + ctrl_ref; % Make controller

    xdot = linear6dof(x, ctrl);
    
    x = x + dt*xdot;
    
    t = t + dt;
    
    
    log_l(:,round((t+dt)/dt)) = [t; x(1:12); ctrl];
    
    end

catch ME
   disp('simulation was abrubtly cancelled')
   disp('Matlab Error Message: ')
   disp(ME)
   
   log_l = log_l(:,1:round((t+dt)/dt)-1);

end



%% Converting log to table format with labels

log_l = array2table(log_l', 'VariableNames', {'t', ...
                                              'u', 'v', 'w', ...
                                              'p', 'q', 'r', ...
                                              'phi', 'theta', 'psi', ...
                                              'xE', 'yE', 'zE', ...
                                              'delta_e', 'delta_a', 'delta_r', 'Tcom'});

%% Plotting signals

% Figure 2: Linear State-Time trajectories
figure('Name', 'Linear State-Time Trajectories')

% Subplot 1: Translational velocities (u, v, w) - Linear
subplot(6,1,1);
hold on
plot(log_l.t(1:maxidx), log_l.u(1:maxidx), 'color', 'r', 'DisplayName', 'u');
plot(log_l.t(1:maxidx), log_l.v(1:maxidx), 'color', 'b', 'DisplayName', 'v');
plot(log_l.t(1:maxidx), log_l.w(1:maxidx), 'color', 'g', 'DisplayName', 'w');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 2: Angular velocities (p, q, r) - Linear
subplot(6,1,2);
hold on
plot(log_l.t(1:maxidx), d2r*log_l.p(1:maxidx), 'color', 'r', 'DisplayName', 'p');
plot(log_l.t(1:maxidx), d2r*log_l.q(1:maxidx), 'color', 'b', 'DisplayName', 'q');
plot(log_l.t(1:maxidx), d2r*log_l.r(1:maxidx), 'color', 'g', 'DisplayName', 'r');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 3: Euler angles (phi, theta, psi) - Linear
subplot(6,1,3);
hold on
plot(log_l.t(1:maxidx), d2r*log_l.phi(1:maxidx), 'color', 'r', 'DisplayName', 'phi');
plot(log_l.t(1:maxidx), d2r*log_l.theta(1:maxidx), 'color', 'b', 'DisplayName', 'theta');
plot(log_l.t(1:maxidx), d2r*log_l.psi(1:maxidx), 'color', 'g', 'DisplayName', 'psi');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 4: Positions (xE, yE, zE) - Linear
subplot(6,1,4);
hold on
plot(log_l.t(1:maxidx), log_l.xE(1:maxidx), 'color', 'r', 'DisplayName', 'x_E');
plot(log_l.t(1:maxidx), log_l.yE(1:maxidx), 'color', 'b', 'DisplayName', 'y_E');
plot(log_l.t(1:maxidx), log_l.zE(1:maxidx), 'color', 'g', 'DisplayName', 'z_E');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 5: Control Inputs (delta_e, delta_a, delta_r, Tcom)
subplot(6,1,5);
hold on
plot(log_l.t(1:maxidx), log_l.Tcom(1:maxidx), 'color', 'k', 'DisplayName', 'Tcom');


subplot(6,1,6);
hold on
plot(log_l.t(1:maxidx), log_l.delta_a(1:maxidx), 'color', 'r', 'DisplayName', 'delta_a');
plot(log_l.t(1:maxidx), log_l.delta_r(1:maxidx), 'color', 'g', 'DisplayName', 'delta_r');
plot(log_l.t(1:maxidx), log_l.delta_e(1:maxidx), 'color', 'b', 'DisplayName', 'delta_e');
hold off
legend('location', 'northwest', 'Box', 'on')
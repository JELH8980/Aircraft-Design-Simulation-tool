clear; close all; clc

%% Loading Model YF-16 + Batch with derivatives + Creating Interpolant

load("aircraft\YF-16.mat", "model");

data = readtable("database\batches\YF-168332513177.csv");

load("database\batches\YF-168332513177.mat");

batch.data = data;
batch.info = info;

interpolant = create_interpolant(batch);

%% Loading a Straight Flight Condition of choice
load('prepared\M075h10000.mat', 'condition');

%% Creating a reference to initalize simulation capabiltiies + solver

latitude0    = 59.3458;
longitude0   = 18.0731;

autopilot.rollpitch = 1;
autopilot.heading = 1;

dt = 10e-4;
tsim = 180;
N = size(0:dt:tsim-dt);

test = 4;


max_time = 180;

% Case 1
switch test
    case 1
        signal_in_delta_e = generate_step(100, dt, 10, 11, deg2rad(21.5), false);
        signal_in_delta_a = generate_step(100, dt, 20, 21, deg2rad(21.5), false);
        signal_in_delta_r = generate_step(tsim, dt, 30, 31, deg2rad(21.5), false);
        signal_in_psi     = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_theta   = timeseries(ones(N)*condition.("theta [rad]"), 0:dt:tsim-dt);
    case 2
        % Case 2
        signal_in_delta_a = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_delta_e = generate_step(100, dt, 10, 11, deg2rad(21.5), false);
        signal_in_delta_r = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_psi     = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_theta   = timeseries(ones(N)*condition.("theta [rad]"), 0:dt:tsim-dt);
    case 3
        % Case 3
        signal_in_delta_e = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_delta_a = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_delta_r = generate_step(100, dt, 10, 11, deg2rad(21.5), false);
        signal_in_psi     = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_theta   = timeseries(ones(N)*condition.("theta [rad]"), 0:dt:tsim-dt);
    case 4
        % Case 4
        signal_in_delta_e = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_delta_a = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_delta_r = timeseries(zeros(N), 0:dt:tsim-dt);
        signal_in_psi     = generate_doublet(tsim, dt, 10, 100, tsim, deg2rad(10), false);
        signal_in_theta   = timeseries(ones(N)*condition.("theta [rad]"), 0:dt:tsim-dt);
end

reference = create_sf_reference(); 

%% Resolving the solution using interpolant + nonlinear dynamic model
result = resolve_sf(model, condition);
condition = result.condition;

derivatives = calculate_sf_derivatives(model, condition);

%% Static Margin approximation

SM = -derivatives.CL_a/derivatives.Cm_a;

%% Forming small disturbance models (Lateral + Longitudinal)

derivatives_sd = calculate_sdsf_derivatives(model, condition);
[Alon, Blon, Alat, Blat] = create_sf_lti(model, condition, derivatives_sd);

%% Performing LTI-analysis to investigate poles and modes

result_longitudinal = sf_lti_analysis(Alon, Blon, 'longitudinal');
result_lateral      = sf_lti_analysis(Alat, Blat, 'lateral');

%% Performing Interactive Root Locus Analysis + Pole placement for SAS


Kq   = 100;
Kp   = -4.7989;
Kr   = 6.9343;
Kyar = 0.6666;

[Kq, Kr, Kp, Kyar] = sf_lti_rlocus(Alon, Blon, Alat, Blat, Kq, Kr, Kp, Kyar);

% Press enter when satisified pole-location.
% Otherwise:
%   Press 1 + Scroll to vary Kp
%   Press 2 + Scroll to vary Kr
%   Press 3 + Scroll to vary Kyar
%   Press 4 + Scroll to vary Kq

%% Performing LTI-analysis to investigate poles and modes

Klat = [0,  0,   Kr, 0;
        0, Kp, Kyar, 0];

Klon = [0, 0, Kq, 0;
        0, 0,  0, 0];


result_longitudinal_ppp = sf_lti_analysis(Alon+Blon*Klon, Blon, 'longitudinal');
result_lateral_ppp      = sf_lti_analysis(Alat+Blat*Klat, Blat, 'lateral');


%% Assigning Control Parameters (PID + Joystick LPF)
% Joystick LP-filters + Gain

%                      _________________
% () stick command --> | s/(s*tau_()+1) |--> Ks_() -->
%                      -----------------
tau_p    = 0.5;
tau_q    = 0.5;  
tau_r    = 0.5;

Ksp      = 0.1;
Ksq      = 0.1;
Ksr      = 0.1; 

% Attitude PID-controllers

% Roll/Pitch controller (Active when not moved manually)
Kp_roll  = -1;
Kd_roll  = -10;
Kp_pitch = -1;
Ki_pitch = -0.1;

% Psi (heading) controller (Active when not moved manually)
Kp_psi = 10;
Kd_psi = 10;
Ki_psi = 0.1;


%% Running Simulation in Simulink

output = sim("Simulation_simulink.slx", reference.tsim);
%%
Xe    = output.Xe;
Vb    = output.Vb;
Wb    = output.Wb;
euler = output.euler;

t     = output.tout;
u     = output.Vb.Data(:,1);
v     = output.Vb.Data(:,2);
w     = output.Vb.Data(:,3);
p     = output.Wb.Data(:,1);
q     = output.Wb.Data(:,2);
r     = output.Wb.Data(:,3);
phi   = output.euler.Data(:,1); 
theta = output.euler.Data(:,2);  
psi   = output.euler.Data(:,3); 
xE    = output.Xe.Data(:,1); 
yE    = output.Xe.Data(:,2); 
zE    = output.Xe.Data(:,3); 

delta_a = output.delta_a.Data(:);
delta_e = output.delta_e.Data(:);
delta_r = output.delta_r.Data(:);
Tcom    = reference.ctrl(4)*ones(size(t));

lat     = output.LLA.Data(:,2);
lon     = output.LLA.Data(:,1);
h       = output.LLA.Data(:,3);

flightdata = array2table([t, delta_e, delta_a, delta_r, u, v, w, p, q, r, phi, theta, psi, xE, yE, zE, lat, lon, h], 'VariableNames', {'t [s]', 'delta_e [rad]', 'delta_a [rad]', 'delta_r [rad]', 'u [m/s]', 'v [m/s]', 'w [m/s]', 'p [rad/s]', 'q [rad/s]', 'r [rad/s]', 'phi [rad]', 'theta [rad]', 'psi [rad]', 'xE [m]', 'yE [m]', 'zE [m]', 'lat [deg]', 'lon [deg]', 'h [m]'});



%% Displaying Simulation Results (Nonlinear 6DOF plant)

max_index = find((t>=max_time), 1, 'first');


% Figure 1: Nonlinear State-Time trajectories
figure_name = append('Input-State trajectories');

% Refine fixed figure size (in pixels)
fig_width = 1600;  
fig_height = 1000;  
screen_size = get(0, 'ScreenSize'); % Get screen size

% Center the figure on screen
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;


% Create fixed-size figure
figure('Name', figure_name, ...
             'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
             'Resize', 'on', ...  % Prevent resizing
             'Toolbar', 'none', 'Menubar', 'figure', ...
             'Color', 'white');



% Subplot 1: Translational velocities (u, v, w) - Nonlinear
subplot(6,1,1);
hold on
plot(t(1:max_index), u(1:max_index), 'color', 'r', 'DisplayName', 'u');
plot(t(1:max_index), v(1:max_index), 'color', 'b', 'DisplayName', 'v');
plot(t(1:max_index), w(1:max_index), 'color', 'g', 'DisplayName', 'w');
hold off
legend('location', 'northwest', 'Box', 'on')
grid minor
xlabel('t [s]')
ylabel('[m/s]')

% Subplot 2: Angular velocities (p, q, r) - Nonlinear
subplot(6,1,2);
hold on
plot(t(1:max_index), 180/pi*p(1:max_index), 'color', 'r', 'DisplayName', 'p');
plot(t(1:max_index), 180/pi*q(1:max_index), 'color', 'b', 'DisplayName', 'q');
plot(t(1:max_index), 180/pi*r(1:max_index), 'color', 'g', 'DisplayName', 'r');
hold off
legend('location', 'northwest', 'Box', 'on')
grid minor
xlabel('t [s]')
ylabel('[deg]')

% Subplot 3: Euler angles (phi, theta, psi) - Nonlinear
subplot(6,1,3);
hold on
plot(t(1:max_index), phi(1:max_index), 'color', 'r', 'DisplayName', 'phi');
plot(t(1:max_index), theta(1:max_index), 'color', 'b', 'DisplayName', 'theta');
plot(t(1:max_index), psi(1:max_index), 'color', 'g', 'DisplayName', 'psi');
plot(signal_in_psi, 'color', 'g', 'LineStyle', '--', 'HandleVisibility', 'off')
hold off
legend('location', 'northwest', 'Box', 'on')
grid minor
xlabel('t [s]')
ylabel('[deg]')

% Subplot 4: Positions (xE, yE, zE) - Nonlinear
subplot(6,1,4);
hold on
plot(t(1:max_index), xE(1:max_index), 'color', 'r', 'DisplayName', 'x_E');
plot(t(1:max_index), yE(1:max_index), 'color', 'b', 'DisplayName', 'y_E');
plot(t(1:max_index), zE(1:max_index), 'color', 'g', 'DisplayName', 'z_E');
hold off
legend('location', 'northwest', 'Box', 'on')
grid minor
xlabel('t [s]')
ylabel('[m]')

% Subplot 5: Control Inputs (delta_e, delta_a, delta_r, Tcom)
subplot(6,1,5);
hold on
plot(t(1:max_index), Tcom(1:max_index), 'color', 'k', 'DisplayName', 'Tcom');
xlabel('t [s]')
ylabel('[%]')
legend('location', 'northwest', 'Box', 'on')
grid minor

subplot(6,1,6);
hold on
plot(t(1:max_index), delta_a(1:max_index), 'color', 'r', 'DisplayName', 'delta_a');
plot(t(1:max_index), delta_e(1:max_index), 'color', 'b', 'DisplayName', 'delta_e');
plot(t(1:max_index), delta_r(1:max_index), 'color', 'g', 'DisplayName', 'delta_r');
hold off
legend('location', 'northwest', 'Box', 'on')
xlabel('t [s]')
ylabel('[deg]')
grid minor

hold off

%% 3D Animation

% Figure 1: Nonlinear State-Time trajectories
figure_name = append('Animation');

% Refine fixed figure size (in pixels)
fig_width = 1300;  
fig_height = 800;  
screen_size = get(0, 'ScreenSize'); % Get screen size

% Center the figure on screen
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;


% Create fixed-size figure
fig = figure('Name', figure_name, ...
             'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
             'Resize', 'on', ...  % Prevent resizing
             'Color', 'white');

hold on

xE =   flightdata.("xE [m]")(1);

yE =   flightdata.("yE [m]")(1);

zE =   flightdata.("zE [m]")(1);

h  = - zE + reference.h;

position = plot3(xE, yE, h, 'Marker','o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'LineWidth',2, 'MarkerSize', 5);

phi   = flightdata.("phi [rad]")(1);
theta = flightdata.("theta [rad]")(1);
psi   = flightdata.("psi [rad]")(1);


R_BV = [cos(theta)*cos(psi), (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)), (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
        cos(theta)*sin(psi), (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)), (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
                -sin(theta),                             sin(phi)*cos(theta),                            cos(phi)*cos(theta)];

R_VE = [1,  0, 0;
        0,  1, 0;
        0,  0,-1];


l = 1000;

x_axis_B = l*[1, 0, 0]';
y_axis_B = l*[0, 1, 0]';
z_axis_B = l*[0, 0, 1]';

x_axis_E = R_VE*R_BV*x_axis_B;
y_axis_E = R_VE*R_BV*y_axis_B;
z_axis_E = R_VE*R_BV*z_axis_B;

x_axis_vector = [xE, xE + x_axis_E(1), yE, yE + x_axis_E(2), h, h + x_axis_E(3)];
y_axis_vector = [xE, xE + y_axis_E(1), yE, yE + y_axis_E(2), h, h + y_axis_E(3)];
z_axis_vector = [xE, xE + z_axis_E(1), yE, yE + z_axis_E(2), h, h + z_axis_E(3)];


x_axis = plot3([x_axis_vector(1), x_axis_vector(2)], [x_axis_vector(3), x_axis_vector(4)], [x_axis_vector(5), x_axis_vector(6)], 'color', 'r');
y_axis = plot3([y_axis_vector(1), y_axis_vector(2)], [y_axis_vector(3), y_axis_vector(4)], [y_axis_vector(5), y_axis_vector(6)], 'color', 'b');
z_axis = plot3([z_axis_vector(1), z_axis_vector(2)], [z_axis_vector(3), z_axis_vector(4)], [z_axis_vector(5), z_axis_vector(6)], 'color', 'g');

trajectory = animatedline('Color', 'k');

ground_track = animatedline('Color', '[0.5, 0.5, 0.5]', 'LineStyle', '--');

view(3)

xlabel('xE [m]')
ylabel('yE [m]')
zlabel('h [m]')

t = 0;

K_clock = 1;

wait_time = dt;

speed = 10;

trel = 0;

tstep = 1;

trel_map = 0;

tstep_map = 1;

kstep = 10;


for k = 1:kstep:max_index

    tic

    xE =   flightdata.("xE [m]")(k);
    yE =   flightdata.("yE [m]")(k);
    h  = - flightdata.("zE [m]")(k) + reference.h;

    u = flightdata.("u [m/s]")(k);
    v = flightdata.("v [m/s]")(k);
    w = flightdata.("w [m/s]")(k);


    phi   = flightdata.("phi [rad]")(k);
    theta = flightdata.("theta [rad]")(k);
    psi   = flightdata.("psi [rad]")(k);

    R_BV = [cos(theta)*cos(psi), (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)), (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
            cos(theta)*sin(psi), (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)), (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
                    -sin(theta),                             sin(phi)*cos(theta),                            cos(phi)*cos(theta)];

    x_axis_E = R_VE*R_BV*x_axis_B;
    y_axis_E = R_VE*R_BV*y_axis_B;
    z_axis_E = R_VE*R_BV*z_axis_B;
    
    x_axis_vector = [xE, xE + x_axis_E(1), yE, yE + x_axis_E(2), h, h + x_axis_E(3)];
    y_axis_vector = [xE, xE + y_axis_E(1), yE, yE + y_axis_E(2), h, h + y_axis_E(3)];
    z_axis_vector = [xE, xE + z_axis_E(1), yE, yE + z_axis_E(2), h, h + z_axis_E(3)];

    
    set(position, 'XData', xE, 'YData', yE, 'ZData', h)
    set(x_axis, 'XData', [x_axis_vector(1), x_axis_vector(2)], 'YData', [x_axis_vector(3), x_axis_vector(4)], 'ZData', [x_axis_vector(5), x_axis_vector(6)])
    set(y_axis, 'XData', [y_axis_vector(1), y_axis_vector(2)], 'YData', [y_axis_vector(3), y_axis_vector(4)], 'ZData', [y_axis_vector(5), y_axis_vector(6)])
    set(z_axis, 'XData', [z_axis_vector(1), z_axis_vector(2)], 'YData', [z_axis_vector(3), z_axis_vector(4)], 'ZData', [z_axis_vector(5), z_axis_vector(6)])
    

    addpoints(trajectory, xE, yE, h);
    addpoints(ground_track, xE, yE, 0);
    
    if trel >= tstep
    
        plot3([x_axis_vector(1), x_axis_vector(2)], [x_axis_vector(3), x_axis_vector(4)], [x_axis_vector(5), x_axis_vector(6)], 'color', [1, 0, 0, 0.5], 'LineWidth', 2);
        plot3([y_axis_vector(1), y_axis_vector(2)], [y_axis_vector(3), y_axis_vector(4)], [y_axis_vector(5), y_axis_vector(6)], 'color', [0, 0, 1, 0.5], 'LineWidth', 2);
        plot3([z_axis_vector(1), z_axis_vector(2)], [z_axis_vector(3), z_axis_vector(4)], [z_axis_vector(5), z_axis_vector(6)], 'color', [0, 1, 0, 0.5], 'LineWidth', 2);

        trel = 0;
    end

    pause(wait_time)


    t = t + dt*kstep;


    loop_time = toc;

    wait_time = wait_time - K_clock*(loop_time - dt*kstep);

    axis equal
    grid on

    trel = trel + dt*kstep;
end


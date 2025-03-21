clear; close all; clc

%% Loading Model YF-16 + Batch with derivatives + Creating Interpolant

load("aircraft\YF-16.mat", "model");

data = readtable("database\batches\YF-161877863644.csv");

load("database\batches\YF-161877863644.mat");

batch.data = data;
batch.info = info;

interpolant = create_interpolant(batch);

%% Solving Straight Flight Conditions (C-L-D) at M = 0.3, h = 1000 m


% Positive gamma sweep
% Initial guess (Tcom, theta, delta_e, alpha)
settings.phase1.x0     =  [1,1,1,10];   
% Logic Control Bit (activate gamma constraint)
settings.phase1.lcb    =  1;   
% Gamma constraint pre sweep
settings.phase1.constr =  0;     
% Nr. of attempts pre sweep;
settings.phase1.tries  =  10;                   

% Negative gamma sweep
% Initial guess (Tcom, theta, delta_e, alpha)
settings.phase2.x0     =  [1,1,1,1];   
% Logic Control Bit (activate gamma constraint)
settings.phase2.lcb    =  1;   
% Gamma constraint pre sweep
settings.phase2.constr =  0;     
% Nr. of attempts pre sweep;
settings.phase2.tries  =  10;                   

solve_sf_conditions(model, batch, settings)

%% Loading a Straight Flight Condition of choice
condition = load_sf_condition(model);           % 1
clc;

%% Creating a reference to initalize simulation capabiltiies + solver

latitude0  = 59.3458;
longitude0 = 18.0731;

autopilot.rollpitch = 1;
autopilot.heading = 1;

dt = 10e-4;
tsim = 100;
N = size(0:dt:tsim-dt);

signal_in_delta_e = timeseries(zeros(N), 0:dt:tsim-dt);
signal_in_delta_a = timeseries(zeros(N), 0:dt:tsim-dt);
signal_in_delta_r = timeseries(zeros(N), 0:dt:tsim-dt);
signal_in_psi     = timeseries(zeros(N), 0:dt:tsim-dt);

reference = create_sf_reference(); 

%% Resolving the solution using interpolant + nonlinear dynamic model
result = resolve_sf(model, condition);
condition = result.condition;

%% Forming small disturbance models (Lateral + Longitudinal)

derivatives_sd = calculate_sdsf_derivatives(model, condition);
[Alon, Blon, Alat, Blat] = create_sf_lti(model, condition, derivatives_sd);

%% Performing LTI-analysis to investigate poles and modes

result_longitudinal = sf_lti_analysis(Alon, Blon, 'longitudinal');
result_lateral      = sf_lti_analysis(Alat, Blat, 'lateral');

%% Performing Interactive Root Locus Analysis + Pole placement for SAS (RESET MODE)

sf_lti_rlocus(Alon, Blon, Alat, Blat, [], [], [], []);

%% Performing Interactive Root Locus Analysis + Pole placement for SAS (RECURSIVE MODE)

Kq   = -100000;
Kp   = -4.7989;
Kr   = 6.9343;
Kyar = 0.6666;

[Kq, Kr, Kp, Kyar] = sf_lti_rlocus(Alon, Blon, Alat, Blat, Kq, Kr, Kp, Kyar);

%% Performing LTI-analysis to investigate poles and modes

Klat = [0,  0,   Kr, 0;
        0, Kp, Kyar, 0];

Klon = [0, 0, Kq, 0;
        0, 0,  0, 0];

%%

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
Kd_roll  =-10;
Kp_pitch = -10;
Ki_pitch = -0.1;

% Psi (heading) controller (Active when not moved manually)
Kp_psi = 10;
Kd_psi = 10;
Ki_psi = 0.1;

%% Simulate in simulink

open("runfg.bat")
open("Simulation_simulink.slx")

%% Simulation settings
dt = 10e-4;

r2d = 180/pi;

%% Running Simulation in MatLab

log_nl = zeros(17, ceil(tsim/dt));

ctrl_ref = reference.ctrl;
x    = reference.x;

t = 0;

try
    while t <= tsim

    ctrl = controller_1(x) + ctrl_ref; % Make controller
    
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

% Formatting Simulation results (states + control inputs) into table

log_nl = array2table(log_nl', 'VariableNames', {'t', ...
                                                'u', 'v', 'w', ...
                                                'p', 'q', 'r', ...
                                                'phi', 'theta', 'psi', ...
                                                'xE', 'yE', 'zE', ...
                                                'delta_e', 'delta_a', 'delta_r', 'Tcom'});


%% Displaying Simulation Results (Nonlinear 6DOF plant)

maxidx = height(log_nl);


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
plot(log_nl.t(1:maxidx), r2d*log_nl.p(1:maxidx), 'color', 'r', 'DisplayName', 'p');
plot(log_nl.t(1:maxidx), r2d*log_nl.q(1:maxidx), 'color', 'b', 'DisplayName', 'q');
plot(log_nl.t(1:maxidx), r2d*log_nl.r(1:maxidx), 'color', 'g', 'DisplayName', 'r');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 3: Euler angles (phi, theta, psi) - Nonlinear
subplot(6,1,3);
hold on
plot(log_nl.t(1:maxidx), r2d*log_nl.phi(1:maxidx), 'color', 'r', 'DisplayName', 'phi');
plot(log_nl.t(1:maxidx), r2d*log_nl.theta(1:maxidx), 'color', 'b', 'DisplayName', 'theta');
plot(log_nl.t(1:maxidx), r2d*log_nl.psi(1:maxidx), 'color', 'g', 'DisplayName', 'psi');
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
    
    ctrl = controller_1(x) + ctrl_ref; % Make controller

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
plot(log_l.t(1:maxidx), r2d*log_l.p(1:maxidx), 'color', 'r', 'DisplayName', 'p');
plot(log_l.t(1:maxidx), r2d*log_l.q(1:maxidx), 'color', 'b', 'DisplayName', 'q');
plot(log_l.t(1:maxidx), r2d*log_l.r(1:maxidx), 'color', 'g', 'DisplayName', 'r');
hold off
legend('location', 'northwest', 'Box', 'on')

% Subplot 3: Euler angles (phi, theta, psi) - Linear
subplot(6,1,3);
hold on
plot(log_l.t(1:maxidx), r2d*log_l.phi(1:maxidx), 'color', 'r', 'DisplayName', 'phi');
plot(log_l.t(1:maxidx), r2d*log_l.theta(1:maxidx), 'color', 'b', 'DisplayName', 'theta');
plot(log_l.t(1:maxidx), r2d*log_l.psi(1:maxidx), 'color', 'g', 'DisplayName', 'psi');
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
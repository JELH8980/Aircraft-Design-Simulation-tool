function conditions = solve_sf_conditions(model, batch, settings)
addpath('aircraft\propulsion\')

T    = batch.data;

info = batch.info;

M = info.M.grid;
h = info.h.grid;

X = [T.alpha, T.delta_e];

% Setting Engine model
engine_model = str2func(strrep(model.parameters.propulsion.method, '.m', ''));  % Convert function name to handle

CL_modelfun  = @(b, x) b(1) + b(2)*x(:,1) + b(3)*sqrt(sin(b(4)*x(:,2)).^2).*x(:,1).^2;

CD_modelfun  = @(b, x) b(1) + b(2)*x(:,1).^2 + b(3)*sqrt(atan(b(4).*x(:,2)).^2).*x(:,1).^2;

Cm_modelfun  = @(b,x) b(1) + b(2)*x(:,1) +b(3)*sqrt(atan(b(4).*x(:,2)).^2).*x(:,1); 


lin_beta0  = [0.1 0.1 0.1 0.1];

quad_beta0 = [0.1 0.1 0.1 0.1];


 % Refine fixed figure size (in pixels)
fig_width = 1200;  
fig_height = 600;  
screen_size = get(0, 'ScreenSize'); % Get screen size

% Center the figure on screen
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;


% Create fixed-size figure
fig = figure('Name', 'Inspection of data', 'Color', 'White', ...
             'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
             'Resize', 'off', ...  % Prevent resizing
             'Color', 'white');

%% Model choice CL (START)

subplot(1,3,1)
hold on
scatter3(T.alpha, T.delta_e, T.CL, 'Marker','.', 'MarkerEdgeColor','k')
xlabel('alpha [deg]')
ylabel('delta_e [deg]')
title('CL(alpha, delta\_e)')

CL_model = fitnlm(X, T.CL, CL_modelfun, lin_beta0);

CL_pred  = predict(CL_model, X);

scatter3(T.alpha, T.delta_e, CL_pred, 'Marker','.', 'MarkerEdgeColor','r')

hold off

view(3)

predictor.CL = CL_model;

subplot(1,3,2)
hold on
scatter3(T.alpha, T.delta_e, T.CD, 'Marker','.', 'MarkerEdgeColor','k')
xlabel('alpha [deg]')
ylabel('delta_e [deg]')
title('CD(alpha, delta\_e)')

% Fitting drag surface

CD_model = fitnlm(X, T.CD, CD_modelfun, quad_beta0);

CD_pred = predict(CD_model, X);

scatter3(T.alpha, T.delta_e, CD_pred, 'Marker','.', 'MarkerEdgeColor','r')

hold off

view(3)

predictor.CD = CD_model;


%% Model choice Cm (START)

subplot(1,3,3)
hold on
scatter3(T.alpha, T.delta_e, T.Cm, 'Marker','.', 'MarkerEdgeColor','k')
xlabel('alpha [deg]')
ylabel('delta_e [deg]')
title('Cm(alpha, delta\_e)')

% Fitting pitch moment surface
Cm_model  = fitnlm(X, T.Cm, Cm_modelfun, lin_beta0);

Cm_pred   = predict(Cm_model, X);

scatter3(T.alpha, T.delta_e, Cm_pred, 'Marker','.', 'MarkerEdgeColor','r')

hold off

view(3)

predictor.Cm = Cm_model;

%% Computational Phase (Uninterrupted)

solutions = [];

xsf = settings.phase1.x0; %[1, 1, 1, 10];     % Tcom, theta, delta_e, alpha   (worth testing different I.C. in terms of alpha depending on low/high mach)

lcb = settings.phase1.lcb;

constr = settings.phase1.constr;           

nr_attempts = settings.phase1.tries;

for attempt = 1:nr_attempts

[xsf, f] = solve_sf_condition(model, predictor, M, h, xsf, lcb, constr);

end


if f < 10e-10

    lcb = 1;
    
    min_gamma = xsf(2) - xsf(4);
    
    max_gamma = xsf(2) - xsf(4);
    
    while f < 10e-10
    
    Tcom      = xsf(1);
    theta     = xsf(2);
    delta_e   = xsf(3);
    alpha     = xsf(4);
    
    gamma = theta - alpha;
    
    gamma = gamma + 1;
    
    constr = gamma;
    
    x0 = xsf;
    
    [xsf, f] = solve_sf_condition(model, predictor, M, h, x0, lcb, constr);
    
    solutions = [solutions; [Tcom, delta_e, alpha, gamma, f]];

    clc;
    disp(append('Increasing gamma, gamma = ', string(gamma), ' [deg]: '));
    disp(' ')
    disp('Solutions found: ')
    disp(' ')
    disp(solutions)
   
    if gamma < min_gamma
        min_gamma = gamma;
    
    elseif gamma > max_gamma
        max_gamma = gamma;
    end
    
    end

end


xsf = settings.phase2.x0; %[1, 1, 1, 10];     % Tcom, theta, delta_e, alpha   (worth testing different I.C. in terms of alpha depending on low/high mach)

lcb = settings.phase2.lcb;

constr = settings.phase2.constr;           

nr_attempts = settings.phase2.tries;

for attempt = 1:nr_attempts

[xsf, f] = solve_sf_condition(model, predictor, M, h, xsf, lcb, constr);

end




if f < 10e-10

    lcb = 1;
    
    min_gamma = xsf(2) - xsf(4);
    
    max_gamma = xsf(2) - xsf(4);
    
    while f < 10e-10
    
    Tcom      = xsf(1);
    theta     = xsf(2);
    delta_e   = xsf(3);
    alpha     = xsf(4);
    
    gamma = theta - alpha;
    
    gamma = gamma - 1;
    
    constr = gamma;
    
    x0 = xsf;
    
    [xsf, f] = solve_sf_condition(model, predictor, M, h, x0, lcb, constr);
    
    solutions = [[Tcom, delta_e, alpha, gamma, f]; solutions];

    clc;
    disp(append('Decreasing gamma, gamma = ', string(gamma), ' [deg]: '));
    disp(' ')
    disp('Solutions found: ')
    disp(' ')
    disp(solutions)
   
    if gamma < min_gamma
        min_gamma = gamma;
    
    elseif gamma > max_gamma
        max_gamma = gamma;
    end
    
    end

end


conditions = array2table(solutions, "VariableNames", {'Tcom [%]', 'delta_e [deg]', 'alpha [deg]', 'gamma [deg]', 'f'});

conditions = conditions(conditions.f < 10e-4, :);

for idx = 1:height(conditions)
    
    condition              = conditions(idx,:);
    
    Tcom                   = condition.("Tcom [%]");
    
    T                      = engine_model(model.parameters.propulsion.Tmax, Tcom, M, h);
    
    delta_e                = pi/180*condition.("delta_e [deg]");
    
    theta                  = pi/180*(condition.("gamma [deg]") + condition.("alpha [deg]"));
    
    alpha                  = pi/180*condition.("alpha [deg]");
    
    [~, ~, rho, ~, ~, a]   = atmosphere_state(h);
    
    V                      = M*a;
    
    u                      = V*cos(alpha);
    
    w                      = tan(alpha).*u;
    
    gamma                  = condition.("gamma [deg]");
    
    
    condition = array2table([M, h, u, w, theta, delta_e, T, Tcom, gamma], 'VariableNames', {'M', 'h', 'u', 'w', 'theta', 'delta_e', 'T', 'Tcom', 'gamma'});
    
    
    save_condition(condition)

end


function [xstar, f] = solve_sf_condition(model, predictor, M, h, x0, lcb, constr)


lcb_gamma = lcb;


c_gamma   = constr;

g = 9.81;

f_CL = @(alpha, delta_e) predict(predictor.CL, [alpha, delta_e]);

f_CD = @(alpha, delta_e) predict(predictor.CD, [alpha, delta_e]);

f_Cm = @(alpha, delta_e) predict(predictor.Cm, [alpha, delta_e]) ;

[~, ~, rho, ~, ~, a] = atmosphere_state(h);

V    = M*a;

Q    = 1/2*rho*V^2;

m    = model.parameters.m;

S    = model.parameters.S;

c    = model.parameters.c;

lB   = model.parameters.lB;

hB   = model.parameters.hB;

lR   = model.parameters.lR;

hR   = model.parameters.hR;

lT   = model.parameters.lT;

hT   = model.parameters.hT;

aT   = model.parameters.aT;

a1 = hB-hR;

a2 = lR-lB;

a3 = cosd(aT)*(hB-hT) + sind(aT)*(lB-lT);

Tmax = model.parameters.propulsion.Tmax;

% Solving Trim conditions

[xstar, f] = fminsearch(@sf, x0, optimset('TolX', 1e-10, 'MaxFunEvals', 10000, 'MaxIter', 10000));

function f = sf(x)

    Tcom      = x(1);
    theta     = x(2);
    delta_e   = x(3);
    alpha     = x(4);

   
    CD = f_CD(alpha, delta_e);

    CL = f_CL(alpha, delta_e);
        
    CX   = -cosd(alpha)*CD + sind(alpha)*CL;

    CZ   = -sind(alpha)*CD - cosd(alpha)*CL;
        
    XA = Q*S*CX;

    ZA = Q*S*CZ;

    Cm = f_Cm(alpha, delta_e);
    
    MA = Q*S*c*Cm;

    T = engine_model(Tmax, Tcom, M, h);

    f1 = XA - m*g*sind(theta) + T*cosd(aT);
    
    f2 = ZA + m*g*cosd(theta) - T*sind(aT);
    
    f3 = MA + a1*XA + a2*ZA + a3*T;

    H = diag([1/(m*g)^2, 1/(m*g)^2, 1/(Q*S*c)^2]);
    f = [f1, f2, f3]*H*[f1; f2; f3];

    if Tcom > 100 || Tcom < 0
        f = f + 1;
    end

    if delta_e > model.parameters.delta_e.max || delta_e < model.parameters.delta_e.min
        f = f + 1;
    end
    
    gamma = theta - alpha;

    f = f + lcb_gamma*(c_gamma-gamma)^2;

end
end


function save_condition(condition)
    
    cd("database\")
    
    if ~isfolder("conditions")
        mkdir("conditions")
    end
    cd("conditions\")
    
    if ~isfolder("SF")
        mkdir("SF")
    end
    cd("SF\")



    if condition.("gamma") < 1 &&  condition.("gamma") > -1

        if ~isfolder("L")
            mkdir("L")
        end

        cd("L\")
        
        filename = strcat("sf_l_", model.name, ".csv");
        
        if isfile(filename)
            existing_data = readtable(filename, 'VariableNamingRule', 'preserve');  % Read existing data
            sf_l_conditions = [existing_data; condition];  % Append new data
            sf_l_conditions = unique(sf_l_conditions, 'rows');

        else
            sf_l_conditions = round(condition, 8); % Adjust precision as needed

        end

        writetable(sf_l_conditions, filename);

    end

    if condition.("gamma") < -1

        if ~isfolder("D")
            mkdir("D")
        end

        cd("D\")
        
        filename = strcat("sf_d_", model.name, ".csv");
        
        if isfile(filename)
            existing_data = readtable(filename, 'VariableNamingRule', 'preserve');  % Read existing data
            sf_d_conditions = [existing_data; condition];  % Append new data
            sf_d_conditions = unique(sf_d_conditions, 'rows');

        else
            sf_d_conditions = round(condition, 8); % Adjust precision as needed

        end

        writetable(sf_d_conditions, filename);

    end


    if condition.("gamma") > 1

        if ~isfolder("C")
            mkdir("C")
        end

        cd("C\")
        
        filename = strcat("sf_c_", model.name, ".csv");
        
        if isfile(filename)
            existing_data = readtable(filename, 'VariableNamingRule', 'preserve');  % Read existing data
            sf_c_conditions = [existing_data; condition];  % Append new data
            sf_c_conditions = unique(sf_c_conditions, 'rows');

        else
            sf_c_conditions = round(condition, 8); % Adjust precision as needed

        end

        writetable(sf_c_conditions, filename);

    end

    
    cd("..\..\..\..\")


end

   

end
function [model] = assign_thrustframe(model)
% ASSIGN_THRUSTFRAME - Function to interactively set the thrust frame position and orientation.
%
% This function allows the user to manually adjust the thrust (T) position
% and angle within the aircraft model. The thrust frame parameters 
% (lT, hT, aT) are visually represented, and the user can interactively 
% modify them through keyboard and mouse inputs.
%
% INPUTS:
%   model - Struct containing the aircraft model with geometry and 
%           parameters fields, including:
%           - model.parameters.lT (Thrust position along X-axis)
%           - model.parameters.hT (Thrust position along Z-axis)
%           - model.parameters.aT (Thrust angle)
%
% OUTPUTS:
%   model - Updated struct with new thrust parameters.
%
% FUNCTIONALITY:
% - Initializes the thrust frame parameters based on model defaults if missing.
% - Plots the aircraft geometry with thrust position markers.
% - Allows interactive adjustment via:
%     - Arrow keys (WASD) to move the thrust position.
%     - Mouse scroll to rotate the thrust angle.
%     - Enter key to confirm selection and close the figure.
% - Ensures real-time visualization of parameter updates.
%
% ERROR HANDLING:
% - The function checks for NaN values in the thrust parameters and assigns
%   default values if necessary.
% - Ensures the model remains within reasonable position limits.
%
% NOTES:
% - The function assumes that the model's geometry (model.geo.mesh) exists.
% - Angle adjustments use degrees, while position changes use meters.
% - The final thrust angle (aT) is negated before returning the model.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


if any(isnan([model.parameters.lT, 0, model.parameters.hT, model.parameters.aT]))
    model.parameters.lT = model.parameters.lR;
    model.parameters.hT = model.parameters.hR;
    model.parameters.aT = 0;
end

arrow_length = 2;

% Changing variables to plot only partition outline

mesh = model.geo.mesh;

% Find the maximum and minimum values along each axis
max_x = max(max(mesh.XYZ(:,:,1)+model.parameters.lR));
min_x = min(min(mesh.XYZ(:,:,1)+model.parameters.lR));
max_y = max(max(mesh.XYZ(:,:,2)));
min_y = min(min(mesh.XYZ(:,:,2)));
max_z = max(max(mesh.XYZ(:,:,3)+model.parameters.hR));
min_z = min(min(mesh.XYZ(:,:,3)+model.parameters.hR));


% Calculate the maximum absolute values for each axis
max_abs_x = max(abs(max_x), abs(min_x));
max_abs_y = max(abs(max_y), abs(min_y));
max_abs_z = max(abs(max_z), abs(min_z));

% Find the largest absolute value across all axes
max_limit = max([max_abs_x, max_abs_y, max_abs_z]);
% Add a 5% margin
margin = 0.1 * max_limit;

figure_name = append('Setting thrust (T) of ', model.name);

% Refine fixed figure size (in pixels)
fig_width = 1000;  
fig_height = 500;  
screen_size = get(0, 'ScreenSize'); % Get screen size

% Center the figure on screen
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;

% Create fixed-size figure
fig = figure('Name', figure_name, ...
             'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
             'Resize', 'off', ...  % Prevent resizing
             'Toolbar', 'none', 'Menubar', 'figure', 'Units', 'normalized', ...
             'Color', 'white');

%% Text boxes

% Information


annotation("textbox", "String", "Thrust (T) ", ... 
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.9, 0.1, 0.1], ...
           "EdgeColor", 'white', 'FontWeight','bold');

XT_info = annotation("textbox", "String", append('lT: ', num2str(round(model.parameters.lT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.84, 0.1, 0.1], ...
           "EdgeColor", 'white');

ZT_info = annotation("textbox", "String", append('hT: ',num2str(round(model.parameters.hT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.76, 0.1, 0.1], ...
           "EdgeColor", 'white');

aT_info = annotation("textbox", "String", append('aT: ', num2str(round(model.parameters.aT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.72, 0.1, 0.1], ...
           "EdgeColor", 'white');

% Instructions

annotation("textbox", "String", " \leftarrow / \rightarrow move (T) left / right", ...
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.1, 0.1, 0.1], ...
           "EdgeColor", 'white');


annotation("textbox", "String", " \uparrow   /  \downarrow  move (T) up / down ", ...
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.05 0.1, 0.1], ...
           "EdgeColor", 'white');


annotation("textbox", "String", " Mouse wheel up / down to increase / decrease \alpha_T", ...
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.3, 0.1 0.1, 0.1], ...
           "EdgeColor", 'white');


%% First view
subplot(1, 2, 1)

fill3(mesh.XYZ(:,:,1)'+model.parameters.lR, mesh.XYZ(:,:,2)', mesh.XYZ(:,:,3)'+model.parameters.hR, 'w', 'LineWidth',  0.5);
axis equal, hold on, view(0, 0)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid minor


plot3(0, 0, 0, 'g+', 'MarkerSize', 10, 'LineWidth', 1);
plot3(0, 0, 0, 'go', 'MarkerSize', 10, 'LineWidth', 1);
plot3(model.parameters.lR, 0, model.parameters.hR, 'r+', 'MarkerSize', 10, 'LineWidth', 1);
plot3(model.parameters.lR, 0, model.parameters.hR, 'ro', 'MarkerSize', 10, 'LineWidth', 1);

line([model.parameters.XM model.parameters.XM + model.parameters.c], [model.parameters.YM model.parameters.YM], [model.parameters.ZM model.parameters.ZM], 'LineWidth', 5, 'Color', 'b');


% Create the black point (marker) and the black arrow (line)
T_1 = plot3(model.parameters.lT, 0, model.parameters.hT, 'ko', 'MarkerSize', 10);
T_a_1 = quiver3(model.parameters.lT, 0, model.parameters.hT, ...
                  -arrow_length * cosd(model.parameters.aT), 0, -arrow_length * sind(model.parameters.aT), ...
                  'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1);

% Set the plot limits for each axis
xlim([min_x - margin, max_x + margin]);
ylim([min_y - margin, max_y + margin]);
zlim([min_z - margin, max_z + margin]);


%% Second view

subplot(1, 2, 2)

g = fill3(mesh.XYZ(:,:,1)'+model.parameters.lR, mesh.XYZ(:,:,2)', mesh.XYZ(:,:,3)'+model.parameters.hR, 'w');
set(g, 'LineWidth', 0.5);
axis equal, hold on, view(0, 90)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid minor


plot3(0, 0, 0, 'g+', 'MarkerSize', 10, 'LineWidth', 1);
plot3(0, 0, 0, 'go', 'MarkerSize', 10, 'LineWidth', 1);
plot3(model.parameters.lR, 0, model.parameters.hR, 'r+', 'MarkerSize', 10, 'LineWidth', 1);
plot3(model.parameters.lR, 0, model.parameters.hR, 'ro', 'MarkerSize', 10, 'LineWidth', 1);

line([model.parameters.XM model.parameters.XM + model.parameters.c], [model.parameters.YM model.parameters.YM], [model.parameters.ZM model.parameters.ZM], 'LineWidth', 5, 'Color', 'b');

% Create the black point (marker) and the black arrow (line)
T_2 = plot3(model.parameters.lT, 0, model.parameters.hT, 'ko', 'MarkerSize', 10);
T_a_2 = quiver3(model.parameters.lT, 0, model.parameters.hT, ...
                 -arrow_length * cos(model.parameters.aT), 0, -arrow_length * sind(model.parameters.aT), ...
                  'Color', 'k', 'LineWidth', 1, 'MaxHeadSize', 1);

% Set the plot limits for each axis
xlim([min_x - margin, max_x + margin]);
ylim([min_y - margin, max_y + margin]);
zlim([min_z - margin, max_z + margin]);


%% Defining callbacks

% Set key press and scroll wheel functions
set(gcf, 'KeyPressFcn', @keyPressCallback);
set(gcf, 'WindowScrollWheelFcn', @scrollWheelCallback);

uiwait(fig);

model.parameters.aT = - model.parameters.aT;


%% Callback functions

function keyPressCallback(~, event)
    % Callback for key press events to move the green point
    step_size = 0.1;
    
    switch event.Key
        case 'w'
            model.parameters.hT = model.parameters.hT + step_size;
        case 's'
            model.parameters.hT = model.parameters.hT - step_size;
        case 'a'
            model.parameters.lT = model.parameters.lT - step_size;
        case 'd'
            model.parameters.lT = model.parameters.lT + step_size;
        case 'return'
            close(fig);
            return;
    end

    set(XT_info, "String", append('lT: ', num2str(round(model.parameters.lT,4))))
    set(ZT_info, "String", append('hT: ', num2str(round(model.parameters.hT,4))))
    

    % Update the green point and the arrow in both subplots
    updatePlot();
end

function scrollWheelCallback(~, event)
    % Callback for scroll wheel events to rotate the arrow in the xz plane
    % Change the angle based on the scroll direction
    rotation_speed = 1;  % Rotation speed in radians per scroll
    model.parameters.aT = model.parameters.aT - event.VerticalScrollCount * rotation_speed;
    set(aT_info, "String", append('aT: ', num2str(round(-model.parameters.aT,4))))
    % Update the green arrow direction in both subplots
    updatePlot();
end

function updatePlot()
    % Update the plot with new green point and arrow position in both subplots
    set(T_1, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT);
    set(T_a_1, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT, ...
                   'UData', -arrow_length * cosd(model.parameters.aT), 'VData', 0, 'WData', arrow_length * -sind(model.parameters.aT));
    
    set(T_2, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT);
    set(T_a_2, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT, ...
                   'UData', -arrow_length * cosd(model.parameters.aT), 'VData', 0, 'WData', arrow_length * -sind(model.parameters.aT));
end



end

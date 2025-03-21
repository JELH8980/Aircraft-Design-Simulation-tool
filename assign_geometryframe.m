function [model] = assign_geometryframe(model)
% ASSIGN_GEOMETRYFRAME - Function to define and visualize the geometry frame of a model.
%
% This function creates an interactive graphical interface to visualize and adjust 
% the geometry frame of a model. It displays the model's mesh, reference points (R), 
% thrust points (T), and mean aerodynamic chord (MAC) in two 3D views. Users can 
% move the geometry frame using arrow keys (W, A, S, D) and confirm changes with 
% the Enter key. The function calculates plot limits with a margin and updates 
% the model parameters dynamically based on user input.
%
% INPUTS:
%   model          - The model struct containing geometry data (model.geo), 
%                    mesh data (model.geo.mesh), and parameters (model.parameters) 
%                    including lR, hR, lT, hT, aT, lB, hB, XM, YM, ZM, and c.
%
% OUTPUT:
%   model          - The updated model with adjusted geometry parameters based on 
%                    user interactions.
%
% FUNCTIONALITY:
% - Extracts the mesh and reference parameters (lR, hR) from the model.
% - Calculates the maximum and minimum bounds of the mesh in X, Y, Z directions 
%   with a 10% margin for visualization.
% - Creates a fixed-size figure with two subplots: a side view (view(0,0)) and a 
%   top view (view(0,90)), displaying the mesh, reference points, thrust points, 
%   and MAC line.
% - Adds text annotations to display reference (R) and thrust (T) coordinates, 
%   along with movement instructions.
% - Implements keypress callbacks (W, A, S, D) to shift the geometry frame and 
%   updates the mesh, points, and MAC line in real-time.
% - Closes the figure and finalizes the model updates when the Enter key is pressed.
%
% NOTES:
% - The figure is centered on the screen with a fixed size (1000x500 pixels) and 
%   resizing disabled.
% - The mesh is plotted with no face color and thin lines for clarity.
% - The function assumes the existence of nested functions `keyPressCallback` and 
%   `update` for interactivity.
% - Units are assumed to be in meters (m) for axes labels.
%
% Author: Ludwig Horvath
% Date: 3/17/2025
mesh = model.geo.mesh;

XR0 = model.parameters.lR;
ZR0 = model.parameters.hR;

% Find the maximum and minimum values along each axis
max_x = max(max(mesh.XYZ(:,:,1)+XR0));
min_x = min(min(mesh.XYZ(:,:,1)+XR0));
max_y = max(max(mesh.XYZ(:,:,2)));
min_y = min(min(mesh.XYZ(:,:,2)));
max_z = max(max(mesh.XYZ(:,:,3)+ZR0));
min_z = min(min(mesh.XYZ(:,:,3)+ZR0));


% Calculate the maximum absolute values for each axis
max_abs_x = max(abs(max_x), abs(min_x));
max_abs_y = max(abs(max_y), abs(min_y));
max_abs_z = max(abs(max_z), abs(min_z));

% Find the largest absolute value across all axes
max_limit = max([max_abs_x, max_abs_y, max_abs_z]);
% Add a 5% margin
margin = 0.1 * max_limit;

figure_name = append('Setting origin (G) of ', model.name);

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

annotation("textbox", "String", "Reference (R) ", ... 
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.9, 0.1, 0.1], ...
           "EdgeColor", 'white', 'FontWeight','bold');

XR_info = annotation("textbox", "String", append('lR: ', num2str(round(model.parameters.lR, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.84, 0.1, 0.1], ...
           "EdgeColor", 'white');

ZR_info = annotation("textbox", "String", append('hR: ',num2str(round(model.parameters.hR, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.76, 0.1, 0.1], ...
           "EdgeColor", 'white');


annotation("textbox", "String", "Thrust (T) ", ... 
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.3, 0.9, 0.1, 0.1], ...
           "EdgeColor", 'white', 'FontWeight','bold');

XT_info = annotation("textbox", "String", append('lT: ', num2str(round(model.parameters.lT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.3, 0.84, 0.1, 0.1], ...
           "EdgeColor", 'white');

ZT_info = annotation("textbox", "String", append('hT: ',num2str(round(model.parameters.hT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.3, 0.76, 0.1, 0.1], ...
           "EdgeColor", 'white');

annotation("textbox", "String", append('aT: ', num2str(round(model.parameters.aT, 4))), ...
           "FontSize", 8, ...
           "FitBoxToText", 'on', ...
           "Position", [0.3, 0.72, 0.1, 0.1], ...
           "EdgeColor", 'white');





% Instructions

annotation("textbox", "String", " \leftarrow / \rightarrow move (G) left / right", ...
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.1, 0.1, 0.1], ...
           "EdgeColor", 'white');


annotation("textbox", "String", " \uparrow   /  \downarrow  move (G) up / down ", ...
           "FontSize", 10, ...
           "FitBoxToText", 'on', ...
           "Position", [0.1, 0.05 0.1, 0.1], ...
           "EdgeColor", 'white');



%% First view
subplot(1, 2, 1)

mesh_1 = fill3(mesh.XYZ(:,:,1)'+model.parameters.lR, mesh.XYZ(:,:,2)'+0, mesh.XYZ(:,:,3)'+model.parameters.hR, 'w');
set(mesh_1, 'FaceColor', 'none', 'LineWidth', 0.1);
axis equal, hold on, view(0, 0)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid minor

G_11 = plot3(model.parameters.lB, 0, model.parameters.hB, 'ko');
set(G_11, 'MarkerSize', 12, 'LineWidth', 1);
G_12 = plot3(model.parameters.lB, 0, model.parameters.hB, 'k+');
set(G_12, 'MarkerSize', 12, 'LineWidth', 1);

T_1 = plot3(model.parameters.lT, 0, model.parameters.hT, 'ko');
set(T_1, 'MarkerSize', 10, 'LineWidth', 1);

R_11 = plot3(model.parameters.lR, 0, model.parameters.hR, 'ro');
set(R_11, 'MarkerSize', 10, 'LineWidth', 1);
R_12 = plot3(model.parameters.lR, 0, model.parameters.hR, 'r+');
set(R_12, 'MarkerSize', 10, 'LineWidth', 1);

MAC_1 = line([model.parameters.XM model.parameters.XM + model.parameters.c], [model.parameters.YM model.parameters.YM], [model.parameters.ZM model.parameters.ZM]);
set(MAC_1, 'LineWidth', 5, 'Color', 'b');

plot3(0, 0, 0, 'g+', 'MarkerSize', 10);
plot3(0, 0, 0, 'go', 'MarkerSize', 10);

% Set the plot limits for each axis
xlim([min_x - margin, max_x + margin]);
ylim([min_y - margin, max_y + margin]);
zlim([min_z - margin, max_z + margin]);


%% Second view

subplot(1, 2, 2)

mesh_2 = fill3(mesh.XYZ(:,:,1)'+model.parameters.lR, mesh.XYZ(:,:,2)', mesh.XYZ(:,:,3)'++model.parameters.hR, 'w');
set(mesh_2, 'FaceColor', 'none', 'LineWidth', 0.1);

axis equal, hold on, view(0, 90)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid minor


G_21 = plot3(model.parameters.lB, 0, model.parameters.hB, 'ko');
set(G_21, 'MarkerSize', 12, 'LineWidth', 1);
G_22 = plot3(model.parameters.lB, 0, model.parameters.hB, 'k+');
set(G_22, 'MarkerSize', 12, 'LineWidth', 1);

T_2 = plot3(model.parameters.lT, 0, model.parameters.hT, 'ko');
set(T_2, 'MarkerSize', 10, 'LineWidth', 1);


R_21 = plot3(model.parameters.lR, 0, model.parameters.hR, 'ro');
set(R_21, 'MarkerSize', 10, 'LineWidth', 1);
R_22 = plot3(model.parameters.lR, 0, model.parameters.hR, 'r+');
set(R_22, 'MarkerSize', 10, 'LineWidth', 1);

MAC_2 = line([model.parameters.XM model.parameters.XM + model.parameters.c], [model.parameters.YM model.parameters.YM], [model.parameters.ZM model.parameters.ZM]);
set(MAC_2, 'LineWidth', 5, 'Color', 'b');


plot3(0, 0, 0, 'g+', 'MarkerSize', 10);
plot3(0, 0, 0, 'go', 'MarkerSize', 10);

% Set the plot limits for each axis
xlim([min_x - margin, max_x + margin]);
ylim([min_y - margin, max_y + margin]);
zlim([min_z - margin, max_z + margin]);


%% Defining callbacks

% Set key press and scroll wheel functions
set(gcf, 'KeyPressFcn', @keyPressCallback);

uiwait(fig);


%% Callback functions

function keyPressCallback(~, event)
    % Callback for key press events to move the green point
    step_size = 0.1;
    
    switch event.Key
        case 'w'
            shift = [0, 0, step_size]; % Move everything down
            update(shift);
        case 's'
            shift = [0, 0, -step_size]; % Move everything up
            update(shift);
        case 'd'
            shift = [step_size, 0, 0]; % Move everything right
            update(shift);
        case 'a'
            shift = [-step_size, 0, 0]; % Move everything left
            update(shift);
        case 'return'
            close(fig);
            return;

    end


end



function update(shift)

    shiftX = shift(1);
    shiftZ = shift(3);

    % Shift all elements (mesh, reference points, MAC line)
    
    % Update mesh positions
    mesh.XYZ(:,:,1) = mesh.XYZ(:,:,1) + shiftX;
    mesh.XYZ(:,:,3) = mesh.XYZ(:,:,3) + shiftZ;

    max_x = max_x + shiftX;
    min_x = min_x + shiftX;
    max_z = max_z + shiftZ;
    min_z = min_z + shiftZ;

    subplot(1,2,1)
    % Set the plot limits for each axis
    xlim([min_x - margin, max_x + margin]);
    ylim([min_y - margin, max_y + margin]);
    zlim([min_z - margin, max_z + margin]);
    
    subplot(1,2,2)
    % Set the plot limits for each axis
    xlim([min_x - margin, max_x + margin]);
    ylim([min_y - margin, max_y + margin]);
    zlim([min_z - margin, max_z + margin]);
        
    % Update mesh plots
    set(mesh_1, 'XData', mesh.XYZ(:,:,1)' + XR0, 'YData', mesh.XYZ(:,:,2)', 'ZData', mesh.XYZ(:,:,3)' + ZR0);
    set(mesh_2, 'XData', mesh.XYZ(:,:,1)' + XR0, 'YData', mesh.XYZ(:,:,2)', 'ZData', mesh.XYZ(:,:,3)' + ZR0);
    
    % Update thrust point
    model.parameters.lT = model.parameters.lT + shiftX;
    model.parameters.hT = model.parameters.hT + shiftZ;
    

    model.parameters.lR = model.parameters.lR + shiftX;
    model.parameters.hR = model.parameters.hR + shiftZ;


    set(T_1, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT);
    set(T_2, 'XData', model.parameters.lT, 'YData', 0, 'ZData', model.parameters.hT);

    % Update reference point
    % model.geo.ref_point = model.geo.ref_point + shift;
    set(R_11, 'XData', model.parameters.lR, 'YData', 0, 'ZData', model.parameters.hR);
    set(R_12, 'XData', model.parameters.lR, 'YData', 0, 'ZData', model.parameters.hR);
    set(R_21, 'XData', model.parameters.lR, 'YData', 0, 'ZData', model.parameters.hR);
    set(R_22, 'XData', model.parameters.lR, 'YData', 0, 'ZData', model.parameters.hR);

    % Update MAC line
    model.parameters.XM = model.parameters.XM + shiftX;
    model.parameters.ZM = model.parameters.ZM + shiftZ;
    set(MAC_1, 'XData', [model.parameters.XM model.parameters.XM + model.parameters.c], 'YData', [model.parameters.YM model.parameters.YM], 'ZData', [model.parameters.ZM model.parameters.ZM]);
    set(MAC_2, 'XData', [model.parameters.XM model.parameters.XM +  model.parameters.c], 'YData', [model.parameters.YM model.parameters.YM], 'ZData', [model.parameters.ZM model.parameters.ZM]);
    

    set(XT_info, "String", append('lT: ', num2str(round(model.parameters.lT,4))))
    set(ZT_info, "String", append('hT: ', num2str(round(model.parameters.hT,4))))

    set(XR_info, "String", append('lR: ', num2str(round(model.parameters.lR,4))))
    set(ZR_info, "String", append('hR: ', num2str(round(model.parameters.hR,4))))
    


end


end
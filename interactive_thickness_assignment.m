function [model] = interactive_thickness_assignment(model, region_name)
% INTERACTIVE_THICKNESS_ASSIGNMENT - Function for interactively modifying thickness and volume of a selected region.
%
% This function allows users to adjust the thickness of a given region (patch) in a 3D model. The user can interactively
% control the thickness by scrolling the mouse wheel and adjust the scaling of the mouse scroll using the up and down 
% arrow keys. The function updates the geometry and visual representation of the region based on the modified thickness.
% The geometry is extruded along the average normal vector of the region, and the resulting volume and thickness are 
% displayed dynamically on the figure. Symmetry is also handled, with the ability to mirror the extrusion on one axis.
%
% INPUTS:
%   model       - A struct containing the model's geometry information. This should include a region (patch) 
%                 specified by `region_name`.
%   region_name - The name of the region (patch) to modify. The region is expected to be found in 
%                 `model.geo.patches` and should have frame information.
%
% OUTPUTS:
%   model       - The updated model struct with modified geometry, volume, and symmetry information for the region.
%
% FUNCTIONALITY:
% - The function creates a fixed-size figure where the user can interactively modify the thickness of the specified 
%   region.
% - The extrusion of the region is controlled via mouse scroll: scrolling up increases thickness, scrolling down decreases it.
% - The function supports symmetric extrusion along the y-axis if `model.geo.patches.(region_name).symmetry` is true.
% - The volume and thickness of the extruded geometry are dynamically displayed in the figure.
% - The user can also adjust the mouse scaling factor using the up and down arrow keys.
% - The figure includes multiple light sources and smooth shading to improve visualization of the extruded geometry.
% - The figure closes when the Enter key is pressed.
%
% NOTES:
% - The function uses `convhull` to compute the convex hull of the extruded points, representing the boundary of the region.
% - The extrusion is performed along the average normal vector of the region's points.
% - The figure is centered on the screen and has a fixed size, with no option to resize it.
% - The function dynamically updates the region's convex hull, volume, and points in `model.geo.bodies` upon each adjustment.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


warning off all

mousescale = 0.0001;
mousescalingparameter = 1;

% Set up the figure

figure_name = append('Extruding ', region_name);

 % Refine fixed figure size (in pixels)
fig_width = 500;  
fig_height = 600;  
screen_size = get(0, 'ScreenSize'); % Get screen size

% Center the figure on screen
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;


% Create fixed-size figure
fig = figure('Name', figure_name, ...
      'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
      'Resize', 'off', ...  % Prevent resizing
      'Toolbar', 'none', 'Menubar', 'figure', 'Color','white');

axis equal

% Add multiple light sources for better illumination
camlight('headlight')  % Light follows camera
camlight left          % Additional light from left
camlight right         % Additional light from right
lighting phong         % Smooth shading
material metal         % Metallic shine



set(fig, 'WindowScrollWheelFcn', @scroll_callback);
set(fig, 'WindowKeyPressFcn', @keypress_callback);

hold on
view(3)

% Create text objects for displaying values
annotation('textbox', [0.1, 0.85, 0.3, 0.05], ...
                                      'String', 'Thickness: ', ...
                                      'FontSize', 10, ...
                                      'FontWeight', 'normal', ...
                                      'EdgeColor', 'k', 'Interpreter','latex');

annotation('textbox', [0.6, 0.85, 0.3, 0.05], ...
                                   'String', 'Volume: ', ...
                                   'FontSize', 10, ...
                                   'FontWeight', 'normal', ...
                                   'EdgeColor', 'k', 'Interpreter','latex');

thick_annotation_value = annotation('textbox', [0.1+0.15, 0.85, 0.3, 0.05], ...
                                      'String', '', ...
                                      'FontSize', 10, ...
                                      'FontWeight', 'normal', ...
                                      'EdgeColor', 'none', 'Interpreter','latex');

vol_annotation_value = annotation('textbox', [0.6+0.12, 0.85, 0.3, 0.05], ...
                                   'String', '', ...
                                   'FontSize', 10, ...
                                   'FontWeight', 'normal', ...
                                   'EdgeColor', 'none', 'Interpreter','latex');


patch = model.geo.patches.(region_name).frames;
points = reshape(patch, [], 3);

avg_normal = calc_avg_normal(points);

function extrude_patch()
    delete(findall(fig, 'Type', 'Patch')); % Delete old shapes
    % Remove existing lights before adding new ones
    delete(findall(gca, 'Type', 'light'))  

    thickness = mousescale; % Store updated thickness
    new_points = [points - avg_normal * thickness; points + avg_normal * thickness]; 

    [K, V] = convhull(new_points);
    trisurf(K, new_points(:,1), new_points(:,2), new_points(:,3), ...
              'FaceColor', [0.7, 0.7, 0.7], ... % Slightly different grey
              'EdgeColor', 'none', ...
              'FaceLighting', 'gouraud', ...
              'SpecularStrength', 1.0, ...
              'SpecularExponent', 50, ...
              'AmbientStrength', 0.5);

    if model.geo.patches.(region_name).symmetry
        new_points_anti = new_points;
        new_points_anti(:,2) = -new_points_anti(:,2);
        [Kanti, ~] = convhull(new_points_anti);
        trisurf(Kanti, new_points_anti(:,1), new_points_anti(:,2), new_points_anti(:,3), ...
                     'FaceColor', [0.7, 0.7, 0.7], ... % Slightly different grey
                     'EdgeColor', 'none', ...
                     'FaceLighting', 'gouraud', ...
                     'SpecularStrength', 1.0, ...
                     'SpecularExponent', 50, ...
                     'AmbientStrength', 0.5);

        symmetry = 1;
    else
        symmetry = 0;
    end

    axis equal

    % Add multiple light sources for better illumination
    camlight('headlight')  % Light follows camera
    camlight left          % Additional light from left
    camlight right         % Additional light from right
    lighting phong         % Smooth shading
    material metal         % Metallic shine

    convexhull = K;
    volume = V;
    
    % Update `geo`
    model.geo.bodies.(region_name).convexhull = convexhull;
    model.geo.bodies.(region_name).volume = volume;
    model.geo.bodies.(region_name).points = new_points;
    model.geo.bodies.(region_name).symmetry = symmetry;
    
    % Update annotation text (Fixed in the grey area)
    set(thick_annotation_value, 'String', string(thickness));
    set(vol_annotation_value, 'String', string(volume));


end

function scroll_callback(~, event)
    step = mousescalingparameter;
    if event.VerticalScrollCount > 0
        mousescale = max(0.01, mousescale - step);
    else
        mousescale = mousescale + step;
    end
    extrude_patch();
end

function avg_normal = calc_avg_normal(points)
    % Initialize normal accumulator
    normal_sum = [0, 0, 0];
    
    % Compute normal using cross product of consecutive triplets
    for i = 1:size(points, 1)-2
        p1 = points(i, :);
        p2 = points(i+1, :);
        p3 = points(i+2, :);

        v1 = p2 - p1;
        v2 = p3 - p1;

        normal = cross(v1, v2);
        normal_sum = normal_sum + normal;
    end
    
    avg_normal = normal_sum / norm(normal_sum);
end

function keypress_callback(~, event)
    if strcmp(event.Key, 'return') % Enter key pressed
        close(fig);
    elseif strcmp(event.Key, 'uparrow')
        mousescalingparameter = min(1, mousescalingparameter * 10); % Increase scale
    elseif strcmp(event.Key, 'downarrow')
        mousescalingparameter = max(0.0001, mousescalingparameter / 10); % Decrease scale
    end
end

extrude_patch(); % Initial extrusion    

% Wait until the figure is closed
waitfor(fig);
warning on all
end

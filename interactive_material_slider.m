function fractions = interactive_material_slider(materials)
% INTERACTIVE_MATERIAL_SLIDER - Function to interactively adjust material fractions.
%
% This function allows the user to adjust the fractions of different materials in a model using a graphical 
% interface. The user can navigate between materials using the left and right arrow keys, and increase or 
% decrease the fraction of the selected material using the up and down arrow keys. The function ensures that the 
% total sum of all material fractions always remains equal to 1, redistributing the fractions of other materials 
% as needed. The selected material is highlighted, and the material colors are shown in a bar plot.
%
% INPUTS:
%   materials    - A struct containing material information. Each field corresponds to a material and includes 
%                  at least a 'color' field with the color of the material (as an RGB string).
%
% OUTPUTS:
%   fractions    - A struct with the updated fractions for each material.
%
% FUNCTIONALITY:
% - The user can interactively adjust the material fractions using keyboard inputs.
% - The 'leftarrow' and 'rightarrow' keys are used to select the previous or next material.
% - The 'uparrow' and 'downarrow' keys are used to increase or decrease the fraction of the selected material.
% - The total fraction of all materials is kept at 1 by redistributing the remaining fractions.
% - The function updates the bar plot in real-time, with each material represented by a colored bar.
% - The selected material is highlighted in the plot.
%
% NOTES:
% - The function uses the `materials` struct, where each material is represented by a field containing a `color` 
%   field (a string in RGB format).
% - The plot is dynamically updated to reflect changes in material fractions.
% - The function ends when the user presses the 'return' key.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

% Extract material names
matNames = fieldnames(materials);
numMats = numel(matNames);

% Initialize fractions (equal distribution)
fractions = struct();
for i = 1:numMats
    fractions.(matNames{i}) = 1 / numMats;
end

% Current material index
currentMatIdx = 1;

% Create figure
fig = figure('KeyPressFcn', @keyPressCallback, 'Name', 'Material Selector', ...
             'NumberTitle', 'off', 'Position', [100, 100, 500, 300]);

% Display initial state
updatePlot();

% Wait until figure is closed
uiwait(fig);

% Key press callback
function keyPressCallback(~, event)
    switch event.Key
        case 'leftarrow' % Select previous material
            currentMatIdx = mod(currentMatIdx - 2, numMats) + 1;
        case 'rightarrow' % Select next material
            currentMatIdx = mod(currentMatIdx, numMats) + 1;
        case 'uparrow' % Increase fraction
            adjustFraction(0.01);
        case 'downarrow' % Decrease fraction
            adjustFraction(-0.01);
        case 'return' % Exit
            close(fig);
    end
    updatePlot();
end

% Adjust fraction while keeping total sum = 1
function adjustFraction(delta)
    matName = matNames{currentMatIdx};
    newFraction = fractions.(matName) + delta;
    
    % Ensure new fraction is within [0,1]
    newFraction = max(0, min(1, newFraction));

    % Compute available space in other materials
    totalOther = sum(struct2array(fractions)) - fractions.(matName);
    if totalOther == 0, return; end % No redistribution possible

    % Scale other materials to maintain total sum = 1
    scaleFactor = (1 - newFraction) / totalOther; % ✅ FIXED
    for i = 1:numMats
        if i ~= currentMatIdx
            fractions.(matNames{i}) = fractions.(matNames{i}) * scaleFactor; % ✅ FIXED
        end
    end
    
    % Apply the new fraction
    fractions.(matName) = newFraction;
end

% Update plot
function updatePlot()
    clf; % Clear figure
    values = struct2array(fractions);
    
    % Create bar plot
    bars = bar(values, 'FaceAlpha', 0.5); % Initialize bars
    
    % Set the face colors
    bars.FaceColor = 'flat'; % Allow individual bar coloring
    
    material_names = fields(materials);
    
    % Assign colors to each bar
    for i = 1:numMats
        color_str = materials.(material_names{i}).color; % Get color as string
        color_rgb = sscanf(color_str, '[%f, %f, %f]')' / 255; % Convert to numeric array & normalize
        bars.CData(i, :) = color_rgb; % Assign to bar
    end
    
    xticklabels(matNames);
    xtickangle(45);
    ylim([0, 1]);
    title(sprintf('Selected: %s', matNames{currentMatIdx}), 'Interpreter', 'latex');
    
    % Highlight selected material
    hold on;
    highlight_bar = bar(currentMatIdx, values(currentMatIdx), 'FaceAlpha', 1);
    highlight_bar.FaceColor = sscanf(materials.(material_names{currentMatIdx}).color, '[%f, %f, %f]')' / 255;
    hold off;

    grid on;
    drawnow;
end
end

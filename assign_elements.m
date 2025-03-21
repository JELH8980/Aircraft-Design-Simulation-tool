function [model] = assign_elements(model)
% ASSIGN_INERTIA - Function to calculate and assign inertia properties to a model.
%
% This function provides an interactive menu for the user to calculate inertia 
% properties of a model using different methods. The user is presented with a series 
% of questions to select a calculation method (e.g., calculate_inertia_1, 
% calculate_inertia_2, or calculate_inertia_3). Based on the selection, the 
% corresponding inertia calculation is applied, and the updated model is returned.
%
% INPUTS:
%   model          - The model struct containing data required for inertia calculations.
%
% OUTPUT:
%   model          - The updated model with calculated inertia properties.
%
% FUNCTIONALITY:
% - Uses a while loop with a switch-case structure to handle user input from a 
%   question prompt (via the `questions` function).
% - Supports three inertia calculation methods (calculate_inertia_1, 
%   calculate_inertia_2, calculate_inertia_3) based on user choice.
% - Allows the user to return to the main menu by selecting option 0.
% - Includes error handling to save the model state in case of an exception using 
%   the `emergency_save` function.
%
% NOTES:
% - If the `questions` function returns an empty response, the loop defaults to an 
%   invalid case (-1), which is handled gracefully.
% - The specific implementation details of the inertia calculation methods 
%   (calculate_inertia_1, etc.) are assumed to be defined elsewhere.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

try
    while true
    q=questions(9);
    
    if isempty(q)
       q=-1; %will go to otherwise section
    end
    
    
    switch (q)
        case 1
            
            options = {};  % Initialize an empty cell array
    
            bodies_names = fields(model.geo.bodies);
            no_bodies = numel(bodies_names);
    
        
            for i = 1:no_bodies
                options = [options, bodies_names{i}];  % Append to the cell array
                disp(append('    [', string(i), ']. ', bodies_names{i}));        
            end
    
            choice = input('    Please enter choice from above: ');
    
            points = model.geo.bodies.(options{choice}).points;
    
            nx = input('    Input no. partitions in x (nx): ');
            ny = input('    Input no. partitions in y (ny): ');
            nz = input('    Input no. partitions in z (nz): ');

            elements = voxelize_convex_hull(points, nx, ny, nz);
            

            model.geo.bodies.(bodies_names{choice}).elements = elements;
    
        case 0
            break;
    
    end

    end

catch
    % Handle errors here
    emergency_save(model)
end

end
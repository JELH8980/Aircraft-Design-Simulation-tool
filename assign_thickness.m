function [model] = assign_thickness(model)
% ========================================================================
% Script Name: assign_thickness.m
% Author: Ludwig Horvath
% Date: 2/11/2025
% Software: Redspot
% Description: 
%   This function allows the user to assign thickness values to specific 
%   patches in the aircraft model's geometry. It provides an interactive 
%   selection process for choosing regions and defining thickness values.
%
%   The function operates as follows:
%     - Lists available geometry patches from model.geo.patches
%     - Prompts the user to select a patch and assign a thickness value
%     - Updates the model structure with the assigned thickness
%     - Ensures error handling with an emergency save mechanism
%
% Notes:
%   - The function assumes the existence of `interactive_thickness_assignment`.
%   - It requires a valid model with `geo.patches` already defined.
% ========================================================================


try
     regions = fields(model.geo.patches);

    while true
        q=questions(7);
    
    if isempty(q)
       q=-1; %will go to otherwise section
    end
    
    switch (q)
        case 1
            options = {};  % Initialize an empty cell array

            for i = 1:numel(regions)
                options = [options, regions{i}];  % Append to the cell array
                disp(append('    [', string(i), ']. ', regions{i}));        
            end
    
            choice = input('    Please enter choice from above: ');

            [model] = interactive_thickness_assignment(model, options{choice});   
    
        case 0
            break;
    
        otherwise
    
    end
    
    end

catch
    % Define Error handling here
    emergency_save(model)
end
    
end
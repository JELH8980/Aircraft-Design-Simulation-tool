function [model] = inptinertia(model)
% INPUTINERTIA - Function to input and handle inertia-related parameters of a model.
%
% This function allows the user to interactively input various inertia parameters for the model, 
% such as the geometry of patches, thickness of regions, material properties, and inertia properties. 
% The user can also display the current inertia results. If an error occurs during the process, 
% the model is automatically saved to avoid data loss. The function continues to prompt the user until 
% the 'exit' option is selected.
%
% INPUTS:
%   model      - The model struct to be updated with inertia-related parameters.
%
% OUTPUTS:
%   model      - The updated model struct after inertia-related parameters have been assigned or displayed.
%
% FUNCTIONALITY:
% - The function continuously prompts the user with a list of questions until the 'exit' option is selected.
% - The user can choose to input patches, assign thickness and materials, assign elements, 
%   or display the inertia results.
% - Each of these steps is performed through corresponding functions (e.g., `assign_patches`, `assign_thickness`, 
%   `assign_material`, `assign_elements`, and `assign_inertia`).
% - The input process is protected with a try-catch block, where any errors encountered will result in an 
%   emergency save of the model to prevent data loss.
%
% NOTES:
% - The questions are defined within the `questions` function (with 4 options).
% - If an error occurs during the input process, the model is saved to the `emergencysave` directory using the `emergency_save` function.
% - If the user chooses option 0, the function breaks the loop and exits.
%
% Author: Ludwig Horvath
% Date: 2/11/2025



try
    while true
    q=questions(4);
    
    if isempty(q)
       q=-1; %will go to otherwise section
    end
    
    switch (q)
        
    
        case 1 % Enter patches                        (create area)
            [model] = assign_patches(model);
    
        case 2 % Enter thickness of different regions (create volume)
            [model] = assign_thickness(model);
       
        case 3 % Enter material of different regions  (create mass) 
            [model] = assign_material(model);
    
        case 4 % Assign elements
            [model] = assign_elements(model);

        case 5 % Assign inertia
            [model] = assign_inertia(model);
    
        case 6 % Display result
            display_inertia(model);
            
        case 0  % Return to main menu
            break;
    
    end
    
    end

catch
    % Handle errors here
    emergency_save(model);

end

end
function [model] = edit_redspot_geo(model)
% EDIT_REDSPOT_GEO - Function to interactively edit the geometry of a model.
%
% This function provides an interactive menu loop for editing a model’s geometry, 
% allowing the user to assign or modify the geometry frame and thrust frame. It 
% uses a question-based interface (via the `questions` function) to prompt user 
% selections and updates the model struct accordingly until the user chooses to exit.
%
% INPUTS:
%   model          - A struct containing the model data to be edited, including 
%                    geometry-related fields (e.g., model.geo, model.parameters).
%
% OUTPUT:
%   model          - The updated model struct with modified geometry frames based 
%                    on user selections.
%
% FUNCTIONALITY:
% - Initializes a loop that continues until the user opts to exit.
% - Displays a menu via `questions(14)` with at least three options:
%   - [1]: Assign or edit the geometry frame using `assign_geometryframe`.
%   - [2]: Assign or edit the thrust frame using `assign_thrustframe`.
%   - [0]: Exit the loop and return the updated model.
% - Executes the chosen action by calling the respective function and updates the 
%   model struct.
% - Includes basic error handling with a try-catch block to catch exceptions 
%   during execution (though no specific error actions are defined).
%
% NOTES:
% - Assumes the existence of external functions: `questions` for menu display, 
%   `assign_geometryframe` for geometry frame editing

try
	loop=1;  
   while loop==1   
      	disp(' ')


        choice=questions(14);
        
      
        switch choice

            case 1
                [model] = assign_geometryframe(model);

            case 2
                [model] = assign_thrustframe(model);

            case 0
                loop = 0;

        end


   end




catch ME
    disp(ME)

end


end

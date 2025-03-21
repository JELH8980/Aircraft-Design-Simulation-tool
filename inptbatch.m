function model = inptbatch(model)
% ASSIGN_STATE_ENVELOPE - Function to assign state envelope data to the model's dynamics.
%
% This function allows the user to define and assign state envelope data for various system states
% within the model's dynamics. The user can select a state variable (e.g., 'alpha', 'beta', 'P', etc.)
% and define its envelope grid values, including the maximum and minimum values for each state variable.
% The user can also visualize the current state envelope through an interactive grid assignment.
%
% INPUTS:
%   model          - The model struct containing the dynamics (model.dynamics), where the state 
%                    envelope data will be assigned. The model must already contain the `dynamics` 
%                    field.
%
% OUTPUTS:
%   model          - The updated model struct with the assigned state envelope data in 
%                    `model.dynamics.envelope`.
%
% FUNCTIONALITY:
% - If the `envelope` field is not already present in `model.dynamics`, it is initialized as an empty 
%   structure.
% - The user is presented with a list of state variables to choose from, including parameters like 
%   'alpha', 'beta', 'P', 'Q', 'R', etc.
% - After selecting a state variable, the user is prompted to define its envelope grid, which represents
%   possible values for that state variable.
% - The grid for each selected state is stored in `model.dynamics.envelope.<field_name>.grid` along with
%   additional information like the number of grid points, maximum, and minimum values.
% - The user can also view the current state envelope via an option to display the envelope for the model.
%
% NOTES:
% - This function assumes that the `interactive_grid_assignment` function will be used to create the
%   grid for each state variable.
% - The envelope for each state is defined by a set of grid points, and the maximum, minimum, and 
%   number of grid points are stored for later use.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

batch.data =  [];
batch.info =  [];

% try
    
    while true
       answ=questions(13);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)

           case 0
               
               break;

           case 1

               batch = create_edit_batch(batch, model);

           case 2

               batch = load_batch(model);

           case 3 
               display_batch(batch);

           case 4
               solve_batch(batch, model);
              
           otherwise
           

       end

    end

% catch
% 
% end

end



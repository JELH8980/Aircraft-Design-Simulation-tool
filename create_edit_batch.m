function batch = create_edit_batch(batch, model)
% CREATE_EDIT_BATCH - Function to create or edit a batch structure interactively.
%
% This function provides an interactive menu for creating or editing a batch 
% structure, allowing the user to assign values to aerodynamic and control state 
% variables. It initializes an empty batch if none exists, displays a list of 
% state labels (e.g., alpha, beta, P, M, delta_a), and prompts the user to input 
% values for each field. The batch can be saved with a unique ID once all required 
% fields are assigned, or canceled to discard changes.
%
% INPUTS:
%   batch          - A struct containing batch information (batch.info) and data 
%                    (batch.data). If batch.info is empty, it is initialized.
%   model          - A model identifier (type not specified) used when assigning 
%                    a batch ID.
%
% OUTPUT:
%   batch          - The updated batch struct with user-assigned fields in 
%                    batch.info and optionally saved data.
%
% FUNCTIONALITY:
% - Defines state labels (with units) and clean field names for 12 variables: 
%   alpha, alpha_dot, beta, beta_dot, P, Q, R, M, h, delta_a, delta_e, delta_r.
% - Initializes batch.info and batch.data if empty.
% - Presents an interactive menu loop:
%   - Displays numbered state labels for user selection.
%   - Option [0] allows exiting to a save/cancel submenu.
%   - Options [1–12] prompt value assignment via `interactive_batch_assignment`.
% - On exit (choice 0):
%   - Checks if all required fields are assigned.
%   - If complete and user selects save (1), assigns a batch ID, sets progress to 1, 
%     and saves the batch using `save_batch`.
%   - If incomplete or canceled (0), clears batch.info or notifies the user.
% - Includes basic error handling to catch and report exceptions.
%
% NOTES:
% - Assumes the existence of external functions: `interactive_batch_assignment` 
%   for field value input, `assign_batch_id` for generating a unique ID, and 
%   `save_batch` for persisting the batch data.
% - Units are specified in state labels: degrees (deg), degrees per second (deg/s), 
%   meters (m), or dimensionless (-).
% - Clears the command window (clc) repeatedly for a clean interface.
% - Does not modify batch.data directly; its purpose is unspecified in this context.
%
% Author: Ludwig Horvath
% Date: 3/17/2025


state_labels = {'alpha [deg]', ...
                'alpha_dot [deg/s]', ...
                 'beta [deg]',  ...
                 'beta_dot [deg/s]', ...
                  'P [deg/s]', ...
                  'Q [deg/s]', ...
                  'R [deg/s]', ...
                      'M [-]', ...
                      'h [m]', ...
              'delta_a [deg]', ...
              'delta_e [deg]', ...
              'delta_r [deg]'};

state_labels_clean = {'alpha', ...
                      'alpha_dot', ...
                      'beta',  ...
                      'beta_dot',  ...
                      'P', ...
                      'Q', ...
                      'R', ...
                      'M', ...
                      'h', ...
                      'delta_a', ...
                      'delta_e', ...
                      'delta_r'};

try
    if isempty(batch.info)
        if isempty(batch.info)
            batch.info = struct();
            batch.data = [];
        end
        
        while true
            clc;
            for i = 1:numel(state_labels)
    
                if i < 10
                    digit_space = ' ';
                else
                    digit_space = '';
                end
    
                disp(append(digit_space, '    [', string(i), ']. ', state_labels{i}));        
            end
    
            disp(' ')
            disp('    [0]. Back / up menu')
            disp(' ')
    
            choice = input('    Please enter choice from above: ');
            disp(' ')
    
            if choice == 0
                answ = 0;
            else
                answ = 1;
            end
    
            if isempty(answ)
                answ = -1;
            end
    
            switch (answ)
    
                case 0 
                    disp(' ')
                    disp('    [1]. Save batch & Back / Up Menu')
                    disp(' ')
                    disp('    [0]. Cancel & Back / Up Menu')
                    disp(' ')
    
                    choice = input('    Please enter choice from above: ');
                    disp(' ')
    
                    % Check if all fields in state_labels exist in batch_info
                    missing_fields = setdiff(state_labels_clean, fieldnames(batch.info));
    
                    if all([isempty(missing_fields), choice == 1])
                        batch.info = assign_batch_id(batch.info, model);
                        batch.info.progress = 1;
                        save_batch(batch);
                        clc;
                        break;
    
                    elseif all([~isempty(missing_fields), choice == 1])
                        disp(' ')
                        disp('    Error: All fields must be assigned before saving.')
                        disp(' ')
                        disp('    Missing fields:')
                        disp(missing_fields)
                        disp(' ')
                        input('    Press Enter to continue...');
                        disp(' ')
                    else
                        batch.info = [];
                        clc;
                        break;
                    end
    
                case 1
                    field_info = split(state_labels{choice}, " ");
                    
    
                    batch.info.(field_info{1}) = interactive_batch_assignment(field_info);   
    
            end % Switch
    
        end % While
    end

catch
    disp('Error in "create_edit_batch.m"')
end % Try

end % Function



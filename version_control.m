function version_control(model)
% VERSION_CONTROL - Function to manage version control and backup operations for a model.
%
% This function provides the user with options to either save a backup of the current model, overwrite an existing model, 
% or return to the previous menu. It is designed to facilitate version management of the model by allowing for backup saving 
% and maintaining comments about changes made to the model.
%
% INPUTS:
%   model      - The current model object that contains the latest version of the model data, including metadata such as the 
%                model name, save comments, and the date of the latest save.
%
% FUNCTIONALITY:
% - The function first checks if a file with the current model name (model.name) exists in the 'aircraft' directory.
% - If the file exists, the user is presented with a menu with the following options:
%   1. Save a backup of the current model.
%   2. Overwrite the current model with the new data, including adding a user-provided comment.
%   3. Return to the previous menu.
% - If the file does not exist, the model is simply saved as a new file.
% - If any errors occur during the process, the `emergency_save` function is called to handle the error and save the model.
%
% OUTPUTS:
%   None. The function performs actions on the filesystem, saving or backing up the model as appropriate.
%
% EXAMPLES:
% - Calling `version_control(model)` with a model that already has an associated `.mat` file will prompt the user to select 
%   an option for version control (backup, overwrite, or return).
% - If the model does not have an associated `.mat` file, it will be saved to the directory.
%
% NOTES:
% - The user is expected to provide a valid integer input (1, 2, or 0) for the menu options. If an invalid option is entered, 
%   the function prompts the user again.
% - The function uses the `datetime` function to track the date of the latest save and stores the comment provided by the user.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


try
    
    disp(' ')
    space = '    ';
    
    cd("aircraft\")
    
    filename = append(model.name, '.mat');
    
    done = false;
    
    while ~done
        clc;
        disp(' ')
        if exist(filename, 'file')
    
            disp(append(space, '[1] Save backup of ', filename))
            disp(append(space, '[2] Overwrite ', filename))
            disp(append(space, '[0] Return'))
           
            answer = input(append(space, 'Select an option: '));
    
            if answer == 1
                
                old_model = load(filename);
    
                cd("backups\")
    
                save(filename, 'old_model', 'old_model')  % This will overwrite 'data.mat' if it exists.
    
                cd("..\")
    
                disp(append(space, 'Backup saved...'))
    
            elseif answer == 2
                comment = input(append(space, 'Provide a comment: '), "s");
    
                model.latestsave.comment = comment;
                model.latestsave.date = datetime;
    
                save(filename, 'model', 'model')  % This will overwrite 'data.mat' if it exists.
    
                
            elseif answer == 0
                break;
    
            else 
                disp('Please provide a valid option..')
    
            end
    
        else
            save(filename, 'model', 'model')
    
        end
    
    end
    
    cd('..\')

catch
    % Handle errors here
    emergency_save(model)

end

end
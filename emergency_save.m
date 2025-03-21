function emergency_save(model)
% EMERGENCY_SAVE - Function to perform an emergency save of a model's data.
%
% This function saves the model's data to a file in a specific directory, with a timestamp and the model's 
% name as part of the filename. The data is saved in the `.mat` format, and the directory is automatically 
% created if necessary. The function ensures that the model data is saved with a timestamp, which prevents 
% overwriting of previous saves and provides an easy way to track different versions of the model.
%
% INPUTS:
%   model          - The model struct, which must include the relevant data to be saved. The function 
%                    checks if the model name is non-empty before attempting to save.
%
% OUTPUTS:
%   None (the function performs a file save operation).
%
% FUNCTIONALITY:
% - The function changes the current directory to `emergencysave\`, ensuring that the save occurs in the 
%   correct location.
% - It generates a timestamped filename by combining the model's name and the current date and time.
% - The model's data is saved using the `save` command with the model struct as input.
% - The function returns to the original directory after saving the file.
%
% NOTES:
% - If the model's name is empty, the function does nothing.
% - The generated filename replaces spaces in the model name with underscores and uses a timestamp format 
%   that is valid for filenames.
% - The function performs the save in the `emergencysave\` directory, which should exist or be created.
%
% Author: Ludwig Horvath
% Date: 2/11/2025



if ~isempty(model.name)
    cd('emergencysave\');
    
    time_stamp = regexprep(string(datetime("now", "Format", "dd-MMM-yyyy_HH-mm-ss")), ':', '-');
    
    filename = append(model.name, '_at_', time_stamp, '.mat');

    filename = regexprep(filename, '\s+', '_');
    
    save(filename, '-struct', 'model');
    
    cd('..')
end

end
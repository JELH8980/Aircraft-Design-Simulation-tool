function [model] = load_model()
% LOAD_MODEL - Function to load a model file and check/upgrade its format.
%
% This function loads a `.mat` file containing a model from the `aircraft` directory. It first lists all available `.mat` files
% and prompts the user to select one. If the selected file's model format is outdated, the user is asked whether to upgrade
% it to a newer format. The function attempts to ensure that the loaded model includes certain fields (e.g., `model.geo.version`, 
% `model.geo.meshtype`). If these fields are missing, the function adds default values and prompts the user to save the model
% in the new format.
%
% INPUTS:
%   None (user input is prompted for the file selection).
%
% OUTPUTS:
%   model       - The loaded model structure, possibly updated to the latest format if the model's format was outdated.
%
% FUNCTIONALITY:
% - The function lists all `.mat` files in the `aircraft` directory and presents them to the user.
% - It then prompts the user to select a file to load.
% - If any necessary fields (such as `model.geo.version` or `model.geo.meshtype`) are missing, the function assigns default values
%   to these fields.
% - If the model file format is outdated, the function asks the user if they want to save the file in the new format.
% - If the user opts to upgrade, the model is saved in the same location with the updated format.
% - If the user declines to upgrade, no changes are made to the file.
%
% NOTES:
% - The function assumes that the `aircraft` directory exists and contains `.mat` files representing models.
% - If the model's format is outdated, the function will only proceed to save the model after user confirmation.
% - The function handles errors gracefully by displaying an error message if something goes wrong during the process.
%
% Author: Ludwig Horvath
% Date: 2/11/2025



try
    clc;
    questions(18);

    cd('aircraft\')
    files = cellstr(ls); % List all files and folders
    files = files(endsWith(files, '.mat'));
    
    
    no_files = numel(files);

    % Display filtered file options
    for file_nr = 1:no_files
        disp(append('    [', string(file_nr),']. ', files{file_nr}))
    end

    disp(' ')
    option = input(' 	Please enter choice from above: ');

    selected_file = files{option};

    load(selected_file, 'model');  


    cd('..')
    
    resave=0;
    
    try model.geo.version;
    
    catch
       resave=1;
       model.geo.name=('Undefined');
       model.geo.project=('Undefined');
       model.geo.version=136;
       model.geo.allmove=zeros(size(model.geo.symetric));
       model.geo.allmove_origin=0;
       model.geo.allmove_axis=0;
       model.geo.allmove_symetric=zeros(size(model.geo.symetric));
       model.geo.allmove_def=zeros(size(model.geo.symetric));
    end
    
    
    try model.geo.meshtype;
    
    catch
        resave=1;
        model.geo.meshtype=ones(size(model.geo.nx));
    end
    
    if resave
       disp(' ')
       disp(' + + + Warn + + +  Old file format detected ') 
       q5=input('Do you wish to rewrite in the new format? [1 0]: ','s');
       disp(' ')
       if q5
           cd('aircraft\')
                save(fname,'model');
           cd('..')
           
       else
           disp(' No file saved. ')
       end
    end
catch
    disp('Error in load_model.m')
end

end
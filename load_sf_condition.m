function condition = load_sf_condition(model)
% LOAD_SF_CONDITION - Loads a specific steady-flight condition from a database for a given model.
%
% This function provides an interactive interface to navigate a directory structure and load steady-flight (SF) 
% condition data from CSV files stored in a database. It allows the user to select a condition file based on the model name, 
% filter files by specific subtypes (e.g., 'C', 'D', 'L'), and choose a specific condition from the displayed data table. 
% The function returns a table row representing the selected condition with parameters such as Mach number, altitude, 
% velocities, and control inputs.
%
% INPUTS:
%   model      - A struct containing at least a 'name' field used to filter condition files specific to the model.
%
% OUTPUTS:
%   condition  - A table row containing the selected steady-flight condition with variables: 
%                'M', 'h [m]', 'u [m/s]', 'w [m/s]', 'theta [rad]', 'delta_e [rad]', 'T [N]', 'Tcom [%]', 'gamma [deg]'.
%
% FUNCTIONALITY:
% - Displays a menu of available condition files in the 'Database\conditions\' directory and allows the user to choose an option.
% - If the user selects the 'SF' directory (case 1), it navigates to subdirectories ('C', 'D', 'L') and filters CSV files 
%   based on the model name.
% - Uses the nested function `get_specific_sf_condition` to handle file selection within a subtype directory, read the CSV file, 
%   and return a specific condition row.
% - Loops until the user chooses to exit (case 0), returning to the parent directory upon exit.
%
% NOTES:
% - The function assumes a directory structure: 'Database\conditions\SF\<subtype>\' where subtype is 'C', 'D', or 'L'.
% - Only CSV files matching the model name are displayed for selection.
% - The nested function `get_specific_sf_condition` handles the detailed file reading and condition selection process.
% - Directory navigation uses relative paths, and the function returns to the original directory upon completion.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

% Nested Function: GET_SPECIFIC_SF_CONDITION - Retrieves a specific steady-flight condition from a subtype directory.
%
% This nested function navigates to a specified subtype directory (e.g., 'C\', 'D\', 'L\'), filters CSV files based on the 
% model name, displays the available files, reads the selected file into a table, and allows the user to pick a specific 
% condition row from the table.
%
% INPUTS:
%   subtype    - A string representing the subdirectory (e.g., 'C\', 'D\', 'L\') within the 'SF' directory.
%
% OUTPUTS:
%   condition  - A table row containing the selected steady-flight condition with named variables.
%
% FUNCTIONALITY:
% - Lists CSV files in the subtype directory that match the model name and end with '.csv'.
% - Reads the selected CSV file into a table, assigns meaningful variable names, and sets row names as indices.
% - Displays the table rounded to 3 decimal places and prompts the user to select a row by index.
% - Returns the selected row as the condition and navigates back to the root directory.
%
% NOTES:
% - The function assumes CSV files contain columns corresponding to the specified variable names.
% - Directory navigation is relative and resets to the original directory ('..\..\..\..\') after completion.

    while true
        clc;
        
        cd('Database\conditions\')
        files = dir(); % Get list of files and folders
        files = {files.name}; % Extract only names
        
        % Remove '.' and '..' from the list
        files = files(~ismember(files, {'.', '..'}));
        
        no_files = numel(files);
        disp(' ')
        for file_nr = 1:no_files
            disp(append('    [', string(file_nr),']. ', files{file_nr}))
        end
        disp(' ')
        disp('    [0]. Back / up menu')
        disp(' ')
        choice1 = input('    Enter a choice from above: ');
        disp(' ')
    
        switch choice1
            
    
            case(1)
                cd("SF\")
                clc;
                files = dir(); % Get list of files and folders
                files = {files.name}; % Extract only names
                
                % Remove '.' and '..' from the list
                files = files(~ismember(files, {'.', '..'}));
                
                no_files = numel(files);
                
                disp(' ')
                for file_nr = 1:no_files
                    disp(append('    [', string(file_nr),']. ', files{file_nr}))
                end
                disp(' ')
                choice2 = input('    Enter a choice from above: ');
                disp(' ')
    
                switch choice2
    
                    case(1)
                
                    condition = get_specific_sf_condition('C\');
    
                    case(2)
    
                    condition = get_specific_sf_condition('D\');
    
                    case(3)
    
                    condition = get_specific_sf_condition('L\');
    
                end
    
        case(0)
            cd('..\..\')
            break;
    
        end
    
    end
    
    
    function condition = get_specific_sf_condition(subtype)
        cd(subtype);
        files = dir(); % Get list of files and folders
        files = {files.name}; % Extract only names
        
        % Remove '.' and '..' from the list
        files = files(~ismember(files, {'.', '..'}));
        
        % Keep only .mat files that match model.name
        files = files(endsWith(files, '.csv') & contains(files, model.name));
    
        no_files = numel(files);
    
        % Display filtered file options
        for file_nr = 1:no_files
            disp(append('    [', string(file_nr),']. ', files{file_nr}))
        end
        disp(' ')
        file_nr = input('    Enter a choice from above: ');
    
        file = files{file_nr};
    
        table = readtable(file);
        
        table.Properties.RowNames = string(1:height(table));
        table.Properties.VariableNames = {'M', 'h [m]', 'u [m/s]', 'w [m/s]', 'theta [rad]', 'delta_e [rad]', 'T [N]', 'Tcom [%]', 'gamma [deg]'};
    
        disp(round(table, 3))
    
        idx = input('    Enter a choice from above: ');
    
        condition = table(idx,:);
    
        cd('..\..\..\..\')
    
    end
    


end
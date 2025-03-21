% ========================================================================
% Script Name: Model Management Script
% Author: Ludwig Horvath
% Date: 3/17/2025
% Description: 
%   This script continuously prompts for user input and performs different 
%   model-related operations based on the selection. It includes options to:
%     1. Assign a new geometry model
%     2. Load an existing model
%     3. Assign a model to the workspace
%     0. Exit the script
% 
%   Error handling is implemented to display appropriate messages if an 
%   operation fails.
% ========================================================================

clear; close all; clc

while true
    q=questions(1);

if isempty(q)
   q=-1; %will go to otherwise section
end


switch (q)
    case 1
        try
            [model] = assign_geometry();
        catch
            disp('Failed creating model...')
        end
    
    case 2
        try
            [model] = load_model();
        catch ME
           disp('Failed loading model...')
        end

    case 3
        try
            [model] = assign_model(model);
        catch ME
            disp('Failed assigning model..')
        end

    case 0
        clc;
        close all;
        break;

end

end




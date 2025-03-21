function [model] = assign_model(model)
% ========================================================================
% Script Name: assign_model.m
% Author: Ludwig Horvath
% Date: 3/17/2025
% Software: Redspot
% Description: 
%   This function manages the assignment and configuration of an aircraft 
%   model within Redspot. It provides options for defining geometry, 
%   estimating inertia, setting up actuation, defining dynamics, version 
%   control, and renaming the model.
%
%   The function interacts with the user through a question-based input 
%   system and includes error handling to ensure model data is saved in 
%   case of unexpected failures.
% ========================================================================



try


    while true
       answ=questions(2);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)
           
    
           case(1)
               
               % Geometry Setup
               [model] = inptgeometry(model);
    
           case(2)
    
               % Center of Gravity & Inertia Estimation
               [model] = inptinertia(model);


           case(3)
               [model] = inptactuation(model);
   
           case(4)
               % Dynamics Setup
	               
               [model] = inptbatch(model);
            
           case(5)
               % Version control
               version_control(model)
    
    
           case(6)
               % Rename model
               model.name = input('    Enter name of model: ', 's');
                
                
           case(0)
    
               break;
               clc;
    
       end
    
    end

catch
    % Handle errors here
    emergency_save(model);

end

end

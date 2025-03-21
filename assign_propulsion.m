function [model] = assign_propulsion(model)
% ========================================================================
% Script Name: assign_propulsion.m
% Author: Ludwig Horvath
% Date: 3/17/2025
% Software: Redspot
% Description: 
%   This function handles the assignment and configuration of propulsion 
%   parameters for an aircraft model within Redspot. It allows the user to:
%     - Select a propulsion model from available files
%     - Define maximum thrust
%     - Set a time constant for propulsion dynamics
%     - Visualize the step response of the propulsion system
%     - Inspect propulsion parameters in a dedicated UI
%
%   The function includes error handling to prevent data loss in case of 
%   unexpected failures.
% ========================================================================


try
    
    while true
       answ=questions(17);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)

           case 0
              break;

           case 1
               cd('aircraft\propulsion\')

                files = dir(); % Get list of files and folders
                files = {files.name}; % Extract only names
               
                cd('..\..\')

                % Remove '.' and '..' from the list
                files = files(~ismember(files, {'.', '..'}));
                
                % Keep only .mat files
                files = files(endsWith(files, '.m'));
                
                no_files = numel(files);

                clc;

                options = {};

                disp(' ')
                for file_nr = 1:no_files
                    disp(append('    [', string(file_nr),']. ', files{file_nr}))
                    options = {options,  files{file_nr}};
            
                end

                disp(' ')
                option = input(' 	Please enter choice from above: ');
                
                selected_file = files{option};
                
                model.parameters.propulsion.method = selected_file;

           case 2
               Tmax = input('    Enter maximum thrust [N]: ');

               if Tmax > 0
                  model.parameters.propulsion.Tmax = Tmax;
               else
                   disp('    Invalid value Tmax must be positive')
               end
           case 3
                tau = input('    Enter time constant of 1st order system: ');
                
                if tau >= 0
                    model.parameters.propulsion.tau = tau;
                else
                    disp('Invalid time constant')
                end

           case 4
                figure_name = append('Name', append('Step response of propulsion control input'));

                 % Refine fixed figure size (in pixels)
                fig_width  = 600;  
                fig_height = 500;  
                screen_size = get(0, 'ScreenSize'); % Get screen size
                
                % Center the figure on screen
                fig_left = (screen_size(3) - fig_width) / 2;
                fig_bottom = (screen_size(4) - fig_height) / 2;
                
                
                % Create fixed-size figure
                fig = figure('Name', figure_name, ...
                      'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
                      'Resize', 'off', ...  % Prevent resizing
                      'Toolbar', 'none', 'Menubar', 'figure', ...
                       'Color', 'white');
                
                
                axes('Parent', fig, 'Position', [0.1, 0.2, 0.7, 0.5]); 
                
                

                tau = model.parameters.propulsion.tau;

                t = linspace(0, 100*tau, 100);
                y = exp(-tau./t);
               
                tau = model.parameters.propulsion.tau;

                plot(t, y, 'DisplayName', append('Step response, propulsion', 'tau =', string(tau)), 'Color', 'k')
                hold on
                
                yline(1, 'LineStyle', '--', 'Color', 'k')

                model.parameters.propulsion.info = stepinfo(y,t); 
                
                ylim([0 1.5])

                xlabel('t [s]')
                ylabel('Tcom [%]')
                grid minor

                annotation("textbox", "String", append('  Time constant:   ',num2str(round(tau, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.80, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('          Rise time:   ',num2str(round(model.parameters.propulsion.info.RiseTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.75, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('     Settling time:   ',num2str(round(model.parameters.propulsion.info.SettlingTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.70, 0.1, 0.1], ...
                           "EdgeColor", 'white');
                            

           case 5
                 figure_name = append('Inspection, Parameters: propulsion');
                
                 % Refine fixed figure size (in pixels)
                fig_width  = 600;  
                fig_height = 150;  
                screen_size = get(0, 'ScreenSize'); % Get screen size
                
                % Center the figure on screen
                fig_left = (screen_size(3) - fig_width) / 2;
                fig_bottom = (screen_size(4) - fig_height) / 2;
                
                
                % Create fixed-size figure
                figure('Name', figure_name, ...
                      'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
                      'Resize', 'off', ...  % Prevent resizing
                      'Toolbar', 'none', 'Menubar', 'figure', ...
                       'Color', 'white');

                

                annotation("textbox", "String", append('  Time constant:   ',num2str(round(model.parameters.propulsion.tau, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.80, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('          Rise time:   ',num2str(round(model.parameters.propulsion.info.RiseTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.50, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('     Settling time:   ',num2str(round(model.parameters.propulsion.info.SettlingTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.20, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('   Tmax: ', num2str(model.parameters.propulsion.Tmax)), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.45, 0.70, 0.1, 0.1], ...
                           "EdgeColor", 'white');
                annotation("textbox", "String", append('method: ', num2str(model.parameters.propulsion.method)), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.45, 0.40, 0.1, 0.1], ...
                           "EdgeColor", 'white');

       end

    end


catch

end


end
function [model] = assign_controlsurface(model, type)


try
    while true
       answ=questions(15);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)

           case 0
                break;

           case 1

                flap_matrix = model.geo.flapped;

                Nsections = width(flap_matrix);
                Nsurfaces = height(flap_matrix);
            
                % Generate row names ("Surface 1", "Surface 2", etc.)
                row_names = arrayfun(@(x) sprintf('Surface %d', x), 1:Nsurfaces, 'UniformOutput', false);
            
                % Generate column names ("Section 1", "Section 2", etc.)
                col_names = arrayfun(@(x) sprintf('Sectione %d', x), 1:Nsections, 'UniformOutput', false);
            
                % Convert matrix to table with proper row and column names
                flap_table = array2table(flap_matrix, 'RowNames', row_names, 'VariableNames', col_names);
            
                disp(append('Matrix "flapped" of ', model.name, ': '))
                disp(' ')
                disp([flap_table, array2table(model.parameters.(type).lmatrix)])
                
                cell = input('Add/Remove corresponding cell [row, column]: ');

                if model.parameters.(type).lmatrix(cell(1), cell(2)) == 0
                   model.parameters.(type).lmatrix(cell(1), cell(2)) = 1;
                else  
                   model.parameters.(type).lmatrix(cell(1), cell(2)) = 0;
                end
                
                disp(append('Matrix "flapped" of ', model.name, ': '))
                disp(' ')
                disp([flap_table, array2table(model.parameters.(type).lmatrix)])

                input('Press Enter to continue.. ')

           case 2
               max = input(append('Enter maximum deflection, ', string(type) ,' [deg]: '));
               
               if all([90 >= max, max >= 0])
                    model.parameters.(type).max = max;
               end
            
               min = input(append('Enter minimum deflection, ', string(type) ,' [deg]: '));
        
               if all([0 >= min, min >= -90])
                    model.parameters.(type).min = min;
               end

           case 3   %time constant
                tau = input('    Enter time constant of 1st order system: ');
                
                if tau >= 0
                    model.parameters.(type).tau = tau;
                else
                    disp('Invalid time constant')
                end

           case 4   % step
                                
                figure_name = append('Step response of ', type);
                
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

                
                tau = model.parameters.(type).tau;

                t = linspace(0, 100*tau, 100);
                y = exp(-tau./t);

                model.parameters.(type).info = stepinfo(y, t); % Compute step response characteristics
                model.parameters.(type).info 

                annotation("textbox", "String", append('  Time constant:   ',num2str(round(tau, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.80, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('          Rise time:   ',num2str(round(model.parameters.(type).info.RiseTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.75, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('     Settling time:   ',num2str(round(model.parameters.(type).info.SettlingTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.70, 0.1, 0.1], ...
                           "EdgeColor", 'white');
                            
                
                plot(t, y, 'DisplayName', append('Step response, ', "delta_", string(type(1)), ', tau =', string(tau)), 'Color', 'k')
                hold on
                yline(1, 'LineStyle', '--', 'Color', 'k')
                
                ylim([0 1.5])

                xlabel('t [s]')
                ylabel(append(type, ' [deg]'))
                grid minor
               

           case 5   % overview

               figure_name = append('Inspection , Parameters: ', type);
                
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
                
                annotation("textbox", "String", append('  Time constant:   ',num2str(round(model.parameters.(type).tau, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.80, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('          Rise time:   ',num2str(round(model.parameters.(type).info.RiseTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.50, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('     Settling time:   ',num2str(round(model.parameters.(type).info.SettlingTime, 4)), ' s'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.1, 0.20, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                annotation("textbox", "String", append('max: ', num2str(model.parameters.(type).max), ' deg'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.45, 0.70, 0.1, 0.1], ...
                           "EdgeColor", 'white');
                annotation("textbox", "String", append('min: ', num2str(model.parameters.(type).min), ' deg'), ...
                           "FontSize", 10, ...
                           "FitBoxToText", 'on', ...
                           "Position", [0.45, 0.40, 0.1, 0.1], ...
                           "EdgeColor", 'white');

                lmatrix = model.parameters.(type).lmatrix;

                matrix_str = "";
                for i = 1:size(lmatrix, 1)
                    matrix_str = matrix_str + join(string(lmatrix(i, :)), '  ') + newline;
                end
                

                annotation("textbox", ...
                           "String", matrix_str, ...
                           "FontSize", 10, ...
                           "FontName", "Courier", ... % Fixed-width font for alignment
                           "FitBoxToText", 'on', ...
                           "Position", [0.7, 0.50, 0.4, 0.4], ...
                           "EdgeColor", 'white', ...
                           "BackgroundColor", 'white');

       end

    end

catch 


end

end

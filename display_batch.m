function display_batch(batch)
% DISPLAY_STATE_BATCH - Visualizes the batch grid points for a given model.
%
% This function creates a fixed-size figure displaying the state batch 
% of various aerodynamic parameters. Each state is represented as a vertical 
% grid line in a subplot, with minimum and maximum values highlighted in 
% red and blue, respectively.
%
% INPUT:
%   model - A struct containing the field `batch`, 
%           which stores state variable limits and grid values.
%
% FUNCTIONALITY:
% - Retrieves and formats state variable names.
% - Sets up a fixed-size figure for visualization.
% - Iterates through available state variables and plots grid lines.
% - Highlights minimum values in red and maximum values in blue.
% - Uses gray for negative values in symmetric state variables.
% - Adds annotations for variable names and units.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

batch_info = batch.info;

state_labels_clean = {'alpha', ...
                      'alpha_dot', ...
                      'beta',  ...
                      'beta_dot', ...
                      'P', ...
                      'Q', ...
                      'R', ...
                      'M', ...
                      'h', ...
                      'delta_a', ...
                      'delta_e', ...
                      'delta_r'};


formatted_names = struct();
formatted_names.alpha       = 'alpha';
formatted_names.alpha_dot   = 'alpha\_dot';
formatted_names.beta        = 'beta';
formatted_names.beta_dot    = 'beta\_dot';
formatted_names.P           = 'P';
formatted_names.Q           = 'Q';
formatted_names.R           = 'R';
formatted_names.M           = 'M';
formatted_names.h           = 'h';
formatted_names.delta_a     = 'delta\_a';
formatted_names.delta_e     = 'delta\_e';
formatted_names.delta_r     = 'delta\_r';


try
    % Define fixed figure size
    fig_width = 800;  
    fig_height = 600;  
    screen_size = get(0, 'ScreenSize'); % Get screen size

    % Center the figure on screen
    fig_left = (screen_size(3) - fig_width) / 2;
    fig_bottom = (screen_size(4) - fig_height) / 2;

    % Create fixed-size figure
    figure('Name', append("Overview ", batch_info.id), ...
          'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
          'Resize', 'off', ...  % Prevent resizing
          'Toolbar', 'none', 'Menubar', 'figure', ...
          'Color', 'white');

    % Get field names from batch_info
    batch_info_fields = state_labels_clean;
    num_states = numel(batch_info_fields);

    % Left side for subplots (only vertical lines)
    for i = 1:num_states
        field_name = batch_info_fields{i};
        if isfield(batch_info, field_name)
            
            state_data = batch_info.(field_name); 
            min_val = state_data.min;
            max_val = state_data.max;

            subplot(num_states, 1, i); % Stack plots vertically
            hold on;
            
            % Plot vertical bars for each grid point
            for j = 1:numel(state_data.grid)

                x = state_data.grid(j);

                if state_data.N == 1
                    xline(x, 'k', 'LineWidth', 1); % Black for other grid points

                else

                    if x == min_val
                        xline(x, 'r', 'LineWidth', 1.5); % Red for min
                    elseif x == max_val
                        xline(x, 'b', 'LineWidth', 1.5); % Blue for max
                    else
                        xline(x, 'k', 'LineWidth', 1); % Black for other grid points
                    end
                    
                end

            end

            if state_data.N == 1
                set(gca, 'XTick', min_val, ...
                         'XTickLabel', string(min_val));

            elseif state_data.N ~= 1
            
                % Set only min and max values as XTicks
                set(gca, 'XTick', [min_val, max_val], ...
                         'XTickLabel', string([min_val, max_val]));
        
             end

            % Hide y-axis ticks
            set(gca, 'YTick', []);
            
            % Remove box
            box off;
            
            corr_factor = 0.9;
    
            % Add annotation to the left of each subplot
            annotation('textbox', [0, 0.95 - corr_factor*(i / (num_states +1)), 0.1, 0.05], ...
                       'String', formatted_names.(field_name), ... % Escape underscore for display
                       'EdgeColor', 'none', ...
                       'FontSize', 10, ...
                       'HorizontalAlignment', 'right');
    
            annotation('textbox', [0.87, 0.95 - corr_factor*(i / (num_states +1)), 0.1, 0.05], ...
                       'String', strrep(state_data.unit, '_', '\_'), ... % Escape underscore for display
                       'EdgeColor', 'none', ...
                       'FontSize', 10, ...
                       'HorizontalAlignment', 'right');
    
            hold off;
        end

    end

catch
    disp('Error in display_batch_info.m');
end
end


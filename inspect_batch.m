function inspect_batch(batch_info)
    % Extract field names (excluding 'id' and 'progress')
    state_fields = setdiff(fieldnames(batch_info), {'id', 'progress'}, 'stable');

    % Define row labels
    row_labels = {'N', 'max', 'cv', 'min', 'unit'};
    
    % Initialize empty array for the table data
    num_rows = length(row_labels);
    num_cols = length(state_fields);
    table_data = cell(num_rows, num_cols);
    
    % Populate the table data
    for col = 1:num_cols
        field_name = state_fields{col};
        for row = 1:num_rows
            subfield_name = row_labels{row};

            % Extract the value
            value = batch_info.(field_name).(subfield_name);


            if iscell(value) && isscalar(value)
                value = value{1};  % Extract number from cell
            end


            if isnumeric(value) && isscalar(value)
                table_data{row, col} = num2str(value);  % Force proper formatting
            elseif ischar(value) || isstring(value)
                table_data{row, col} = char(value);  % Convert string to char array
            else
                table_data{row, col} = '';  % Fallback for unknown types
            end
        end
    end
    

    T = cell2table(table_data, 'VariableNames', state_fields, 'RowNames', row_labels);
    

    disp(T);
end

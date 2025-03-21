function batch = load_batch(model)
    clc;
    cd('Database\Batches\')
    
    files = dir(); % Get list of files and folders
    files = {files.name}; % Extract only names
    
    % Remove '.' and '..' from the list
    files = files(~ismember(files, {'.', '..'}));
    
    % Keep only .mat files that match model.name
    files = files(endsWith(files, '.mat') & contains(files, model.name));
    
    no_files = numel(files);

    if no_files == 0
        disp('No matching batch files found.');
        cd('..')
        cd('..')
        batch = struct();
        batch.info = [];
        batch.data = [];
        return;
    end

    while true
        try
            clc;
            disp(' ')
            disp('    [1]. Inspect a batch')
            disp('    [2]. Load a batch')
            disp(' ')
            disp('    [0]. Back / up menu')
            disp(' ')
            choice = input(' 	Please enter choice from above: ');
        
            if choice ~= 0
                clc;
                disp(' ')

                % Display filtered file options
                for file_nr = 1:no_files
                    disp(append('    [', string(file_nr),']. ', files{file_nr}))
                end

                disp(' ')
                option = input(' 	Please enter choice from above: ');

                selected_file = files{option};

                if choice == 1
                    load(selected_file, 'info');
                    cd('..')
                    cd('..')
                    inspect_batch(info)
                    input('Press Enter to continue ...')
                    clear('batch')
                    cd('Database\Batches\')
        
                elseif choice == 2
                    info_struct = load(selected_file);
                    batch.info = info_struct.info;
                    try
                        batch.data = readtable(strrep(selected_file, '.mat', '.csv'));
                    catch 
                        disp(' ')
                        disp('    Batch.data is empty (data non-existent)')
                        disp(' ')
                    end
                    cd('..')
                    cd('..')
                    break;
                end
                
            else
                cd('..')
                cd('..')
                batch = struct();
                batch.info = [];
                batch.data = [];
                break;
            end

        catch 
            disp('Error in "load_batch.m"')
        end
    end
end

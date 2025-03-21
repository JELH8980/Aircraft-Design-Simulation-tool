function solve_batch(batch, model)
% SOLVE_BATCH - Solves a batch of aerodynamic conditions and stores results in a CSV file with real-time progress monitoring.
%
% This function processes a batch of aerodynamic conditions defined by a grid of parameters (e.g., angle of attack, Mach number, 
% control surface deflections) for a given model. It performs vortex lattice method (VLM) calculations to compute aerodynamic 
% coefficients and their derivatives, appends results to a CSV file, and updates a progress monitoring figure with estimated 
% completion time. The function handles both static and dynamic variables, supports error handling, and saves progress incrementally.
%
% INPUTS:
%   batch      - A struct containing batch information:
%                - batch.info.id: unique identifier for the batch
%                - batch.info.progress: starting index of the batch (for resuming)
%                - batch.info.<parameter>.grid: grid values for parameters (e.g., alpha, beta, M, etc.)
%                - batch.info.<parameter>.N: number of grid points for each parameter
%   model      - A struct containing model data:
%                - model.geo: geometry data for VLM calculations
%                - model.parameters: physical parameters (e.g., delta_a.lmatrix for control surfaces)
%
% OUTPUTS:
%   None       - Results are written to a CSV file at 'database\batches\<batch.info.id>.csv', and progress is saved to a .mat file.
%
% FUNCTIONALITY:
% - Creates or appends to a CSV file with headers for aerodynamic coefficients (CD, CL, etc.) and their derivatives.
% - Defines a grid of conditions based on batch parameters, separating constants (N=1) and variables (N>1).
% - Iterates over the grid starting from batch.info.progress, computing aerodynamic forces and moments using VLM.
% - Updates a real-time progress figure with estimated duration, worst-case/best-case completion times, and current iteration status.
% - Saves progress to a .mat file after each iteration to allow resuming if interrupted.
% - Handles errors gracefully with a try-catch block, displaying an error message if execution fails.
%
% NOTES:
% - External functions (fLattice_setup, solver, coeff_create, funsteady, ISAtmosphere, set_flaps) are assumed to be defined elsewhere.
% - The progress figure includes plots for predicted duration, sweeper lines, and annotations for estimated time left and completion date.
% - Aerodynamic coefficients are scaled (e.g., CL from CZ) and adjusted for unsteady effects (e.g., alpha_dot, beta_dot).
% - Units are converted from degrees to radians for state variables (e.g., alpha, P, Q, R) during computation.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

try

    batch_data_path =append('database\batches\', batch.info.id, '.csv');

    batch_info_path   = append('database\batches\', batch.info.id, '.mat');
    


    if ~exist(batch_data_path, 'file')

        fileID = fopen(batch_data_path, 'w');  
        fprintf(fileID, ['alpha,alpha_dot,beta,beta_dot,P,Q,R,M,h,delta_a,delta_e,delta_r,CD,CC,CL,Cl,Cm,Cn,' ...
                         'CD_a,CD_a_dot,CD_b,CD_P,CD_Q,CD_R,CD_delta_a,CD_delta_e,CD_delta_r,' ...
                         'CC_a,CC_b,CC_b_dot,CC_P,CC_Q,CC_R,CC_delta_a,CC_delta_e,CC_delta_r,' ...
                         'CL_a,CL_a_dot,CL_b,CL_P,CL_Q,CL_R,CL_delta_a,CL_delta_e,CL_delta_r,' ...
                         'Cl_a,Cl_b,Cl_P,Cl_Q,Cl_R, Cl_delta_a,Cl_delta_e,Cl_delta_r,' ...
                         'Cm_a,Cm_a_dot,Cm_b,Cm_P,Cm_Q,Cm_R,Cm_delta_a,Cm_delta_e,Cm_delta_r,' ...
                         'Cn_a,Cn_b,Cn_b_dot,Cn_P,Cn_Q,Cn_R,Cn_delta_a,Cn_delta_e,Cn_delta_r\n']);  

        fclose(fileID);  % Close file

    end


    % Get the lengths of each grid variable
    alpha_vals       = batch.info.alpha.grid;
    alpha_dot_vals   = batch.info.alpha_dot.grid;
    beta_vals        = batch.info.beta.grid;
    beta_dot_vals    = batch.info.beta_dot.grid;
    P_vals           = batch.info.P.grid;
    Q_vals           = batch.info.Q.grid;
    R_vals           = batch.info.R.grid;
    M_vals           = batch.info.M.grid;
    h_vals           = batch.info.h.grid;
    delta_e_vals     = batch.info.delta_e.grid;
    delta_a_vals     = batch.info.delta_a.grid;
    delta_r_vals     = batch.info.delta_r.grid;

    Nstart = batch.info.progress;

    % Calculate total number of grid points
    N = numel(M_vals)*numel(h_vals)*numel(alpha_vals)*numel(alpha_dot_vals)*numel(beta_vals)*numel(beta_dot_vals)*numel(P_vals)*numel(Q_vals)*numel(R_vals) * numel(delta_e_vals) * numel(delta_a_vals) * numel(delta_r_vals) - Nstart + 1;
    
    % Initialize Figure
    figure_name = append('Monitor');
    fig_width = 600;
    fig_height = 500;
    screen_size = get(0, 'ScreenSize');
    fig_left = (screen_size(3) - fig_width) / 2;
    fig_bottom = (screen_size(4) - fig_height) / 2;
    
    % Create fixed-size figure
    fig = figure('Name', figure_name, ...
                 'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
                 'Resize', 'off', 'Toolbar', 'none', 'Menubar', 'figure', ...
                 'Color', 'white');
    
    axes('Parent', fig, 'Position', [0.1, 0.55, 0.8, 0.4]); 
    
    hold on;
   
    duration_plot = plot(nan, nan, 'k', 'LineWidth', 0.5);
   
    xlabel('N [-]')
    ylabel('T [h]')
    
    sweeper_vaxis = xline(0, 'Color', 'b', 'LineWidth', 0.5);

    sweeper_haxis = line([0 N], [0, 0], 'Color', 'b', 'LineWidth', 0.5);
    
    xline(N, 'Color', 'k', 'LineWidth', 0.5);
    xline(N*1.01, 'Color', 'k', 'LineWidth', 0.5);

    end_dot      = plot(0, 0, 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    

    wc_endtime = 0;
    bc_endtime = inf;

    exp_endtime = 0;

    % Create horizontal lines
    wc_axis  = line([0 N], [wc_endtime, wc_endtime]  , 'Color', 'r', 'LineWidth', 0.5);
    bc_axis  = line([0 N], [bc_endtime, bc_endtime]  , 'Color', 'g', 'LineWidth', 0.5);
    exp_axis = line([0 N], [exp_endtime, exp_endtime], 'Color', 'k', 'LineWidth', 0.5);

    exp_timeiteration = 0;

    exp_timeleft = 0;

    exp_date      = datetime('now');
    
    annotation("textbox", "String", "Est time left: ", ...
               "FontSize", 10, ...
               "FitBoxToText", 'on', ...
               "Position", [0.1, 0.1, 0.1, 0.1], ...
               "EdgeColor", 'white');

    exp_time_annotation = annotation("textbox", "String", string(exp_timeleft), ...
                                      "FontSize", 10, ...
                                      "FitBoxToText", 'on', ...
                                      "Position", [0.4, 0.1, 0.1, 0.1], ...
                                      "EdgeColor", 'white');

    annotation("textbox", "String", "Est time of completion: ", ...
               "FontSize", 10, ...
               "FitBoxToText", 'on', ...
               "Position", [0.1, 0.2, 0.1, 0.1], ...
               "EdgeColor", 'white');

    exp_date_annotation = annotation("textbox", "String", string(exp_date), ...
                                      "FontSize", 10, ...
                                      "FitBoxToText", 'on', ...
                                      "Position", [0.4, 0.2, 0.1, 0.1], ...
                                      "EdgeColor", 'white');

    annotation("textbox", "String", "Time/Iteration: ", ...
               "FontSize", 10, ...
               "FitBoxToText", 'on', ...
               "Position", [0.1, 0.3, 0.1, 0.1], ...
               "EdgeColor", 'white');

    timeiteration_annotation = annotation("textbox", "String", string(exp_timeiteration), ...
                                      "FontSize", 10, ...
                                      "FitBoxToText", 'on', ...
                                      "Position", [0.4, 0.3, 0.1, 0.1], ...
                                      "EdgeColor", 'white');

        
    parameters = {'alpha', 'alpha_dot', 'beta', 'beta_dot', 'P', 'Q', 'R', 'M', 'h', 'delta_a', 'delta_e', 'delta_r'};

    variables = {};
    constants = struct();
    grid_input_string = "";
    grid_output_string = "";

    for p = 1:numel(parameters)
        np = batch.info.(parameters{p}).N;
        if np == 1
            constants.(parameters{p}) = batch.info.(parameters{p}).grid;
        else
            variables{end+1} = parameters{p};

            if ~strcmp(grid_input_string, "")
                input_suffix = strcat(',', parameters{p}, '_vals');
                output_suffix = strcat(',', parameters(p), '_Grid');
            else
                input_suffix = strcat(parameters{p}, '_vals');
                output_suffix = strcat(parameters(p), '_Grid');
            end

            grid_input_string = grid_input_string + input_suffix;
            grid_output_string = grid_output_string + output_suffix;

        end
    end

    function_call_1 = strcat('[', grid_output_string, ']=ndgrid(',grid_input_string ,');');
    
    eval(function_call_1);

    for v = 1:numel(variables)
        grid_flattener_string_rhs = strcat(variables{v}, '_Grid', '(:);');
        grid_flattener_string_lhs = strcat(variables{v}, '_Flat');
        function_call_2 = strcat(grid_flattener_string_lhs, '=', grid_flattener_string_rhs);
        eval(function_call_2);
    end

    for idx = Nstart:N

        tic;

        if isfield(constants, 'alpha')
            alpha = constants.alpha;
        else
            alpha = alpha_Flat(idx);
        end

        if isfield(constants, 'alpha_dot')
            alpha_dot = constants.alpha_dot;
        else
            alpha_dot = alpha_dot_Flat(idx);
        end

        if isfield(constants, 'beta')
            beta = constants.beta;
        else
            beta = beta_Flat(idx);
        end

        if isfield(constants, 'beta_dot')
            beta_dot = constants.beta_dot;
        else
            beta_dot = beta_dot_Flat(idx);
        end

        if isfield(constants, 'P')
            P = constants.P;
        else
            P = P_Flat(idx);
        end

        if isfield(constants, 'Q')
            Q = constants.Q;
        else
            Q = Q_Flat(idx);
        end

        if isfield(constants, 'R')
            R = constants.R;
        else
            R = R_Flat(idx);
        end
    
        if isfield(constants, 'M')
            M = constants.M;
        else
            M = M_Flat(idx);
        end

        if isfield(constants, 'h')
            h = constants.h;
        else
            h = h_Flat(idx);
        end

        if isfield(constants, 'delta_a')
            delta_a = constants.delta_a;
        else
            delta_a = delta_a_Flat(idx);
        end

        if isfield(constants, 'delta_e')
            delta_e = constants.delta_e;
        else
            delta_e = delta_e_Flat(idx);
        end

        if isfield(constants, 'delta_r')
            delta_r = constants.delta_r;
        else
            delta_r = delta_r_Flat(idx);
        end

                

        % Access the corresponding grid values
        state.alpha     = alpha      * pi / 180;
        state.alphadot  = alpha_dot  * pi / 180;
        state.bethadot  = beta_dot   * pi / 180;
        state.betha     = beta       * pi / 180;
        state.P         = P          * pi / 180;
        state.Q         = Q          * pi / 180;
        state.R         = R          * pi / 180;
        
        speed = M;
        state.ALT = h;
        [state.rho, a] = ISAtmosphere(state.ALT);
        state.AS = speed * a;
        
        
        if  M > 0.3
            state.pgcorr = 1;                           % Correction for viscous effects (Prandtl-Gaertner correction factor)
        else
            state.pgcorr = 0;    
        end

        
        % Setting Control Surface Deflections
        model = set_flaps(model, "delta_a", delta_a);
        model = set_flaps(model, "delta_e", delta_e);
        model = set_flaps(model, "delta_r", delta_r);

        [delta_a_indx, ~] = find(model.parameters.delta_a.lmatrix == 1, 1, 'first');
        [delta_e_indx, ~] = find(model.parameters.delta_e.lmatrix == 1, 1, 'first');
        [delta_r_indx, ~] = find(model.parameters.delta_r.lmatrix == 1, 1, 'first');
        
        delta_indices = [delta_a_indx, delta_e_indx, delta_r_indx];
    
        [~, indices] = sort(delta_indices);
    
        delta_a_indx = indices(delta_indices == delta_a_indx);
        delta_e_indx = indices(delta_indices == delta_e_indx);
        delta_r_indx = indices(delta_indices == delta_r_indx);

        % Performing Vortex Lattice Calculations
        [lattice, ref] = fLattice_setup(model.geo, state, 0);
        resultVLM.dwcond=0;
        [results] = solver(resultVLM, state, model.geo, lattice, ref);
        [results] = coeff_create(results, lattice, state, ref, model.geo);
        
        A=funsteady(lattice,model.geo,ref,M);
        
        results.CZ_a_dot=-A(1);
        results.Cm_a_dot=-A(2);
        
        %Flipping the lattice to get the sideslipe solved.
        lattice.COLLOC=[lattice.COLLOC(:,1) lattice.COLLOC(:,3) lattice.COLLOC(:,2)];
        lattice.N=[lattice.N(:,1) lattice.N(:,3) lattice.N(:,2)];
    
        lattice.VORTEX2(:,:,1)=lattice.VORTEX(:,:,1);
        lattice.VORTEX2(:,:,3)=lattice.VORTEX(:,:,2);
        lattice.VORTEX2(:,:,2)=lattice.VORTEX(:,:,3);
        lattice.VORTEX=lattice.VORTEX2;
    
        lattice.XYZ2(:,:,1)=lattice.XYZ(:,:,1);
        lattice.XYZ2(:,:,3)=lattice.XYZ(:,:,2);
        lattice.XYZ2(:,:,2)=lattice.XYZ(:,:,3);
        lattice.XYZ=lattice.XYZ2;
    
        A=funsteady(lattice,model.geo,ref,M);
        
        results.CY_b_dot=-A(1);
        results.Cn_b_dot=-A(2);
        
        factor1=results.CL/results.CZ;  %L scaling factor  
        if isnan(factor1)
            factor1=1;
        end
        
        factor2=results.CD/results.CZ;  %D scaling factor  
        if isnan(factor2)
            factor2=0;
        end
     
        results.CL_a_dot=results.CZ_a_dot*factor1;
        results.CD_a_dot=results.CZ_a_dot*factor2;
        results.CC_b_dot=results.CY_b_dot;
           
        % Storing derivatives into result_matrix
        result_stat = [results.CD  , results.CC  , results.CL  , results.Cl  , results.Cm  , results.Cn];

        result_dyn =  [results.CD_a, results.CD_a_dot, results.CD_b,                   results.CD_P, results.CD_Q, results.CD_R, results.CD_d(delta_a_indx), results.CD_d(delta_e_indx), results.CD_d(delta_r_indx),...
                       results.CC_a,                   results.CC_b, results.CC_b_dot, results.CC_P, results.CC_Q, results.CC_R, results.CC_d(delta_a_indx), results.CC_d(delta_e_indx), results.CC_d(delta_r_indx),...
                       results.CL_a, results.CL_a_dot, results.CL_b,                   results.CL_P, results.CL_Q, results.CL_R, results.CL_d(delta_a_indx), results.CL_d(delta_e_indx), results.CL_d(delta_r_indx),...
                       results.Cl_a,                   results.Cl_b,                   results.Cl_P, results.Cl_Q, results.Cl_R, results.Cl_d(delta_a_indx), results.Cl_d(delta_e_indx), results.Cl_d(delta_r_indx),...
                       results.Cm_a, results.Cm_a_dot, results.Cm_b,                   results.Cm_P, results.Cm_Q, results.Cm_R, results.Cm_d(delta_a_indx), results.Cm_d(delta_e_indx), results.Cm_d(delta_r_indx)...
                       results.Cn_a,                   results.Cn_b, results.Cn_b_dot, results.Cn_P, results.Cn_Q, results.Cn_R, results.Cn_d(delta_a_indx), results.Cn_d(delta_e_indx), results.Cn_d(delta_r_indx)];
                            

        % Specify conditions first
        writematrix([alpha, alpha_dot, beta, beta_dot, P, Q, R,  M, h, delta_a, delta_e, delta_r, result_stat, result_dyn], batch_data_path, 'WriteMode', 'append');

        batch.info.progress     = batch.info.progress + 1;

        info = batch.info;

        save(batch_info_path, "info", "-mat")

        % Update progress in real-time

        time = toc / 3600;
       
        duration_model = fitlm([0, 1], [0, time]);
        horizon = (1:N)';
        pred_times = predict(duration_model, horizon);
        set(duration_plot, 'XData',  horizon, 'YData', pred_times);
        
        pred_endtime = pred_times(end);

        if pred_endtime >= wc_endtime
            wc_endtime = pred_endtime;

        end
        
        if pred_endtime <= bc_endtime
            bc_endtime = pred_endtime;
        end

        % Update the YData property to move the horizontal lines
        set(wc_axis,  'YData', [wc_endtime, wc_endtime]);
        set(bc_axis,  'YData', [bc_endtime, bc_endtime]);

        if idx == 1
            exp_endtime = pred_endtime;
        end

        exp_endtime = (exp_endtime*idx + pred_endtime)/(idx+1);

        set(exp_axis, 'YData', [exp_endtime, exp_endtime]);
        
        set(sweeper_vaxis, 'Value', idx)
        set(sweeper_haxis, 'YData', [pred_endtime, pred_endtime]);

        set(end_dot, 'XData', N, 'YData', pred_endtime);
         
        drawnow;

        % Annotations    
     
        pred_totaltime_left = pred_endtime*(1-idx/N);

        pred_enddate_total = datetime('now') + hours(pred_totaltime_left);
        
        exp_time_annotation.set('String', append(string(pred_totaltime_left), ' h'))
        exp_date_annotation.set('String', append(string(pred_enddate_total)))
        timeiteration_annotation.set('String', append(string(time*3600), ' s'))
        

end



catch
    disp('Error in "solve_batch.m"')
end
end


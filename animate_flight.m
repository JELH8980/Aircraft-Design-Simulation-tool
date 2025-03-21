function animate_flight(flightdata, reference, speed)
    

    % Figure 1: Nonlinear State-Time trajectories
    figure_name = append('Animation');
    
    % Refine fixed figure size (in pixels)
    fig_width = 1300;  
    fig_height = 800;  
    screen_size = get(0, 'ScreenSize'); % Get screen size
    
    % Center the figure on screen
    fig_left = (screen_size(3) - fig_width) / 2;
    fig_bottom = (screen_size(4) - fig_height) / 2;
    
    
    % Create fixed-size figure
    fig = figure('Name', figure_name, ...
                 'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
                 'Resize', 'on', ...  % Prevent resizing
                 'Color', 'white');
    
    hold on

    dt = mean(diff(flightdata.("t [s]")));
    
    xE =   flightdata.("xE [m]")(1);
    
    yE =   flightdata.("yE [m]")(1);
    
    zE =   flightdata.("zE [m]")(1);
    
    h  = - zE + reference.h;
    
    position = plot3(xE, yE, h, 'Marker','o', 'MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'LineWidth',2, 'MarkerSize', 5);
    
    phi   = flightdata.("phi [rad]")(1);
    theta = flightdata.("theta [rad]")(1);
    psi   = flightdata.("psi [rad]")(1);
    
    
    R_BV = [cos(theta)*cos(psi), (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)), (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
            cos(theta)*sin(psi), (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)), (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
                    -sin(theta),                             sin(phi)*cos(theta),                            cos(phi)*cos(theta)];
    
    R_VE = [1,  0, 0;
            0, -1, 0;
            0,  0,-1];
    
    
    l = 1000;
    
    x_axis_B = l*[1, 0, 0]';
    y_axis_B = l*[0, 1, 0]';
    z_axis_B = l*[0, 0, 1]';
    
    x_axis_E = R_VE*R_BV*x_axis_B;
    y_axis_E = R_VE*R_BV*y_axis_B;
    z_axis_E = R_VE*R_BV*z_axis_B;
    
    x_axis_vector = [xE, xE + x_axis_E(1), yE, yE + x_axis_E(2), h, h + x_axis_E(3)];
    y_axis_vector = [xE, xE + y_axis_E(1), yE, yE + y_axis_E(2), h, h + y_axis_E(3)];
    z_axis_vector = [xE, xE + z_axis_E(1), yE, yE + z_axis_E(2), h, h + z_axis_E(3)];
    
    
    x_axis = plot3([x_axis_vector(1), x_axis_vector(2)], [x_axis_vector(3), x_axis_vector(4)], [x_axis_vector(5), x_axis_vector(6)], 'color', 'r');
    y_axis = plot3([y_axis_vector(1), y_axis_vector(2)], [y_axis_vector(3), y_axis_vector(4)], [y_axis_vector(5), y_axis_vector(6)], 'color', 'b');
    z_axis = plot3([z_axis_vector(1), z_axis_vector(2)], [z_axis_vector(3), z_axis_vector(4)], [z_axis_vector(5), z_axis_vector(6)], 'color', 'g');
    
    trajectory = animatedline('Color', 'k');
    
    ground_track = animatedline('Color', '[0.5, 0.5, 0.5]', 'LineStyle', '--');
    
    view(3)
    
    xlabel('xE [m]')
    ylabel('yE [m]')
    zlabel('h [m]')
    
    t = 0;
    
    K_clock = 1;
    
    wait_time = dt;
    
    trel = 0;
    
    tstep = 1;
    
    trel_map = 0;
    
    tstep_map = 1;
    
    kstep = 10;
    
    
    uif = uifigure("Name", 'Geographical view');
    g = geoglobe(uif);
    geobasemap(g, 'satellite');
    
    max_index = height(flightdata);

    
    
    for k = 1:kstep:max_index
    
        tic
    
        xE =   flightdata.("xE [m]")(k);
        yE =   flightdata.("yE [m]")(k);
        h  = - flightdata.("zE [m]")(k) + reference.h;
    
        u = flightdata.("u [m/s]")(k);
        v = flightdata.("v [m/s]")(k);
        w = flightdata.("w [m/s]")(k);
    
    
        phi   = flightdata.("phi [rad]")(k);
        theta = flightdata.("theta [rad]")(k);
        psi   = flightdata.("psi [rad]")(k);
    
        R_BV = [cos(theta)*cos(psi), (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)), (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
                cos(theta)*sin(psi), (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)), (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
                        -sin(theta),                             sin(phi)*cos(theta),                            cos(phi)*cos(theta)];
    
        x_axis_E = R_VE*R_BV*x_axis_B;
        y_axis_E = R_VE*R_BV*y_axis_B;
        z_axis_E = R_VE*R_BV*z_axis_B;
        
        x_axis_vector = [xE, xE + x_axis_E(1), yE, yE + x_axis_E(2), h, h + x_axis_E(3)];
        y_axis_vector = [xE, xE + y_axis_E(1), yE, yE + y_axis_E(2), h, h + y_axis_E(3)];
        z_axis_vector = [xE, xE + z_axis_E(1), yE, yE + z_axis_E(2), h, h + z_axis_E(3)];
    
        
        set(position, 'XData', xE, 'YData', yE, 'ZData', h)
        set(x_axis, 'XData', [x_axis_vector(1), x_axis_vector(2)], 'YData', [x_axis_vector(3), x_axis_vector(4)], 'ZData', [x_axis_vector(5), x_axis_vector(6)])
        set(y_axis, 'XData', [y_axis_vector(1), y_axis_vector(2)], 'YData', [y_axis_vector(3), y_axis_vector(4)], 'ZData', [y_axis_vector(5), y_axis_vector(6)])
        set(z_axis, 'XData', [z_axis_vector(1), z_axis_vector(2)], 'YData', [z_axis_vector(3), z_axis_vector(4)], 'ZData', [z_axis_vector(5), z_axis_vector(6)])
        
    
        addpoints(trajectory, xE, yE, h);
        addpoints(ground_track, xE, yE, 0);
        
        
        
        if trel >= tstep
        
            plot3([x_axis_vector(1), x_axis_vector(2)], [x_axis_vector(3), x_axis_vector(4)], [x_axis_vector(5), x_axis_vector(6)], 'color', [1, 0, 0, 0.5], 'LineWidth', 2);
            plot3([y_axis_vector(1), y_axis_vector(2)], [y_axis_vector(3), y_axis_vector(4)], [y_axis_vector(5), y_axis_vector(6)], 'color', [0, 0, 1, 0.5], 'LineWidth', 2);
            plot3([z_axis_vector(1), z_axis_vector(2)], [z_axis_vector(3), z_axis_vector(4)], [z_axis_vector(5), z_axis_vector(6)], 'color', [0, 1, 0, 0.5], 'LineWidth', 2);
    
            trel = 0;
        end
    
        if trel_map >= tstep_map
    
            geoplot3(g, flightdata.("lat [deg]")(1:k), flightdata.("lon [deg]")(1:k), flightdata.("h [m]")(1:k), "y")
            
            trel_map = 0;
        end
    
        pause(wait_time)
    
    
        t = t + dt*kstep;
    
    
        loop_time = toc;
    
        wait_time = wait_time - K_clock*(loop_time - dt*kstep);
    
        axis equal
        grid on
    
        trel = trel + dt*kstep;
    
        trel_map = trel_map + dt*kstep;
    
    
    end

end
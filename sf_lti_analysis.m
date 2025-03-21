function result = sf_lti_analysis(A, B, type)
    labels_lateral      = {'v [m/s]', 'p [rad/s]', 'r [rad/s]', '\phi [rad]'};
    labels_longitudinal = {'u [m/s]', 'w [m/s]', 'q [rad/s]', '\theta [rad]'};
    
    [V, D] = eig(A);

    v1 = V(:,1);
    v2 = V(:,2);
    v3 = V(:,3);
    v4 = V(:,4);

    lam1 = D(1,1);
    lam2 = D(2,2);
    lam3 = D(3,3);
    lam4 = D(4,4);

    lambdas = [lam1, lam2, lam3, lam4];

    unf1  = abs(lam1);
    unf2  = abs(lam2);
    unf3  = abs(lam3);
    unf4  = abs(lam4);

    damp1 = -real(lam1)/unf1;
    damp2 = -real(lam2)/unf2;
    damp3 = -real(lam3)/unf3;
    damp4 = -real(lam4)/unf4;

    tau1 = 1/abs(real(lam1));
    tau2 = 1/abs(real(lam2));
    tau3 = 1/abs(real(lam3));
    tau4 = 1/abs(real(lam4));

    t1 = 0:0.01:approp_time(lam1);
    t2 = 0:0.01:approp_time(lam2);
    t3 = 0:0.01:approp_time(lam3);
    t4 = 0:0.01:approp_time(lam4);

    x1 = exp(lam1*t1).*v1; 
    x2 = exp(lam2*t2).*v2;
    x3 = exp(lam3*t3).*v3; 
    x4 = exp(lam4*t4).*v4;

    figure('Name', 'Modal analysis', 'Color', 'white')

    % put text above each column to convey mode dep. on lateral/long.

    for i = 1:4
        subplot(4,4,1+4*(i-1))
        plot(t1, real(x1(i,:)), 'k') % Black plot
        xlabel('t [s]')
        ylabel(labels_lateral(i))
        if strcmp(type, 'lateral')
            ylabel(labels_lateral(i))
        elseif strcmp(type, 'longitudinal')
            ylabel(labels_longitudinal(i))
        end
        ylim([-1 1])
        grid minor

        if i == 1
            title(strcat('Mode nr. 1: \lambda_1 = ', string(lam1)), 'FontSize', 10, 'FontWeight','normal')
        end
    end

    for i = 1:4
        subplot(4,4,2+4*(i-1))
        plot(t2, real(x2(i,:)), 'k') % Black plot
        xlabel('t [s]')
        if strcmp(type, 'lateral')
            ylabel(labels_lateral(i))
        elseif strcmp(type, 'longitudinal')
            ylabel(labels_longitudinal(i))
        end
        ylim([-1 1])
        grid minor

        if i == 1
            title(strcat('Mode nr. 2: \lambda_2 = ', string(lam2)), 'FontSize', 10, 'FontWeight','normal')
        end
    end

    for i = 1:4
        subplot(4,4,3+4*(i-1))
        plot(t3, real(x3(i,:)), 'k') % Black plot
        xlabel('t [s]')
        ylabel(labels_lateral(i))
        if strcmp(type, 'lateral')
            ylabel(labels_lateral(i))
        elseif strcmp(type, 'longitudinal')
            ylabel(labels_longitudinal(i))
        end
        ylim([-1 1])
        grid minor

        if i == 1
            title(strcat('Mode nr. 3: \lambda_3 = ', string(lam3)), 'FontSize', 10, 'FontWeight','normal')
        end
    end

    for i = 1:4
        subplot(4,4,4+4*(i-1))
        plot(t4, real(x4(i,:)), 'k') % Black plot
        xlabel('t [s]')
        ylabel(labels_lateral(i))
        if strcmp(type, 'lateral')
            ylabel(labels_lateral(i))
        elseif strcmp(type, 'longitudinal')
            ylabel(labels_longitudinal(i))
        end
        ylim([-1 1])
        grid minor

        if i == 1
            title(strcat('Mode nr. 4: \lambda_4 = ', string(lam4)), 'FontSize', 10, 'FontWeight','normal')
        end
    end
    
    % Pole-Zero map

    figure('Name',append('Pole-Zero Map: ', type), 'Color','white')
    hold on

    colors = {'r', 'g', 'b', 'm'};

    realdistinct = real(lambdas);
    
    for i=1:4
        lam = lambdas(i);
        re  = real(lam);
        im  = imag(lam);
        idx = find((realdistinct == re), 1, "first");
        plot(re, im, 'Marker','x', 'Color', colors{idx}, 'MarkerSize',12, 'LineWidth', 2)
    end

    xline(0, 'LineStyle',':')
    yline(0, 'LineStyle',':')
    axis equal
    grid on
    sgrid
    xlabel('Re')
    ylabel('Im')

    % Controllability & Stabilizability

    C = ctrb(A, B);
    
    controllable = (rank(C) == 4);

    unstable_indices = (diag(D) >= 0);

    Vunst   = V(:,unstable_indices);

    stabilizable = (rank([C, Vunst]) == rank(C));

    result.ctrb    = controllable;
    result.stbl    = stabilizable;
    result.modes   = V;
    result.lambda  = [lam1, lam2, lam3, lam4];
    result.nfrequency = [unf1, unf2, unf3, unf4];
    result.damping = [damp1, damp2, damp3, damp4];
    result.timeconstant = [tau1, tau2, tau3, tau4];
    
end

function tend = approp_time(lambda)
    % Rule of thumb: T = 5 / |Re(lambda)|
    if real(lambda) < 0
        tend = 5 / abs(real(lambda)); % Ensure decay is visible
    elseif real(lambda) > 0
        tend = 3 / real(lambda); % Prevent runaway instability
    else
        tend = 10; % Arbitrary for purely imaginary cases
    end
end


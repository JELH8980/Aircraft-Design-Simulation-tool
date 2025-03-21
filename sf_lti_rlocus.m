function [Kq, Kr, Kp, Kyar] = sf_lti_rlocus(Alon, Blon, Alat, Blat, Kq, Kr, Kp, Kyar)
% Initial figure setup
fig = figure('Name', 'Pole-Zero Map', 'Color', 'white');
hold on;

% Initial gains
if isempty(Kq)
    Kq   = 0;
end
if isempty(Kr)
    Kr   = 0;
end
if isempty(Kp)
    Kp   = 0;
end
if isempty(Kyar)
    Kyar   = 0;
end

% Initial Klat matrix (fixing duplicate definition)
Klat = [0,  0,   Kr, 0;
        0, Kp, Kyar, 0];

% Store initial plot handle
selected_gain = 0; % 0: none, 1: Kp, 2: Kr, 3: Kyar, 4: Kq

% Set callback functions
set(fig, 'KeyPressFcn', @keyPressCallback);
set(fig, 'WindowScrollWheelFcn', @scrollCallback);
set(fig, 'CloseRequestFcn', @closeCallback); % Add this to handle figure close

% Initial plot
update_longitudinal();
update_lateral();

% Wait for the figure to close (optional, depending on desired behavior)
uiwait(fig); % Pauses execution until figure is closed

% After figure closes, return the final values
% (This line is reached only if uiwait is used)

    function update_longitudinal()
        subplot(2,1,1)
        Klon = [0, 0, Kq, 0;
                0, 0,  0, 0];
        [~, D] = eig(Alon + Blon*Klon);
        lambdas = diag(D);
        colors = {'r', 'g', 'b', 'm'};
        realdistinct = real(lambdas);

        for i = 1:4
            lam = lambdas(i);
            re = real(lam);
            im = imag(lam);
            idx = find(realdistinct == re, 1, "first");
            plot(re, im, 'Marker', 'x', ...
                         'Color', colors{idx}, ...
                         'MarkerSize', 12, ...
                         'LineWidth', 2);
            if i == 1
                hold on
            end
        end
        hold off
        title(sprintf('Kq=%.2f', Kq));
        grid on;
        xlabel('Re');
        ylabel('Im');
        axis equal
        grid on
        sgrid
        drawnow;
    end

    function update_lateral()
        subplot(2,1,2)
        Klat = [0,  0,   Kr, 0;
                0, Kp, Kyar, 0];
        [~, D] = eig(Alat + Blat*Klat);
        lambdas = diag(D);
        colors = {'r', 'g', 'b', 'm'};
        realdistinct = real(lambdas);

        for i = 1:4
            lam = lambdas(i);
            re = real(lam);
            im = imag(lam);
            idx = find(realdistinct == re, 1, "first");
            plot(re, im, 'Marker', 'x', ...
                         'Color', colors{idx}, ...
                         'MarkerSize', 12, ...
                         'LineWidth', 2);
            if i == 1
                hold on
            end
        end
        hold off
        title(sprintf('Kp=%.2f, Kr=%.2f, Kyar=%.2f', Kp, Kr, Kyar));
        grid on;
        xlabel('Re');
        ylabel('Im');
        axis equal
        grid on
        sgrid
        drawnow;
    end

    function keyPressCallback(~, event)
        switch event.Key
            case '1'
                selected_gain = 1; % Kp
                disp('Selected Kp');
            case '2'
                selected_gain = 2; % Kr
                disp('Selected Kr');
            case '3'
                selected_gain = 3; % Kyar
                disp('Selected Kyar');
            case '4'
                selected_gain = 4; % Kq
                disp('Selected Kq');
            case 'return'
                close(fig); % This will trigger closeCallback
            otherwise
                selected_gain = 0; % None
                disp('No gain selected');
        end
    end

    function scrollCallback(~, event)
        if selected_gain > 0
            if event.VerticalScrollCount > 0 % Scroll down
                scale = -0.05;
            else % Scroll up
                scale = 0.05;
            end
            switch selected_gain
                case 1 % Kp
                    Kp = Kp + scale*(1+abs(Kp));
                    update_lateral();
                case 2 % Kr
                    Kr = Kr + scale*(1+abs(Kr));
                    update_lateral();
                case 3 % Kyar
                    Kyar = Kyar + scale*(1+abs(Kyar));
                    update_lateral();
                case 4 % Kq
                    Kq = Kq + scale*(1+abs(Kq));
                    update_longitudinal();
            end
        end
    end

    function closeCallback(~, ~)
        % When figure is closed, delete it and resume execution
        delete(fig);
        % Optionally: uiresume(fig) if using uiwait
    end
end


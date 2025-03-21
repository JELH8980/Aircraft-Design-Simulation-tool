function signal = generate_doublet(tmax, dt, t_start, t_mid, t_end, amplitude, display)
% GENERATE_DOUBLET - Function to generate a doublet signal as a time series.
%
% This function creates a doublet signal—a positive pulse followed by a negative 
% pulse of equal amplitude—over a specified time range. The signal is defined by 
% a start time, midpoint, and end time, and is returned as a MATLAB timeseries 
% object. An optional plot can be displayed to visualize the signal.
%
% INPUTS:
%   tmax           - Maximum simulation time [s].
%   dt             - Time step [s].
%   t_start        - Start time of the positive pulse [s].
%   t_mid          - Midpoint time where the signal switches from positive to 
%                    negative [s].
%   t_end          - End time of the negative pulse [s].
%   amplitude      - Magnitude of the positive and negative pulses (dimensionless).
%   display        - Boolean flag: true to plot the signal, false to skip plotting.
%
% OUTPUT:
%   signal         - A timeseries object containing the doublet signal values and 
%                    corresponding time vector.
%
% FUNCTIONALITY:
% - Generates a time vector `t` from 0 to tmax-dt with step size dt.
% - Initializes a zero signal array of the same length as `t`.
% - Constructs the doublet:
%   - Positive pulse: amplitude from t_start to t_mid.
%   - Negative pulse: -amplitude from t_mid to t_end.
%   - Zero elsewhere.
% - Converts the signal array and time vector into a timeseries object.
% - If display is true:
%   - Creates a figure titled "Time-series 'generate_doublet'" with a white background.
%   - Plots the signal versus time with a black line, minor grid, and labeled axes 
%     (t [s], s [-]).
%
% NOTES:
% - The signal is dimensionless (indicated by 's [-]' in the plot), suitable for 
%   applications like control inputs or perturbations.
% - Time inputs (tmax, dt, t_start, t_mid, t_end) must be in seconds and satisfy 
%   0 ≤ t_start < t_mid < t_end ≤ tmax; no validation is performed.
% - The plot is not returned; it remains open in the MATLAB environment.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

t = 0:dt:tmax-dt;  
signal = zeros(size(t));

for i = 1:length(t)
    % Positive pulse (rising then falling)
    if t(i) >= t_start && t(i) <t_mid
        signal(i) = amplitude;
    elseif t(i) >= t_mid && t(i) <= t_end
        signal(i) = -amplitude;
    end
end

% Convert to timeseries object
signal = timeseries(signal, t);

if display
    figure('Name', 'Time-series "generate_doublet"', 'Color', 'white');
    plot(t, signal.Data(:), 'Color', 'k');
    xlabel('t  [s]');
    ylabel('s [-]');
    grid minor;
end

end
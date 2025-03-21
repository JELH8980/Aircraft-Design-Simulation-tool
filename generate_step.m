function signal = generate_step(tmax, dt, t_start, t_end, amplitude, display)
% GENERATE_STEP - Function to generate a step signal as a time series.
%
% This function creates a step signal—a constant amplitude applied over a specified 
% time interval—within a given time range. The signal rises to the specified amplitude 
% at t_start and falls back to zero after t_end, and is returned as a MATLAB timeseries 
% object. An optional plot can be displayed to visualize the signal.
%
% INPUTS:
%   tmax           - Maximum simulation time [s].
%   dt             - Time step [s].
%   t_start        - Start time of the step [s].
%   t_end          - End time of the step [s].
%   amplitude      - Magnitude of the step (dimensionless).
%   display        - Boolean flag: true to plot the signal, false to skip plotting.
%
% OUTPUT:
%   signal         - A timeseries object containing the step signal values and 
%                    corresponding time vector.
%
% FUNCTIONALITY:
% - Generates a time vector `t` from 0 to tmax-dt with step size dt.
% - Initializes a zero signal array of the same length as `t`.
% - Constructs the step:
%   - Sets the signal to amplitude for t_start ≤ t ≤ t_end.
%   - Remains zero outside this interval.
% - Converts the signal array and time vector into a timeseries object.
% - If display is true:
%   - Creates a figure titled "Time-series 'generate_step'" with a white background.
%   - Plots the signal versus time with a black line, minor grid, and labeled axes 
%     (t [s], s [-]).
%
% NOTES:
% - The signal is dimensionless (indicated by 's [-]' in the plot), suitable for 
%   applications like control inputs or step response testing.
% - Time inputs (tmax, dt, t_start, t_end) must be in seconds and satisfy 
%   0 ≤ t_start ≤ t_end < tmax; no validation is performed.
% - The plot is not returned; it remains open in the MATLAB environment.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

t = 0:dt:tmax-dt;  
signal = zeros(size(t));

for i = 1:length(t)
    if t(i) >= t_start && t(i) <= t_end
        signal(i) = amplitude;
    end
end

% Convert to timeseries object
signal = timeseries(signal, t);

if display
    figure('Name', 'Time-series "generate_step"', 'Color', 'white');
    plot(t, signal.Data(:), 'Color', 'k');
    xlabel('t  [s]');
    ylabel('s [-]');
    grid minor;
end
end
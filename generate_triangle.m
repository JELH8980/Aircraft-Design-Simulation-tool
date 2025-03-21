function signal = generate_triangle(tmax, dt, t_start, t_peak, t_end, amplitude, display)
% GENERATE_TRIANGLE - Function to generate a triangular signal as a time series.
%
% This function creates a triangular signal—a linear rise to a peak amplitude followed 
% by a linear decline—over a specified time range. The signal starts at zero, rises 
% to the amplitude at t_peak, and falls back to zero by t_end. It is returned as a 
% MATLAB timeseries object, with an optional plot for visualization.
%
% INPUTS:
%   tmax           - Maximum simulation time [s].
%   dt             - Time step [s].
%   t_start        - Start time of the triangular rise [s].
%   t_peak         - Time at which the signal reaches its peak amplitude [s].
%   t_end          - End time of the triangular decline [s].
%   amplitude      - Peak magnitude of the triangle (dimensionless).
%   display        - Boolean flag: true to plot the signal, false to skip plotting.
%
% OUTPUT:
%   signal         - A timeseries object containing the triangular signal values and 
%                    corresponding time vector.
%
% FUNCTIONALITY:
% - Generates a time vector `t` from 0 to tmax-dt with step size dt.
% - Initializes a zero signal array of the same length as `t`.
% - Constructs the triangular signal:
%   - Linear rise: From 0 at t_start to amplitude at t_peak, scaled by time fraction.
%   - Linear decline: From amplitude at t_peak to 0 at t_end, scaled by time fraction.
%   - Zero outside the interval [t_start, t_end].
% - Converts the signal array and time vector into a timeseries object.
% - If display is true:
%   - Creates a figure titled "Time-series 'generate_triangle'" with a white background.
%   - Plots the signal versus time with a black line, minor grid, and labeled axes 
%     (t [s], s [-]).
%
% NOTES:
% - The signal is dimensionless (indicated by 's [-]' in the plot), suitable for 
%   applications like control inputs or testing transient responses.
% - Time inputs (tmax, dt, t_start, t_peak, t_end) must be in seconds and satisfy 
%   0 ≤ t_start < t_peak < t_end ≤ tmax; no validation is performed.
% - The plot is not returned; it remains open in the MATLAB environment.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

t = 0:dt:tmax-dt;  

signal = zeros(size(t));

for i = 1:length(t)
    if t(i) >= t_start && t(i) < t_peak
        signal(i) = amplitude * (t(i) - t_start) / (t_peak - t_start);
    elseif t(i) >= t_peak && t(i) <= t_end
        signal(i) = amplitude * (t_end - t(i)) / (t_end - t_peak);
    end
end

% Convert to timeseries object
signal = timeseries(signal, t);

if display
    figure('Name', 'Time-series "generate_triangle"', 'Color', 'white');
    plot(t, signal.Data(:), 'Color', 'k');
    xlabel('t  [s]');
    ylabel('s [-]');
    grid minor;
end 
end
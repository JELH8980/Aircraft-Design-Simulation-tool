function signal = generate_impulse(tmax, dt, t_center, amplitude, display)
% GENERATE_IMPULSE - Function to generate an impulse signal as a time series.
%
% This function creates an impulse signal—a single spike of specified amplitude 
% at a given time—over a time range. The impulse is approximated as a single non-zero 
% value at the time step closest to the specified center time, and the signal is 
% returned as a MATLAB timeseries object. An optional plot can be displayed to 
% visualize the signal.
%
% INPUTS:
%   tmax           - Maximum simulation time [s].
%   dt             - Time step [s].
%   t_center       - Time at which the impulse occurs [s].
%   amplitude      - Magnitude of the impulse (dimensionless).
%   display        - Boolean flag: true to plot the signal, false to skip plotting.
%
% OUTPUT:
%   signal         - A timeseries object containing the impulse signal values and 
%                    corresponding time vector.
%
% FUNCTIONALITY:
% - Generates a time vector `t` from 0 to tmax-dt with step size dt.
% - Initializes a zero signal array of the same length as `t`.
% - Identifies the index closest to t_center and sets the signal value at that 
%   index to the specified amplitude, creating a discrete impulse.
% - Converts the signal array and time vector into a timeseries object.
% - If display is true:
%   - Creates a figure titled "Time-series 'generate_impulse'" with a white background.
%   - Plots the signal versus time with a black line, minor grid, and labeled axes 
%     (t [s], s [-]).
%
% NOTES:
% - The signal is dimensionless (indicated by 's [-]' in the plot), suitable for 
%   applications like control inputs or perturbations.
% - The impulse is approximated as a single non-zero value; its width is determined 
%   by dt, and its exact timing depends on the nearest time step to t_center.
% - Time inputs (tmax, dt, t_center) must be in seconds with 0 ≤ t_center < tmax; 
%   no validation is performed.
% - The plot is not returned; it remains open in the MATLAB environment.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

t = 0:dt:tmax-dt;  
signal = zeros(size(t));

% Find the index closest to t_center
[~, idx] = min(abs(t - t_center));
signal(idx) = amplitude;

% Convert to timeseries object
signal = timeseries(signal, t);

if display
    figure('Name', 'Time-series "generate_impulse"', 'Color', 'white');
    plot(t, signal.Data(:), 'Color', 'k');
    xlabel('t  [s]');
    ylabel('s [-]');
    grid minor;
end

end
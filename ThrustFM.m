function FMThrust = ThrustFM(Tcom, M, h)
% THRUSTFM - Computes thrust forces and moments for a propulsion system.
%
% This function calculates the forces and moments generated by a propulsion system based on a thrust command, Mach number, 
% and altitude. It uses a specified engine model to determine the thrust magnitude, then resolves the thrust into body-frame 
% forces (X, Y, Z) and moments (L, M, N) based on the engine's position and orientation relative to the aircraft's reference point.
%
% INPUTS:
%   Tcom       - Thrust command (typically in percentage or normalized units, e.g., 0 to 100).
%   M          - Mach number (dimensionless).
%   h          - Altitude (meters).
%
% OUTPUTS:
%   FMThrust   - A 6-element vector of forces and moments in the body frame: [X, Y, Z, L, M, N]', where:
%                - X, Y, Z: forces in body frame (N)
%                - L, M, N: moments in body frame (N·m)
%
% FUNCTIONALITY:
% - Initializes persistent variables (engine parameters and model) on the first call from a 'model' struct in the base workspace.
% - Computes the thrust magnitude (T) using a dynamically loaded engine model function.
% - Resolves the thrust into forces (X, Z) based on the thrust angle (aT), with Y = 0 (assuming symmetry).
% - Calculates the pitching moment (M) based on the thrust vector and its lever arms (lT, hT) relative to the reference point (lB, hB).
% - Assumes zero rolling (L) and yawing (N) moments due to symmetry and alignment with the body axis.
% - Returns the combined force and moment vector for use in 6-DOF dynamics.
%
% NOTES:
% - Persistent variables (lT, hT, aT, lB, hB, Tmax, engine_model) are initialized once from 'model.parameters'.
% - The engine model is loaded as a function handle from 'model.parameters.propulsion.method' (e.g., 'engine_function.m').
% - Thrust angle (aT) is assumed to be in degrees and converted to radians for trigonometric calculations.
% - The function adds the 'aircraft\propulsion\' directory to the MATLAB path during initialization.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

    % Define persistent variables
    persistent initialized_flag lT hT aT lB hB Tmax engine_model

    % Check if initialization is needed
    if isempty(initialized_flag)
        % Expecting the second argument to be a struct called 'parameters'
        addpath("aircraft\propulsion\")
        
        
        model               = evalin("base", 'model');
        lT                  = model.parameters.lT;
        hT                  = model.parameters.hT;
        aT                  = model.parameters.aT;
        lB                  = model.parameters.lB;
        hB                  = model.parameters.hB;
        Tmax                =  model.parameters.propulsion.Tmax;
        engine_model        = str2func(strrep(model.parameters.propulsion.method, '.m', ''));  % Convert function name to handle

        initialized_flag    = 1;
       
    end

    T = engine_model(Tmax, Tcom, M, h);

    X =  T*cosd(aT);
    Y =  0;
    Z = -T*sind(aT);
    
    
    L = 0;
    M = cosd(aT)*(hB-hT) + sind(aT)*(lB-lT)*T;
    N = 0;

    FMThrust = [X, Y, Z, L, M, N]';

end
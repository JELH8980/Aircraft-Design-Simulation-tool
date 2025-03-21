function [T] = MilEngine(Tmax,Tcom,M,Hp_m)
%
% Function:
% ---------
% A generic model of thrust for a single military aircraft engine. By
% giving the maximum thrust for one engine of the project this can be used.
% The comanded thrust setting is used as a throttle from 0-100% comanded
% thrust. The data is given for different Mach (0.0-0.9) and pressure
% altitudes (0-40kft). The input altitide, given in meters is converted
% into kft to use the data. The output is given in Newton.
% 
% Input:
% ------
%   Tmax    Max Thrust (N)
%   Tcom    Comanded Thrust (% of Tmax, 0-100%)
%   M       Mach number ( )
%   Hp_m    Pressure altitude (m)
%
% Output:
% -------
%   T   Thrust (N)
%
% Version:
% -------
% 1. 2025-01-27 First version by Roger Larsson.
%

%--------------------------------------------------------------------------
% Code
%--------------------------------------------------------------------------
%
% Constants
%
m2kft=3.281; % Meter-2-kfeet.
%
% Data
%
% Mach vector in data, unit (   ).
EngineMach = repmat([0.0 0.2 0.4 0.6 0.8 1.0],6,1);
%
% Pressure altitude in data, unit (kft).
EngineHp = repmat([0 10 20 30 40 50]',1,6);

%
% Thrust matrix in data T/Tmax for Tcom at 100%, unit (   ).
EngineData = [...
    1.0000 1.0000 0.9945 0.9968 0.9771 0.9211
    0.7216 0.7216 0.7344 0.7759 0.8025 0.7767
    0.4890 0.4979 0.5213 0.5591 0.6112 0.6349
    0.3115 0.3186 0.3383 0.3675 0.4196 0.4811
    0.1932 0.1948 0.2050 0.2240 0.2563 0.2997
    0.1104 0.1104 0.1230 0.1309 0.1522 0.1822];
%
% Converting Hp from meter to kfeet
Hp_kft = m2kft*Hp_m/1000;
%
% Thrust interpolation in general thrust data
Tgen = interp2(EngineMach,EngineHp,EngineData,M,Hp_kft);
%
% Trust output for specific engine and comanded thrust setting
T = Tgen*Tmax*(Tcom/100);
%
% End of MilEngine
%
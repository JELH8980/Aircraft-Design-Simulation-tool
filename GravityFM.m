function FGravity = GravityFM(e1, e2)
% GRAVITYFM - Computes gravitational force and moment vector in body frame.
%
% Calculates the gravitational force vector (XG, YG, ZG) in the body frame given 
% roll (e1) and pitch (e2) angles, using mass (m) and gravity (g) from the base 
% workspace. Returns a 6x1 vector with forces and zero moments.
%
% INPUTS: e1 - Roll angle [rad], e2 - Pitch angle [rad].
% OUTPUT: FGravity - [XG, YG, ZG, 0, 0, 0]' [N, N, N, N*m, N*m, N*m].
%
% Author: Ludwig Horvath
% Date: 3/17/2025


 % Define persistent variables
    persistent initialized_flag m g


    if isempty(initialized_flag)
        
        model = evalin('base', 'model');
        reference = evalin('base', 'reference');
        m                   = model.parameters.m;
        g                   = reference.g;
        initialized_flag    = 1;
       
    end

    XG = -sin(e2) * m * g;            % Forward (positive when nose-up)
    YG =  sin(e1) * cos(e2) * m * g;  % Right (negative when rolled right)
    ZG =  cos(e1) * cos(e2) * m * g;  % Down (positive downward)

    FGravity = [XG, YG, ZG, 0, 0, 0]';

end
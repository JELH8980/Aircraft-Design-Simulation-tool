function [model] = calculate_inertia_1(model)
% CALCULATE_parameters_1 - Prompt user for parameters properties and update model.
%
% This function interacts with the user to manually input mass, moments of 
% parameters, and center of gravity for a given model. The values are stored 
% in the model's `parameters` field.
%
% INPUT:
%   model - A struct containing model properties.
%
% OUTPUT:
%   model - Updated struct with calculated parameters values.
%
% FUNCTIONALITY:
% - Requests user input for mass, parameters values, and center of gravity.
% - Converts input strings to numerical values.
% - Stores the parameters properties in `model.parameters`.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

model.parameters.m = input('    Enter mass: ');

model.parameters.Ixx = input('    Enter moment of parameters Ixx: ');
model.parameters.Iyy = input('    Enter moment of parameters Iyy: ');
model.parameters.Izz = input('    Enter moment of parameters Izz: ');
model.parameters.Ixz = input('    Enter moment of parameters Ixz: ');

model.parameters.lB = input('    Enter center of gravity in X (lB): ');
model.parameters.hB = input('    Enter center of gravity in Z (hB): ');

model.parameters.method.inertia = 'calculate_inertia_1';

model.geo.CG(1) = model.parameters.lB - model.parameters.lR;
model.geo.CG(2) = 0;
model.geo.CG(3) = model.parameters.lB - model.parameters.lB;

end
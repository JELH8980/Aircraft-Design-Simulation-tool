function interpolant = create_interpolant(batch)
% CREATE_INTERPOLANT - Function to create interpolation models for aerodynamic coefficients.
%
% This function generates a struct of gridded interpolants for non-dimensional 
% aerodynamic force and moment coefficients (CD, CC, CL, Cl, Cm, Cn) and their 
% derivatives, based on batch data. It uses alpha (angle of attack) and delta_e 
% (elevator deflection) as input grids, converting them from degrees to radians, 
% and reshapes coefficient data into 2D grids for interpolation.
%
% INPUTS:
%   batch          - A struct containing batch information (batch.info) and data 
%                    (batch.data). batch.info includes alpha.grid and delta_e.grid 
%                    (in degrees); batch.data includes coefficient arrays (CD, CC, 
%                    CL, Cl, Cm, Cn, and their derivatives).
%
% OUTPUT:
%   interpolant    - A struct containing griddedInterpolant objects for each 
%                    aerodynamic coefficient and its derivatives, plus metadata 
%                    (id, input labels, and comment).
%
% FUNCTIONALITY:
% - Extracts data and info from the batch struct.
% - Creates 2D grids for alpha and delta_e (converted to radians) using ndgrid.
% - Reshapes coefficient arrays (e.g., CD, CC, CL) and their derivatives (e.g., 
%   CD_a_dot, CC_b) into the grid size defined by alpha and delta_e.
% - Builds griddedInterpolant objects for:
%   - Primary coefficients: CD, CC, CL, Cl, Cm, Cn.
%   - Derivative coefficients with respect to alpha_dot, beta, beta_dot, P, Q, R, 
%     delta_a, delta_r (where applicable).
% - Adds metadata: batch ID, input labels ('alpha [rad]', 'delta_e [rad]'), and a 
%   comment describing the interpolant’s purpose.
%
% NOTES:
% - Assumes batch.data contains flat arrays for coefficients and derivatives, 
%   with lengths matching the product of alpha.grid and delta_e.grid sizes.
% - Input grids (alpha.grid, delta_e.grid) are in degrees and converted to radians 
%   for interpolation.
% - Uses MATLAB’s `griddedInterpolant` for efficient 2D interpolation; default 
%   method is linear unless specified otherwise.
% - CC may represent a custom coefficient (e.g., side force coefficient CY); its 
%   exact meaning depends on the batch data context.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

data = batch.data;

info = batch.info;

[alpha_grid, delta_e_grid] = ndgrid(deg2rad(info.alpha.grid), deg2rad(info.delta_e.grid));

grid_size = size(delta_e_grid);

CD_grid = reshape(data.CD, grid_size);
interpolant.CD = griddedInterpolant(alpha_grid, delta_e_grid, CD_grid);

CC_grid = reshape(data.CC, grid_size);
interpolant.CC = griddedInterpolant(alpha_grid, delta_e_grid, CC_grid);

CL_grid = reshape(data.CL, grid_size);
interpolant.CL = griddedInterpolant(alpha_grid, delta_e_grid, CL_grid);

Cl_grid = reshape(data.Cl, grid_size);
interpolant.Cl = griddedInterpolant(alpha_grid, delta_e_grid, Cl_grid);

Cm_grid = reshape(data.Cm, grid_size);
interpolant.Cm = griddedInterpolant(alpha_grid, delta_e_grid, Cm_grid);

Cn_grid = reshape(data.Cn, grid_size);
interpolant.Cn = griddedInterpolant(alpha_grid, delta_e_grid, Cn_grid);

% Additional CD coefficients
CD_a_dot_grid = reshape(data.CD_a_dot, grid_size);
interpolant.CD_a_dot = griddedInterpolant(alpha_grid, delta_e_grid, CD_a_dot_grid);

CD_b_grid = reshape(data.CD_b, grid_size);
interpolant.CD_b = griddedInterpolant(alpha_grid, delta_e_grid, CD_b_grid);

CD_P_grid = reshape(data.CD_P, grid_size);
interpolant.CD_P = griddedInterpolant(alpha_grid, delta_e_grid, CD_P_grid);

CD_Q_grid = reshape(data.CD_Q, grid_size);
interpolant.CD_Q = griddedInterpolant(alpha_grid, delta_e_grid, CD_Q_grid);

CD_R_grid = reshape(data.CD_R, grid_size);
interpolant.CD_R = griddedInterpolant(alpha_grid, delta_e_grid, CD_R_grid);

CD_delta_a_grid = reshape(data.CD_delta_a, grid_size);
interpolant.CD_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, CD_delta_a_grid);

CD_delta_r_grid = reshape(data.CD_delta_r, grid_size);
interpolant.CD_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, CD_delta_r_grid);

% Additional CC coefficients
CC_b_grid = reshape(data.CC_b, grid_size);
interpolant.CC_b = griddedInterpolant(alpha_grid, delta_e_grid, CC_b_grid);

CC_b_dot_grid = reshape(data.CC_b_dot, grid_size);
interpolant.CC_b_dot = griddedInterpolant(alpha_grid, delta_e_grid, CC_b_dot_grid);

CC_P_grid = reshape(data.CC_P, grid_size);
interpolant.CC_P = griddedInterpolant(alpha_grid, delta_e_grid, CC_P_grid);

CC_Q_grid = reshape(data.CC_Q, grid_size);
interpolant.CC_Q = griddedInterpolant(alpha_grid, delta_e_grid, CC_Q_grid);

CC_R_grid = reshape(data.CC_R, grid_size);
interpolant.CC_R = griddedInterpolant(alpha_grid, delta_e_grid, CC_R_grid);

CC_delta_a_grid = reshape(data.CC_delta_a, grid_size);
interpolant.CC_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, CC_delta_a_grid);

CC_delta_r_grid = reshape(data.CC_delta_r, grid_size);
interpolant.CC_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, CC_delta_r_grid);

% Additional CL coefficients
CL_a_dot_grid = reshape(data.CL_a_dot, grid_size);
interpolant.CL_a_dot = griddedInterpolant(alpha_grid, delta_e_grid, CL_a_dot_grid);

CL_b_grid = reshape(data.CL_b, grid_size);
interpolant.CL_b = griddedInterpolant(alpha_grid, delta_e_grid, CL_b_grid);

CL_P_grid = reshape(data.CL_P, grid_size);
interpolant.CL_P = griddedInterpolant(alpha_grid, delta_e_grid, CL_P_grid);

CL_Q_grid = reshape(data.CL_Q, grid_size);
interpolant.CL_Q = griddedInterpolant(alpha_grid, delta_e_grid, CL_Q_grid);

CL_R_grid = reshape(data.CL_R, grid_size);
interpolant.CL_R = griddedInterpolant(alpha_grid, delta_e_grid, CL_R_grid);

CL_delta_a_grid = reshape(data.CL_delta_a, grid_size);
interpolant.CL_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, CL_delta_a_grid);

CL_delta_r_grid = reshape(data.CL_delta_r, grid_size);
interpolant.CL_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, CL_delta_r_grid);

% Additional Cl coefficients
Cl_b_grid = reshape(data.Cl_b, grid_size);
interpolant.Cl_b = griddedInterpolant(alpha_grid, delta_e_grid, Cl_b_grid);

Cl_P_grid = reshape(data.Cl_P, grid_size);
interpolant.Cl_P = griddedInterpolant(alpha_grid, delta_e_grid, Cl_P_grid);

Cl_Q_grid = reshape(data.Cl_Q, grid_size);
interpolant.Cl_Q = griddedInterpolant(alpha_grid, delta_e_grid, Cl_Q_grid);

Cl_R_grid = reshape(data.Cl_R, grid_size);
interpolant.Cl_R = griddedInterpolant(alpha_grid, delta_e_grid, Cl_R_grid);

Cl_delta_a_grid = reshape(data.Cl_delta_a, grid_size);
interpolant.Cl_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, Cl_delta_a_grid);

Cl_delta_r_grid = reshape(data.Cl_delta_r, grid_size);
interpolant.Cl_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, Cl_delta_r_grid);

% Additional Cm coefficients
Cm_a_dot_grid = reshape(data.Cm_a_dot, grid_size);
interpolant.Cm_a_dot = griddedInterpolant(alpha_grid, delta_e_grid, Cm_a_dot_grid);

Cm_b_grid = reshape(data.Cm_b, grid_size);
interpolant.Cm_b = griddedInterpolant(alpha_grid, delta_e_grid, Cm_b_grid);

Cm_P_grid = reshape(data.Cm_P, grid_size);
interpolant.Cm_P = griddedInterpolant(alpha_grid, delta_e_grid, Cm_P_grid);

Cm_Q_grid = reshape(data.Cm_Q, grid_size);
interpolant.Cm_Q = griddedInterpolant(alpha_grid, delta_e_grid, Cm_Q_grid);

Cm_R_grid = reshape(data.Cm_R, grid_size);
interpolant.Cm_R = griddedInterpolant(alpha_grid, delta_e_grid, Cm_R_grid);

Cm_delta_a_grid = reshape(data.Cm_delta_a, grid_size);
interpolant.Cm_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, Cm_delta_a_grid);

Cm_delta_r_grid = reshape(data.Cm_delta_r, grid_size);
interpolant.Cm_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, Cm_delta_r_grid);

% Additional Cn coefficients
Cn_a_grid = reshape(data.Cn_a, grid_size);
interpolant.Cn_a = griddedInterpolant(alpha_grid, delta_e_grid, Cn_a_grid);

Cn_b_grid = reshape(data.Cn_b, grid_size);
interpolant.Cn_b = griddedInterpolant(alpha_grid, delta_e_grid, Cn_b_grid);

Cn_b_dot_grid = reshape(data.Cn_b_dot, grid_size);
interpolant.Cn_b_dot = griddedInterpolant(alpha_grid, delta_e_grid, Cn_b_dot_grid);

Cn_P_grid = reshape(data.Cn_P, grid_size);
interpolant.Cn_P = griddedInterpolant(alpha_grid, delta_e_grid, Cn_P_grid);

Cn_Q_grid = reshape(data.Cn_Q, grid_size);
interpolant.Cn_Q = griddedInterpolant(alpha_grid, delta_e_grid, Cn_Q_grid);

Cn_R_grid = reshape(data.Cn_R, grid_size);
interpolant.Cn_R = griddedInterpolant(alpha_grid, delta_e_grid, Cn_R_grid);

Cn_delta_a_grid = reshape(data.Cn_delta_a, grid_size);
interpolant.Cn_delta_a = griddedInterpolant(alpha_grid, delta_e_grid, Cn_delta_a_grid);

Cn_delta_r_grid = reshape(data.Cn_delta_r, grid_size);
interpolant.Cn_delta_r = griddedInterpolant(alpha_grid, delta_e_grid, Cn_delta_r_grid);


interpolant.id = batch.info.id;

interpolant.input = {'alpha [rad]', 'delta_e [rad]'};

interpolant.comment = 'returns non-dimensional aerodynamic coefficients as function of inputs.';

end

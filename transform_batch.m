function batch = transform_batch(batch)
% TRANSFORM_BATCH - Transforms aerodynamic coefficients from wind frame to body frame for a batch of data.
%
% This function converts aerodynamic coefficients (drag, side force, lift) and their derivatives from the wind frame to the 
% body frame based on the angle of attack (alpha) and sideslip angle (beta). It processes each row of the batch data table, 
% applying a rotation transformation to the coefficients, and renames the variables (CD to CX, CC to CY, CL to CZ) to reflect 
% the body-frame convention. The transformed data is then stored back into the batch struct.
%
% INPUTS:
%   batch      - A struct containing a data table:
%                - batch.data: A table with columns including 'alpha', 'beta', 'CD', 'CC', 'CL', and their derivatives 
%                              (e.g., 'CD_a', 'CL_delta_e'), where angles are in degrees and coefficients are in wind frame.
%
% OUTPUTS:
%   batch      - The updated struct with transformed data:
%                - batch.data: The table with coefficients transformed to body frame and renamed (e.g., 'CX', 'CY', 'CZ').
%
% FUNCTIONALITY:
% - Checks for the existence of 'alpha_dot' in the table to conditionally transform related derivatives.
% - Iterates through each row of the data table, extracting alpha and beta.
% - Applies a 3D rotation transformation to convert wind-frame coefficients (CD, CC, CL) to body-frame coefficients (CX, CY, CZ).
% - Transforms derivatives with respect to state variables (e.g., alpha, beta, P, Q, R) and control inputs (e.g., delta_a, delta_e, delta_r).
% - Renames table headers to reflect body-frame notation (e.g., 'CD' to 'CX', 'CC' to 'CY', 'CL' to 'CZ').
% - Updates the batch struct with the transformed data table.
%
% NOTES:
% - Angles (alpha, beta) are assumed to be in degrees and are converted to radians using cosd and sind for trigonometric calculations.
% - The transformation assumes a standard aerospace coordinate system where:
%   - CD (drag) is along the wind velocity vector.
%   - CC (side force) is perpendicular to the lift plane in the wind frame.
%   - CL (lift) is normal to the velocity vector in the wind frame.
% - Body-frame coefficients (CX, CY, CZ) align with the aircraft's longitudinal, lateral, and vertical axes, respectively.
% - Derivatives with respect to alpha_dot and beta_dot are conditionally processed based on their presence in the table.
%
% Author: Ludwig Horvath
% Date: 2/11/2025

    T = batch.data;
    
    alpha_dot_existenceflag = any(ismember(T.Properties.VariableNames, 'alpha_dot'));

    for i = 1:height(T)
    
        t = T(i,:);
       
    
        alpha           = t.alpha;
        beta            = t.beta;
        
        T.CD(i)         = -cosd(alpha)*cosd(beta)*t.CD - cosd(alpha)*sind(beta)*t.CC + sind(alpha)*t.CL;
        T.CC(i)         = -sind(beta)*t.CD + cosd(beta)*t.CC;
        T.CL(i)         = -sind(alpha)*cosd(beta)*t.CD - sind(alpha)*sind(beta)*t.CC - cosd(alpha)*t.CL;
        
        T.CD_a(i)       = -cosd(alpha)*cosd(beta)*t.CD_a - cosd(alpha)*sind(beta)*t.CC_a + sind(alpha)*t.CL_a;
        T.CC_a(i)       = -sind(beta)*t.CD_a + cosd(beta)*t.CC_a;
        T.CL_a(i)       = -sind(alpha)*cosd(beta)*t.CD_a - sind(alpha)*sind(beta)*t.CC_a - cosd(alpha)*t.CL_a;
        
        if alpha_dot_existenceflag
            T.CD_a_dot(i)   = -cosd(alpha)*cosd(beta)*t.CD_a_dot + sind(alpha)*t.CL_a_dot;
            T.CL_a_dot(i)   = -sind(alpha)*cosd(beta)*t.CD_a_dot - cosd(alpha)*t.CL_a_dot;
        end

        T.CD_b(i)       = -cosd(alpha)*cosd(beta)*t.CD_b - cosd(alpha)*sind(beta)*t.CC_b + sind(alpha)*t.CL_b;
        T.CC_b(i)       = -sind(beta)*t.CD_b + cosd(beta)*t.CC_b;
        T.CL_b(i)       = -sind(alpha)*cosd(beta)*t.CD_b - sind(alpha)*sind(beta)*t.CC_b - cosd(alpha)*t.CL_b;

        T.CD_P(i)       = -cosd(alpha)*cosd(beta)*t.CD_P - cosd(alpha)*sind(beta)*t.CC_P + sind(alpha)*t.CL_P;
        T.CC_P(i)       = -sind(beta)*t.CD_P + cosd(beta)*t.CC_P;
        T.CL_P(i)       = -sind(alpha)*cosd(beta)*t.CD_P - sind(alpha)*sind(beta)*t.CC_P - cosd(alpha)*t.CL_P;
        
        T.CD_Q(i)       = -cosd(alpha)*cosd(beta)*t.CD_Q - cosd(alpha)*sind(beta)*t.CC_Q + sind(alpha)*t.CL_Q;
        T.CC_Q(i)       = -sind(beta)*t.CD_Q + cosd(beta)*t.CC_Q;
        T.CL_Q(i)       = -sind(alpha)*cosd(beta)*t.CD_Q - sind(alpha)*sind(beta)*t.CC_Q - cosd(alpha)*t.CL_Q;
        
        T.CD_R(i)       = -cosd(alpha)*cosd(beta)*t.CD_R - cosd(alpha)*sind(beta)*t.CC_R + sind(alpha)*t.CL_R;
        T.CC_R(i)       = -sind(beta)*t.CD_R + cosd(beta)*t.CC_R;
        T.CL_R(i)       = -sind(alpha)*cosd(beta)*t.CD_R - sind(alpha)*sind(beta)*t.CC_R - cosd(alpha)*t.CL_R;
        
        T.CD_delta_a(i) = -cosd(alpha)*cosd(beta)*t.CD_delta_a - cosd(alpha)*sind(beta)*t.CC_delta_a + sind(alpha)*t.CL_delta_a;
        T.CC_delta_a(i) = -sind(beta)*t.CD_delta_a + cosd(beta)*t.CC_delta_a;
        T.CL_delta_a(i) = -sind(alpha)*cosd(beta)*t.CD_delta_a - sind(alpha)*sind(beta)*t.CC_delta_a - cosd(alpha)*t.CL_delta_a;
        
        T.CD_delta_e(i) = -cosd(alpha)*cosd(beta)*t.CD_delta_e - cosd(alpha)*sind(beta)*t.CC_delta_e + sind(alpha)*t.CL_delta_e;
        T.CC_delta_e(i) = -sind(beta)*t.CD_delta_e + cosd(beta)*t.CC_delta_e;
        T.CL_delta_e(i) = -sind(alpha)*cosd(beta)*t.CD_delta_e - sind(alpha)*sind(beta)*t.CC_delta_e - cosd(alpha)*t.CL_delta_e;
        
        T.CD_delta_r(i) = -cosd(alpha)*cosd(beta)*t.CD_delta_r - cosd(alpha)*sind(beta)*t.CC_delta_r + sind(alpha)*t.CL_delta_r;
        T.CC_delta_r(i) = -sind(beta)*t.CD_delta_r + cosd(beta)*t.CC_delta_r;
        T.CL_delta_r(i) = -sind(alpha)*cosd(beta)*t.CD_delta_r - sind(alpha)*sind(beta)*t.CC_delta_r - cosd(alpha)*t.CL_delta_r;
    
    
    end
    %%
    
    headers = T.Properties.VariableNames;
    
    updated_headers = strrep(headers, 'CD', 'CX');
    updated_headers = strrep(updated_headers, 'CC', 'CY');
    updated_headers = strrep(updated_headers, 'CL', 'CZ');
    
    T.Properties.VariableNames = updated_headers;

    batch.data = T;

end


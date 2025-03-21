function [model] = assign_patches(model)
% ASSIGN_PATCHES - Function to assign patches to regions of the model's geometry.
%
% This function prompts the user to define and assign patches to specified regions of the model.
% The user can specify the symmetry type, the number of patches, and the patch locations within 
% the geometry. The patch details are stored in the `model.geo.patches` structure for further 
% analysis or processing.
%
% INPUTS:
%   model          - The model struct containing the geometry (model.geo), where the patches will 
%                    be assigned. The model must already contain the `geo` field.
%
% OUTPUTS:
%   model          - The updated model struct with the assigned patches in the `geo.patches` field.
%
% FUNCTIONALITY:
% - The function initializes the `model.geo.patches` structure if it doesn't exist.
% - The user is prompted with options to define regions (symmetry and number of patches).
% - The user specifies patch numbers, and the corresponding coordinates are extracted from the 
%   mesh to define each patch's position in 3D space.
% - A figure window is opened to display the geometry and facilitate easy selection of patches.
% - Error handling is provided for invalid patch numbers, and the process repeats until the user 
%   selects the option to close the figure and exit.
%
% NOTES:
% - This function allows the user to define the number of patches and their associated symmetry 
%   types within each region. 
% - The patch data, including coordinates and symmetry information, is saved in the model for 
%   future use or analysis.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


try 
    
    if ~isfield(model.geo, 'patches')
        model.geo.patches = struct();
    end
    
    % Changing variables to plot only partition outline
    g = model.geo;
    g.nx = double(g.nx > 0);
    g.ny = double(g.ny > 0);
    g.fnx = double(g.fnx > 0);
    s.AS = 1;
    s.alpha = 0;
    s.betha = 0;
    s.P = 0;
    s.Q = 0;
    s.R = 0;
    s.ALT = 0;
    s.rho = 1;
    s.pgcorr = 0;
    
    [mesh, ~] = fLattice_setup(g, s, 1);
    
    fig = display_geometry(model, true); % Enable labeling of patches for easy selection

    while true
    q=questions(5);
    
    if isempty(q)
       q=-1; %will go to otherwise section
    end
    
    
    switch (q)
        case 1
        
        region_name = questions(6);
         
        % Replace spaces with underscores for robustness
        region_name = regexprep(region_name, '\s+', '_');
    
        % Get symmetry type
        disp(append('  [1] Symmetric ', region_name))
        disp(append('  [0] Asymmetric ', region_name))
        symmetry_type = input('  Enter a choice from above: ');
    
        % Get number of patches
        region_prompt = append('  Enter number of patches in ', region_name, ': ');
        no_patches = input(region_prompt);
    
        % Initialize patch storage
        model.geo.patches.(region_name).N = no_patches;
        model.geo.patches.(region_name).frames = zeros(no_patches, 5, 3);
        model.geo.patches.(region_name).symmetry = symmetry_type;
    
        % Loop over patches
        for i = 1:no_patches
            patch_prompt = sprintf('  Input Patch Number (%d/%d): ', i, no_patches);
            patch_nr = input(patch_prompt);
    
            % Validate index
            if patch_nr < 1 || patch_nr > size(mesh.XYZ, 1)
                disp('  Invalid patch number, skipping.');
                continue;
            end
    
            % Extract coordinates
            X = mesh.XYZ(patch_nr,:,1);
            Y = mesh.XYZ(patch_nr,:,2); % Fix: was incorrectly using X
            Z = mesh.XYZ(patch_nr,:,3); % Fix: was incorrectly using X
    
            % Store in structure
            model.geo.patches.(region_name).frames(i,:,1) = X;
            model.geo.patches.(region_name).frames(i,:,2) = Y;
            model.geo.patches.(region_name).frames(i,:,3) = Z;
        end
    
        case 0
            close(fig);
            break;
    
    end

    end

catch
     % Input error handling here.
    disp('Error occured')


end

   
end

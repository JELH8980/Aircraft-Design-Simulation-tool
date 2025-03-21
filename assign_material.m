function [model] = assign_material(model)
% ASSIGN_MATERIAL - Function to assign and manage materials for bodies in the model.
%
% This function handles the assignment and management of materials for bodies within
% the geometry of the model. It provides a user interface to either add new materials, 
% view existing materials, or define a composition of materials for each body. It supports
% the addition of materials, modification of material properties (density and color), 
% and updating material fractions in case of multiple materials for a body.
%
% INPUTS:
%   model          - The model struct containing geometry and body data, including
%                    material properties to be assigned or modified.
%
% OUTPUT:
%   model          - The updated model with assigned or modified materials.
%
% FUNCTIONALITY:
% - The function displays a menu of options for the user to either add new materials,
%   modify material composition, or view existing material details.
% - New materials can be added with a name, density, and color.
% - For bodies with existing materials, their properties (density, color) are displayed,
%   and the user can modify material fractions interactively.
% - The user can cancel or exit the menu, returning control to the main program.
% - If an error occurs during execution, the `emergency_save` function is called to save
%   the current state of the model.
%
% NOTES:
% - The function loops, allowing the user to continue assigning materials until they 
%   choose to exit (option 0).
% - Material names are sanitized by replacing spaces with underscores.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


try
    while true
    q=questions(8);
    
    if isempty(q)
       q=-1;  %will go to otherwise section
    end
    
    
    switch (q)     
        case 1
            
            options = {};  % Initialize an empty cell array
    
            bodies_names = fields(model.geo.bodies);
            no_bodies = numel(bodies_names);
            
        
            for i = 1:no_bodies
                options = [options, bodies_names{i}];  % Append to the cell array
                disp(append('    [', string(i), ']. ', bodies_names{i}));        
            end
            
            disp('    [0]. Return');  

            disp(' ')

            choice = input('    Please enter choice from above: ');
            disp(' ')
        
            if choice == 0
                break;
            end

            if isfield(model.geo.bodies.(options{choice}), 'materials')
                materials_names = fields(model.geo.bodies.(options{choice}).materials);
                num_materials = numel(materials_names);

                disp('    Existing materials: ')
                for i = 1:num_materials
    
                    name =  materials_names{i};
                    color = model.geo.bodies.(options{choice}).materials.(name).color;
                    density = model.geo.bodies.(options{choice}).materials.(name).density;
    
                    disp(append('    - Name: ', name, ', Color: ', num2str(color), ', Density: ', num2str(density)));        
                end
                disp(' ')

            else
                disp('    No existing materials...')
                disp(' ')

            end
    
            name = input('    Enter additional material name: ', 's');

            name = regexprep(name, '\s+', '_');
            
            disp(' ')

            material = struct();
        
            density = str2double(input(append('    Enter ', name, ' density [kg/m^3]: '), 's'));
            
            material.density = density;

            color = input(append('    Enter a color for ', name, ': '), 's');
            
            material.color = color;

            model.geo.bodies.(options{choice}).materials.(name) = material;
    

        case 2
            options = {};  % Initialize an empty cell array
    
            bodies_names = fields(model.geo.bodies);
            no_bodies = numel(bodies_names);
            
        
            for i = 1:no_bodies
                options = [options, bodies_names{i}];  % Append to the cell array
                disp(append('    [', string(i), ']. ', bodies_names{i}));        
            end
            disp(' ')

            choice = input('    Please enter choice from above: ');
            disp(' ')


            if isfield(model.geo.bodies.(options{choice}), 'materials')
               materials_names = fields(model.geo.bodies.(options{choice}).materials);
               num_materials = numel(materials_names);
    
               disp('    Existing materials: ')
               for i = 1:num_materials
    
                   name =  materials_names{i};
                   color = model.geo.bodies.(options{choice}).materials.(name).color;
                   density = model.geo.bodies.(options{choice}).materials.(name).density;
    
                   disp(append('    - Name: ', name, ', Color: ', num2str(color), ', Density: ', num2str(density)));            
               end
               disp(' ')


               composition = interactive_material_slider(model.geo.bodies.(options{choice}).materials);
               
               for i = 1:num_materials
                    model.geo.bodies.(options{choice}).materials.( materials_names{i}).fraction = composition.(materials_names{i});
                end


            end

            

        case 0


            break;

            
    
    end
    
    end

catch
    % Handle errors here
    emergency_save(model)
end


end

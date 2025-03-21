function [model] = set_flaps(model, type, angle)
    
    flap_matrix = model.geo.flap_vector;

    % Convert angle to radians
    angle_rad = angle * pi / 180;
    
    % Logical indexing using lmatrix
    lmatrix = model.parameters.(type).lmatrix;
    flap_matrix(lmatrix == 1) = angle_rad;
    
    % Assign back to model
    model.geo.flap_vector = flap_matrix;

end
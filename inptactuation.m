function [model] = inptactuation(model)
% INPTACTUATION - Interactively assigns actuation components to a model.
%
% Provides a menu loop to assign control surfaces (delta_e, delta_a, delta_r) or 
% propulsion to a model using a question-based interface. Updates the model struct 
% based on user selections until exit is chosen.
%
% INPUTS: model - Struct with model data to be updated.
% OUTPUT: model - Updated model with assigned actuation components.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

try
    
    while true
       answ=questions(16);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)

           case 0
               
               break;

           case 1
                
               model = assign_controlsurface(model, 'delta_e');
                
           case 2

               model = assign_controlsurface(model, 'delta_a');

           case 3 

               model = assign_controlsurface(model, 'delta_r');

            case 4 

               model = assign_propulsion(model);


           otherwise
           

       end

    end

catch

end

end
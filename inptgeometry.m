function [model]=inptgeometry(model)
% INPTGEOMETRY - Interactively edits or displays a model’s geometry.
%
% Offers a menu loop to edit geometry using Tornado or Redspot methods, or display 
% it, via a question-based interface. Updates the model struct based on user choices 
% until exit is selected.
%
% INPUTS: model - Struct with model data to be modified or displayed.
% OUTPUT: model - Updated model with edited geometry.
%
% Author: Ludwig Horvath
% Date: 3/17/2025

try

loop=1;

    while loop==1
        q=questions(3);
        
        if isempty(q)
           q=-1; %will go to otherwise section
        end
        
        
            switch q
                  
                case(1)   
                   [model]=edit_tornado_geo(model);
                
                case(2)
                   [model]=edit_redspot_geo(model);
                   
                case(3) 
                   display_geometry(model, false);
                
                case 0
                  loop=0;
                
                otherwise
	                terror(9);	
            end
    end
catch

end
end


    
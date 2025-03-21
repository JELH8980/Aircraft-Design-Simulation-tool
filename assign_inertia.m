function [model] = assign_inertia(model)


try
    while true
    q=questions(10);
    
        if isempty(q)
           q=-1; %will go to otherwise section
        end
        
        switch (q)
    
            case 1      % Input

                model = calculate_inertia_1(model);

            case 2      % Calculate

                model = calculate_inertia_2(model);

            case 3

                model = calculate_inertia_3(model);
                
    
            case 0  % Return to main menu
                break;
    
        end
    
    end

catch
    % Handle errors here
    emergency_save(model);

end


end
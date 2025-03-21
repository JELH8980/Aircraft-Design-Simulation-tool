function [model] = simulate(model)

try
    
    while true
       answ=questions(20);       %question string generator function
       if isempty(answ)
          answ=-1;
       end
       
       switch (answ)

           case 0
               
               break;

           case 1

               [model] = simulate(model);

           case 2
               break;

               %[model] = display_results(model);

           otherwise

       end

    end

catch
    disp('Error in "inptsimulation.m"')
end



end
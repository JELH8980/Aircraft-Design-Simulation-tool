function [output]=questions(no)
% QUESTIONS - Function to display menus and capture user input for different model setup options.
%
% This function provides a series of menus for the user to interact with, specifically related to a model setup process.
% The user is prompted with a choice of options depending on the input argument `no`. The menus cover a range of operations,
% including geometry setup, inertia setup, dynamics setup, and material assignment.
%
% INPUTS:
%   no           - A scalar integer representing the menu number to display. The corresponding menu options will be displayed 
%                  based on this number. The menu structure is hierarchical, with different levels of options depending on the
%                  value of `no`.
%
% OUTPUTS:
%   output       - The user's choice from the displayed menu options. This is returned as the function output.
%
% FUNCTIONALITY:
% - Based on the value of `no`, a menu is displayed, and the user is prompted to enter a choice.
% - The function handles various stages of the model setup process, including geometry, inertia, and dynamics setups.
% - Each menu displays a set of options relevant to the current stage, such as defining new geometry, assigning materials,
%   or editing the wing structure.
% - The user’s input is captured and returned for further use in the model setup process.
% 
% EXAMPLES:
% - Calling `questions(1)` will display the main menu with options to create or load a model, or exit the program.
% - Calling `questions(12)` will display the geometry editor menu, where the user can add or remove wings and partitions.
%
% NOTES:
% - The function uses `input` to capture user choices. The user should enter a valid choice based on the options shown.
% - The function does not handle incorrect user inputs directly; it assumes the user will enter a valid integer from the list.
%
% Author: Ludwig Horvath
% Date: 2/11/2025


clc;

switch no
    case(1)
        disp('______________________________________________________')
        disp('  Red Spot  Version 1.0 Alpha version                 ')
        disp('  build 2025 3 2                                      ')
        disp('  Main Menu                                           ')
        disp('______________________________________________________')
        disp(' ')
        disp(' Input operations. ');
        disp('	[1]. Create a new model')
        disp('	[2]. Load an existing model')
        disp('	[3]. Model setup menu')
        disp('	[0]. Exit Red Spot')
        disp(' ')
        
        output=input('	Please enter choice from above: ');


    case(2)
    
        disp('______________________________________________________')
        disp('  Red Spot  Version 1.0 Alpha version                 ')
        disp('  build 2025 3 2                                      ')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('______________________________________________________')
        disp(' ')
        disp(' Input operations. ');
        disp('	[1]. Geometry setup menu')
        disp('	[2]. Inertia setup menu')
        disp('	[3]. Actuation setup menu')
        disp('	[4]. Aerodynamics setup menu')
        disp('	[5]. Version control menu')
        disp('	[6]. Rename model')
        disp('	[0]. Back / up menu')
        disp(' ')
        
        output=input('	Please enter choice from above: ');

    case(3)
        
        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Geometry setup menu                 ')
        disp('                                                      ')
        disp('______________________________________________________')
        
        disp(' ')
        disp('	[1]. Edit Tornado geometry')
        disp('	[2]. Edit Redspot geometry')
        disp('    [3]. Display geometry     ')
        disp(' ')

        disp('	[0]. Back / up menu')
        disp(' ')

        output=input(' 	Please enter choice from above: ');

    case(4)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                                                      ')
        disp('______________________________________________________')
        
        disp(' ')
        disp('	[1]. Assign massive patches ')   
        disp('	[2]. Assign patch thickness ')
        disp('	[3]. Assign material ')
        disp('	[4]. Assign elements ')
        disp('	[5]. Assign inertia ')
        disp('	[6]. Display massive model ')
        
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')

        output=input(' 	Please enter choice from above: ');

    case(5)
        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign massive patches        ')
        disp('______________________________________________________')
        
        disp(' ')
        disp('	[1]. Continue assignment')   
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')

        output=input(' 	Please enter choice from above: ');

    
    case(6)
        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign patch region           ')
        disp('______________________________________________________')

        disp(' ')

        output=input('  Please enter name of region: ', 's');

    case(7)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign patch thickness        ')
        disp('______________________________________________________')

        disp(' ')
        disp('	[1]. Continue assignment')   
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')


    output=input('    Please enter choice from above: ');

    case(8)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign material               ')
        disp('______________________________________________________')

        disp(' ')
        disp('	[1]. Assign material') 
        disp('	[2]. Assign material composition') 
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')
       
        output=input(' Please enter choice from above: ');
    
    case(9)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign elements               ')
        disp('______________________________________________________')

        disp(' ')
        disp('	[1]. Continue assignment') 
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')
       
        output=input('    Please enter choice from above: ');
    
    case(10)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Inertia setup menu                  ')
        disp('                  |---->Assign inertia                ')
        disp('______________________________________________________')

        disp(' ')
        disp('	[1]. Return [m, Ixx, Iyy, Izz](m, Ixx, Iyy, Izz)')
        disp('	[2]. Return [m, Ixx, Iyy, Izz](m, model)')
        disp('	[3]. Return [m, Ixx, Iyy, Izz](material, composition, model)')
        disp(' ')
        disp('	[0]. Back / up menu')
        disp(' ')

    output=input('    Please enter choice from above: ');

    case(12)
           disp(' ')
	       disp('______________________________________________________')
           disp('                                                      ')
           disp('  Main Menu                                           ')
           disp('      |---->Model setup menu                          ')
           disp('           |---->Geometry setup menu                  ')
           disp('                |---->Geometry editor menu            ')
           disp('                                                      ')
	       disp('______________________________________________________')
           disp(' ')
           disp('[1] Add Wing')
           disp('[2] Remove Wing')
           disp(' ')
           disp('[3] Add partition to a wing')
           disp('[4] Remove partition from a wing')
           disp(' ')
           disp('[5] View wing data')
           disp(' ')
           disp('[6] Edit wing/partition data')   
           disp(' ')
           disp('[0] Back / up menu')
           disp(' ')
           output=input(' Please enter choice from above: ');


   case(13)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Batch setup menu                ')
        disp('______________________________________________________')

        disp(' ')
        disp('	[1]. Define/Edit batch') 
        disp('    [2]. Load batch')
        disp('	[3]. Display batch')
        disp('    [4]. Solve batch')
        disp('	[0]. Back / up menu')
        disp(' ')
        output=input(' Please enter choice from above: ');


        case(14)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Geometry setup menu                 ')
        disp('                  |---->Edit Redspot geometry         ')
        disp('______________________________________________________')

        disp(' ')
        disp('    [1]. Edit Origin   (G)       ')
        disp('    [2]. Edit Thrust   (T, aT)   ')
        disp(' ')
        disp('    [0]. Back / up menu           ')
        disp(' ')
        output=input(' Please enter choice from above: ');

        case(15)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Actuation setup menu                ')
        disp('                  |---->Assign/Edit control surface      ')
        disp('______________________________________________________')

        disp(' ')
        disp('    [1]. Assign/Edit logic flap matrix   ')
        disp('    [2]. Assign/Edit deflection boounds  ')
        disp('    [3]. Assign/Edit time constant       ')
        disp('    [4]. Simulate step response          ')
        disp('    [5]. Inspect actuator model           ')
        disp(' ')
        disp('    [0]. Back / up menu           ')
        disp(' ')
        output=input('    Please enter choice from above: ');


        case(16)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Actuation setup menu                ')
        disp('______________________________________________________')

        disp(' ')
        disp('    [1]. Assign elevator    (delta_e)                    ')
        disp('    [2]. Assign aileron     (delta_a)                    ')
        disp('    [3]. Assign rudder      (delta_r)                    ')
        disp('    [4]. Assign propulsion                               ')
        disp(' ')
        disp('    [0]. Back / up menu           ')
        disp(' ')
        output=input('    Please enter choice from above: ');


        case(17)

        disp(' ')
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Model setup menu                          ')
        disp('            |---->Actuation setup menu                ')
        disp('______________________________________________________')

        disp(' ')
        disp('    [1]. Assign propulsion model                         ')
        disp('    [2]. Assign maximum thrust (Tmax)                    ')
        disp('    [3]. Assign/Edit time constant                       ')
        disp('    [4]. Simulate step response                          ')
        disp('    [5]. Inspect actuator model                           ')
        disp(' ')
        disp('    [0]. Back / up menu           ')
        disp(' ')
        output=input('    Please enter choice from above: ');

        case(18)
        disp('______________________________________________________')
        disp('  Main Menu                                           ')
        disp('      |---->Load an existing model                    ')
        disp('______________________________________________________')
        disp(' ')
        



    otherwise

end
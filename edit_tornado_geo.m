function [model]=edit_tornado_geo(model)


if model.geo.nwing==0
   disp(' ')
   disp(' No geometry loaded')
   disp(' ')
else
    
try
    loop=1;  

   while loop==1 

        emergency_save(model);



         disp(' ')
         disp(strcat('Number of wings are :  ',num2str(model.geo.nwing)))
         disp(strcat('Number of partition per wing are :  ',num2str(model.geo.nelem)))
        
         choice1=questions(12);
         
         if isempty(choice1)
                choice1=10;
                terror(9)
         end
            disp(' ')
      
        switch choice1
      
                case 1
	                disp('New wing will be added as higest wing number. ')
	                disp(' ')
	                disp('Dont forget to set the new wings properties')
		                
                model.geo.nwing=model.geo.nwing+1;  
		                model.geo.nelem(model.geo.nwing)=1;            
                model.geo.startx(model.geo.nwing)=0;
	                model.geo.starty(model.geo.nwing)=0;
	                model.geo.startz(model.geo.nwing)=0;
	                model.geo.symetric(model.geo.nwing)=0;              
                model.geo.c(model.geo.nwing)=1;
                model.geo.foil(model.geo.nwing,:,1)={'0'};  %inboard profile
	                model.geo.foil(model.geo.nwing,:,2)={'0'};  %outboard profile              
                model.geo.nx(model.geo.nwing,:)=1;
                model.geo.TW(model.geo.nwing,:,1)=0;        %inboard twist
	                model.geo.TW(model.geo.nwing,:,2)=0;        %outboard twist
                model.geo.dihed(model.geo.nwing,:)=0;
                model.geo.ny(model.geo.nwing,:)=1;
                model.geo.b(model.geo.nwing,:)=1;
                model.geo.T(model.geo.nwing,:)=1;
                model.geo.SW(model.geo.nwing,:)=0;                
                model.geo.meshtype(model.geo.nwing,:)=1;    
	                model.geo.fc(model.geo.nwing,:)=0;
                model.geo.fnx(model.geo.nwing,:)=0;
                model.geo.fsym(model.geo.nwing,:)=0;
                model.geo.flap_vector(model.geo.nwing,:)=0;
	                model.geo.flapped(model.geo.nwing,:)=0;
                %T136 additions
                model.geo.allmove(model.geo.nwing)=0;
                model.geo.allmove_origin(model.geo.nwing,:)=[0 0 0];
                model.geo.allmove_axis(model.geo.nwing,:)=[0 0 0];
                model.geo.allmove_symetric(model.geo.nwing,:)=1;
	                 
                
                case 2
	                disp(' ')
	                choice2=input('Remove wing no: ');
	                    model.geo.nwing=model.geo.nwing-1;
                model.geo.nelem(choice2)=[];
	                    model.geo.SW(choice2,:)=[];
	                    model.geo.T(choice2,:)=[];
                model.geo.flap_vector(choice2,:)=[];
		                model.geo.TW(choice2,:,:)=[];
		                model.geo.b(choice2,:)=[];
    	                model.geo.c(choice2)=[];
		                model.geo.dihed(choice2,:)=[];
		                model.geo.fc(choice2,:)=[];
    	                model.geo.flapped(choice2,:)=[];
	                    model.geo.fnx(choice2,:)=[];
	                    model.geo.foil(choice2,:,:)=[];
	                model.geo.fsym(choice2,:)=[];
	                
	                model.geo.nx(choice2,:)=[];
	                model.geo.ny(choice2,:)=[];
	                model.geo.startx(choice2)=[];
	                model.geo.starty(choice2)=[];
	                model.geo.startz(choice2)=[];
	                model.geo.symetric(choice2)=[];
                model.geo.meshtype(choice2,:)=[];
                
                %T136 additions
                model.geo.allmove(choice2)=[];
                model.geo.allmove_origin(choice2,:)=[];
                model.geo.allmove_axis(choice2,:)=[];
                model.geo.allmove_symetric(model.geo.nwing,:)=[];

	                case 3
	                disp(' ')
	                disp('Partition will be added at wingtip ')
	                choice3=input('Add partition to wing no: ');
      
         
	                model.geo.nelem(choice3)=model.geo.nelem(choice3)+1;
	                i=model.geo.nelem(choice3);
	                model.geo.b(choice3,i)=0;
                model.geo.SW(choice3,i)=0;
	                model.geo.T(choice3,i)=0;
	                model.geo.TW(choice3,i,1)=0;
	                model.geo.TW(choice3,i,2)=0;
	                model.geo.dihed(choice3,i)=0;
	                model.geo.fc(choice3,i)=0;
	                model.geo.flapped(choice3,i)=0;
                model.geo.flap_vector(choice3,i)=0;
                model.geo.fnx(choice3,i)=0;
	                model.geo.foil(2*choice3-1,i,:)={'0'};
	                model.geo.foil(2*choice3,i,:)={'0'};
	                model.geo.fsym(choice3,i)=0;
	                model.geo.nx(choice3,i)=0;
	                model.geo.ny(choice3,i)=0;
                model.geo.meshtype(choice3,i)=1;
         
                case 4
	                disp(' ')
      
	                choice3=input('Remove partition from wing no: ');
	                choice4=input('Remove partition no: ');
	                no_of_elements=model.geo.nelem(choice3);
      
	                for step=choice4:no_of_elements %step from selected element outwards.   
		                if step<model.geo.nelem(choice3)
			                model.geo.SW(choice3,step)=model.geo.SW(choice3,step+1);	%moving properties one step inwards
			                model.geo.T(choice3,step)=model.geo.T(choice3,step+1);
                        model.geo.flap_vector(choice3,step)=model.geo.flap_vector(choice3,step+1);
			                model.geo.TW(choice3,step,1)=model.geo.TW(choice3,step+1,1);
			                model.geo.TW(choice3,step,2)=model.geo.TW(choice3,step+1,2);
			                model.geo.b(choice3,step)=model.geo.b(choice3,step+1);
		                    model.geo.dihed(choice3,step)=model.geo.dihed(choice3,step+1);
			                model.geo.fc(choice3,step)=model.geo.fc(choice3,step+1);
	                    model.geo.flapped(choice3,step)=model.geo.flapped(choice3,step+1);
			                model.geo.fnx(choice3,step)=model.geo.fnx(choice3,step+1);
                        model.geo.foil(choice3,step,1)=model.geo.foil(choice3,step+1,1);
			                model.geo.foil(choice3,step,2)=model.geo.foil(choice3,step+1,2);
                        model.geo.fsym(choice3,step)=model.geo.fsym(choice3,step+1);
			                model.geo.nx(choice3,step)=model.geo.nx(choice3,step+1);
			                model.geo.ny(choice3,step)=model.geo.ny(choice3,step+1);
                        
                        model.geo.meshtype(choice3,step)=model.geo.meshtype(choice3,step+1);
		                end
	                end
       	                model.geo.SW(choice3,no_of_elements)=0;%Deleting properties on most outboard element		
				                model.geo.T(choice3,no_of_elements)=0;
                        model.geo.flap_vector(choice3,no_of_elements)=0;         
			                    model.geo.TW(choice3,no_of_elements,1)=0;
				                model.geo.TW(choice3,no_of_elements,2)=0;  
				                model.geo.b(choice3,no_of_elements)=0	;
			                    model.geo.dihed(choice3,no_of_elements)=0;
				                model.geo.fc(choice3,no_of_elements)=0;
			                model.geo.flapped(choice3,no_of_elements)=0;
				                model.geo.fnx(choice3,no_of_elements)=0;
				                model.geo.foil(choice3,no_of_elements,1)={'0'};
				                model.geo.foil(choice3,no_of_elements,2)={'0'};
				                model.geo.fsym(choice3,no_of_elements)=0;
				                model.geo.nx(choice3,no_of_elements)=0;
   		                    model.geo.ny(choice3,no_of_elements)=0;
                        model.geo.meshtype(choice3,no_of_elements)=0;
                        
                
                
                
      
	                if sum(model.geo.b(:,no_of_elements))==0
	                    model.geo.SW(:,no_of_elements)=[];%Deleting zeroz column[]		
	 		                model.geo.T(:,no_of_elements)=[];
                    model.geo.flap_vector(:,no_of_elements)=[];    
	 		                model.geo.TW(:,no_of_elements,:)=[];
	 		                model.geo.b(:,no_of_elements)=[];
 		                model.geo.dihed(:,no_of_elements)=[];
			                model.geo.fc(:,no_of_elements)=[];
 			                model.geo.flapped(:,no_of_elements)=[];
	 		                model.geo.fnx(:,no_of_elements)=[];
	 	                model.geo.foil(:,no_of_elements,:)=[];
			                model.geo.fsym(:,no_of_elements)=[];
			                model.geo.nx(:,no_of_elements)=[];
   	                        model.geo.ny(:,no_of_elements)=[];
                    model.geo.meshtype(:,no_of_elements)=[];
	                end
      
         
               model.geo.nelem(choice3)=model.geo.nelem(choice3)-1;
               
               
            case 5   
	                disp(' ')
	                wedit=input('View wing no: ');   
	                disp(' ')
		                disp('______________________________________________________')
	                disp('                                                    ')
	                disp('   WING DATA                                            ')
                disp(strcat('     Wing number:  ',num2str(wedit),('                  ')))
		                disp('______________________________________________________')
   
	                disp(strcat('   Number of Partitions	    	:  ',num2str(model.geo.nelem(wedit))))
	                disp(' ')   
	                disp('Global entries ')
                disp(' ')
                disp(strcat('   Model name       	         :  ',num2str(model.name)))
                disp(' ');
                disp(strcat('   Reference point position     :  ',num2str(model.geo.ref_point)))
               
                
                
                disp(' ')
                
                disp('Wing specific entries ')   

                disp(' ')
	                disp(strcat('   Wing Symmetric		    	:  ',num2str(model.geo.symetric(wedit))))
                    disp(strcat('   Apex coordinates	        	:  ',num2str([model.geo.startx(wedit) model.geo.starty(wedit) model.geo.startz(wedit) ])))
                    disp(strcat('   Base chord 			        :  ',num2str(model.geo.c(wedit))))
                disp(strcat('   All moving surface	        :  ',num2str(model.geo.allmove(wedit))))
                disp(strcat('       Origin	                :  ',num2str(model.geo.allmove_origin(wedit,:))))
                disp(strcat('       Axis         	        :  ',num2str(model.geo.allmove_axis(wedit,:))))
                disp(strcat('       Deflection symmetry     :  ',num2str(model.geo.allmove_symetric(wedit))))
                    disp(' ')    
                    disp('Partition specific entries')   
                    disp(' ')   
  
	                disp(strcat('   Partitions half-span      	:  ',num2str(model.geo.b(wedit,:))))
	                disp(strcat('   Partitions sweep 		    :  ',num2str(model.geo.SW(wedit,:))))
                disp(strcat('   Partitions Dihedral          :  ',num2str(model.geo.dihed(wedit,:))))

	                disp(strcat('   Partitions taper 		    :  ',num2str(model.geo.T(wedit,:))))
                
	                disp('   Partitions inner airfoil 	:  '),disp(model.geo.foil(wedit,:,1))
	                disp('   Partitions inner airfoil 	:  '),disp(model.geo.foil(wedit,:,2))
            
                disp(strcat('   Partition inner twists      	:  ',num2str(model.geo.TW(wedit,:,1))))
	                disp(strcat('   Partition outer twists 	    :  ',num2str(model.geo.TW(wedit,:,2))))
	                disp(' ')
	                disp(strcat('   Partition flapped  		    :  ',num2str(model.geo.flapped(wedit,:))))
	                disp(strcat('   Flap chords (Parts)		    :  ',num2str(model.geo.fc(wedit,:))))
	                disp(strcat('   Flaps deflect symmetric	    :  ',num2str(model.geo.fsym(wedit,:))))
	                disp(' ') 
	                disp(strcat('   No. Chord-wise panels  	    :  ',num2str(model.geo.nx(wedit,:))))
	                disp(strcat('   No. Span-wise panels         :  ',num2str(model.geo.ny(wedit,:))))
	                disp(strcat('   No. Flap-chord panels  	    :  ',num2str(model.geo.fnx(wedit,:))))
                disp(strcat('   Panel distribution            :  ',num2str(model.geo.meshtype(wedit,:))))
                disp(strcat('   Flap setting           	    :  ',num2str(model.geo.flap_vector(wedit,:))))
	                disp(' ') 
	                disp(' Paused, press space to continue. ') 
               pause
               
            case 6
	                disp(' ')
	                wedit=input('Edit wing no: ');  
	                eedit=input('Edit partition no: '); 
	                loop2=1;
            while loop2==1
		                disp(' ')
		                disp('_________________________________________________')
		                disp('                                                 ')
		                disp('   Wing Data                                     ')
		                disp('_________________________________________________')
        
		                disp(strcat('    Number of Partitions       :',num2str(model.geo.nelem(wedit))))
		                disp(' ')
                    disp('Global entries ')
                    disp(' ')
                    disp(strcat(' [1] Reference point position       :',num2str(model.geo.ref_point)))      
                    disp(' ')
		                disp('Wing specific entries ')
		                disp(' ')
		                    disp(strcat(' [2] Wing Symmetric                 :',num2str(model.geo.symetric(wedit))))
		                disp(strcat(' [3] Apex coordinates               :',num2str([model.geo.startx(wedit) model.geo.starty(wedit) model.geo.startz(wedit) ])))
		                disp(strcat(' [4] Base chord                     :',num2str(model.geo.c(wedit))))
                    disp(strcat(' [5] All moving surface	            :  ',num2str(model.geo.allmove(wedit))))
                    disp(strcat(' [6]      Origin	                :  ',num2str(model.geo.allmove_origin(wedit,:))))
                    disp(strcat(' [7]     Axis         	            :  ',num2str(model.geo.allmove_axis(wedit,:))))
                    
		                    disp(' ')
		                    disp('Partition specific entries ')
		                    disp(' ')
                    disp(strcat(' [8] Dihedral                      :',num2str(model.geo.dihed(wedit,eedit))))
	                    disp(strcat(' [9] Partition half-span           :',num2str(model.geo.b(wedit,eedit))))
		                disp(strcat(' [10] Partition sweep               :',num2str(model.geo.SW(wedit,eedit))))
		                disp(strcat(' [11] Partition taper               :',num2str(model.geo.T(wedit,eedit))))
		                   
                    disp(' [12] Partition inner airfoil       :'),disp((model.geo.foil(wedit,eedit,1)))
                    disp(' [13] Partition outer airfoil       :'),disp((model.geo.foil(wedit,eedit,2)))

                    disp(strcat(' [14] Partition inner twist         :',num2str(model.geo.TW(wedit,eedit,1))))
		                disp(strcat(' [15] Partition outer twist         :',num2str(model.geo.TW(wedit,eedit,2))))
		                disp(' ')
		                disp(strcat(' [16] Partition flapped             :',num2str(model.geo.flapped(wedit,eedit))))
		                disp(strcat(' [17] Flap chords (Parts)           :',num2str(model.geo.fc(wedit,eedit))))
		                disp(strcat(' [18] Flaps deflect symmetric       :',num2str(model.geo.fsym(wedit,eedit))))
		                disp(' ') 
		                disp(strcat(' [19] No. of chord-wise panels      :',num2str(model.geo.nx(wedit,eedit))))
		                disp(strcat(' [20] No. of span-wise panels       :',num2str(model.geo.ny(wedit,eedit))))
		                disp(strcat(' [21] No. of flap-chord panels      :',num2str(model.geo.fnx(wedit,eedit))))
                    
                    
                    disp(strcat(' [22] Panel Distribution            :',num2str(model.geo.meshtype(wedit,eedit))))
                    

                    disp(' ')
		                disp('[0]	EXIT ') 
		                disp(' ') 	

	                edit2=input('Edit Menu Item [1-22]: ');
                    if isempty(edit2)
			                edit2=18;
			                terror(9)
                    end
         
                    switch edit2


                        case 1
                            model.geo.ref_point=str2num(input('Reference Point [x y z]: ','s'));
  
                        case 2
                            model.geo.symetric(wedit)=input('Wing symmetry bit [1 0]: ');
                            
                            
                        case 3
				                model.geo.startx(wedit)=input('Apex x-coord: ');
				                model.geo.starty(wedit)=input('Apex y-coord: ');
				                model.geo.startz(wedit)=input('Apex z-coord: ');
                            
			                case 4
                            model.geo.c(wedit)=input('Base chord: ');
                            
                        case 5
                            model.geo.allmove(wedit)=input('All Moving Surface bit [1 0]: ');
                            
                        case 6
                            model.geo.allmove_origin(wedit,:)=str2num(input('All Moving Surface origin [x y z]: ','s'));
                            
                        case 7
                            model.geo.allmove_axis(wedit,:)=str2num(input('All Moving Surface axis [x y z]: ','s'));

                        case 8
			                    model.geo.dihed(wedit,eedit)=input('Partition dihedral: ');	
                            
                        case 9
                            model.geo.b(wedit,eedit)=input('Partition Span: ');
                            
			                case 10
				                model.geo.SW(wedit,eedit)=input('Partition sweep [rad]: ');
                            
			                case 11
				                model.geo.T(wedit,eedit)=input('Partition taper : ');
                            
                        case 12
				                data=input('Partition inner foil : ','s');
                            model.geo.foil(wedit,eedit,1)={data};
                            disp('Remember to change adjacent partition to ensure continuity');
                            
                        case 13
                            data=input('Partition outer foil : ','s');
                            model.geo.foil(wedit,eedit,2)={data};
                            disp('Remember to change adjacent partition to ensure continuity');
                            
			                case 14  
				                model.geo.TW(wedit,eedit,1)=input('Partition inner twist [rad]: ');
                            disp('Remember to change adjacent partition to ensure continuity');
                            
 		                    case 15
				                model.geo.TW(wedit,eedit,2)=input('Partition outer twist [rad]: ');
                            disp('Remember to change adjacent partition to ensure continuity');
                            
			                case 16
				                 model.geo.flapped(wedit,eedit)=input('Is Partition flapped : ');
                             
			                case 17
				                 model.geo.fc(wedit,eedit)=input('Flap chord (parts) : ');
                             
			                case 18
				                 model.geo.fsym(wedit,eedit)=input('Symmetric deflection : ');
                             
		                    case 19
				                 model.geo.nx(wedit,eedit)=input('No: Chordwise panels  : ');
                             
                        case 20
                            model.geo.ny(wedit,eedit)=input('No: spanwise panels  : ');
                            
			                case 21
				                model.geo.fnx(wedit,eedit)=input('No: Chordwise panels on flap  : ');
              
                        case 22
                            model.geo.meshtype(wedit,eedit)=input('Meshtype  [1 - 7]: ');
                            
                        case 0
  			                loop2=0; 
			                otherwise
   			                terror(9)   
                    end
            end
                
            case 8
                keyboard

            case 0
		                loop=0;
                case -1
                    %donothing
                    terror(9) 
            case 10
                
            otherwise
                terror(9)     
        end



   end
    


catch 
    
end

end


g = model.geo;

state.AS=1;					%airspeed
state.alpha=0;				%angle of attack
state.betha=0;				%angle of sideslip
state.P=0;					%roll angluar rate	
state.Q=0;					%pitch angular rate
state.R=0;					%Yaw angular rate
state.alphadot=0;           %Angle of attack time derivative
state.bethadot=0;           %Angle of sidesliptime derivative
state.ALT=0;                %Altitude, meters.
state.rho=0;                %Air density, kg/m^3.
state.pgcorr=0;             %Prandtl-Glauert compressibillity correction.

g.nx=double(g.nx>0);
g.ny=double(g.ny>0); 
g.fnx=double(g.fnx>0); 

[mesh, ~] = fLattice_setup(g, state, 1);

model.geo.mesh = mesh;

end



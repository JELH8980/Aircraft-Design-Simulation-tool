function [model] = assign_geometry()

try
        settings=config('startup'); %setting directories
        
        model = struct(); 
         
        model.name = input('    Enter name of model: ', 's');  
      
        disp(' ')
        %ok1=0;
        geo=[];
        geo.version=136;        
        
        disp('_______________________________________________________')
        disp('                                                      ')
        disp('              New Aircraft Definition                 ')
        disp('                -Global Parameters-                   ')
        disp('_______________________________________________________')
        disp(' ');
        disp('    b - Back one question. ')
        disp('    q - Abort input sequence. ')
        disp(' ')
        disp('-------------------------------------------------------')
        
       stepper=2;
        
       geo.ref_point(1)=0;
       geo.ref_point(2)=0;
       geo.ref_point(3)=0;

        while (stepper < 99)
    
           switch stepper 
   

                case 2
                    disp (' ')
                    data=input('Number of Wings: ','s');
                    type=1;
                    if isinptok(data,type)==1
                       geo.nwing=str2double(data);
                       stepper=99;
                    end
                    disp (' ')
                    
           end  %case  
           
           if isinptok(data,type)==1
               stepper=stepper+1;
           elseif isinptok(data,type)==-1
               stepper=stepper-1;
           elseif isinptok(data,type)==-2
               stepper=99;
               return
           end %If     
        end  %while
    
     
	    for s=1:geo.nwing			%Stepping over the wings.
		    %home	
		 	    disp('______________________________________________________')
		    disp('                                                      ')
            disp(strcat('  Data regarding wing number :', num2str(s),'  '))
		    disp('______________________________________________________')
            disp(' ');
         
         stepper3=1;
         while stepper3 <99
             switch stepper3
             
                 case 1     
                 
                    data=input('Number of semispanwise partitions for this wing: ','s');
                    type=1;
                    if isinptok(data,type)==1
                        geo.nelem(s)=str2double(data);     
                    end
                    
                 case 2
            
                   data=input('Is the wing mirrored in the xz-plane [1 0]: ','s');
                   type=1;
                   if isinptok(data,type)==1
                    geo.symetric(s)=str2double(data);
                   end
                   
                 
                   
                 case 3

                    if s == 1
                        geo.startx(s)=0;
                    else
                        disp(' ')
                        data=input('Apex x-coordinate: ','s');
                        type=1;
                        if isinptok(data,type)==1
                            geo.startx(s)=str2double(data);
                        end
                    end
                 case 4
                    
                    if s == 1
                        geo.starty(s)=0;
                    else
                        data=input('Apex y-coordinate: ','s');
                        type=1;
                        if isinptok(data,type)==1
                            geo.starty(s)=str2double(data);
                        end
                    end
             
                 case 5
                    
                    if s == 1
                        geo.startz(s)=0;
                    else
                        data=input('Apex z-coordinate: ','s'); 
                        type=1;
                        if isinptok(data,type)==1
                            geo.startz(s)=str2double(data);
                        end  
                    end
                   
                case 6
                    disp(' ') 
                    data=input('Is this wing an all-moving control surface [1 0]: ','s');
                    type=1;
                    if isinptok(data,type)==1
                        geo.allmove(s)=str2num(data);
                    end
                    if geo.allmove==0
                         stepper=100;
                    end
                  
                 case 7
                     
                        if geo.allmove(s)==1
        				    data=input('Hinge Line Origin [x y z]: ','s');
                            type=1;
                            if isinptok(data,type)==1
                               geo.allmove_origin(s,:)=str2num(data);
                            end
                        else
                            geo.allmove_origin(s,:)=[0 0 0];
                        end
                 
                case 8
                     
                 if geo.allmove(s)==1
        			    data=input('Hinge Line Axis [x y z]: ','s');
                        type=1;
                        if isinptok(data,type)==1
                           geo.allmove_axis(s,:)=str2num(data);
                        end
                 else
                     geo.allmove_axis(s,:)=[0 1 0];
                 end
                 geo.allmove_def(s,:)=0;
                    
             
                 case 9
                    if geo.allmove(s)==1
                        if geo.symetric(s)==1
            			    data=input('Symmetric deflection [0 1]: ','s');
                            type=1;
                            if isinptok(data,type)==1
                               geo.allmove_symetric(s,:)=str2num(data);
                            else
                               geo.allmove_symetric(s,:)=str2num(data);
                            end
                        end
                    else
                     geo.allmove_symetric(s,:)=0;
                    end
             stepper3=99;
             end
    
             
             
             
                       
            if isinptok(data,type)==1
               stepper3=stepper3+1;
            elseif isinptok(data,type)==-1
               stepper3=stepper3-1;
            elseif isinptok(data,type)==-2
               stepper3=99;
               return
            else
            end
         end
        
         
		    t=0;
        for s2=1:geo.nelem(s)
     	    t=t+1;
         
     	    disp(' ')
		 	       disp('_____________________________________________________')
		       disp('                                                    ')
                disp(strcat(' Wing number:  ',num2str(s),('                  ')))
               disp(strcat(' Data entry for partition number:  ',num2str(t),('                  ')))
            
     	    disp(' ')   
        
            stepper=4;
            while stepper<99
                 switch stepper 
          
                  case 4
                     if t==1  %Only define inboard profile on innermost partition  
                        data=input('Root chord: ','s');
                        type=1;
                        if isinptok(data,type)==1
                           geo.c(s)=str2num(data);
                        end
                     end
                  
                  case 5
                     if t==1  %Only define inboard profile on innermost partition      
                        ok=0;  
                        while ok==0
                            try
                                cd(settings.afdir)
      
    
                                disp(' ')
                                disp('____________________________ ')
                                disp(' ')
                                disp(' AVAILABLE AIRFOILS: ')
                                ls
                                disp('____________________________ ')
                                disp(' ')
                                disp('Enter profile filename from the list above (ex CLARKY.DAT) ')
                                disp('OR any NACA four digits series numer (ex: 2412)')
                                disp('0 (zero) for a flat plate. ')
                                disp(' ')
                                
                                
                                data=input('Base chord airfoil: ','s');
                                
                                
                                if isempty(str2num((cell2mat({data}))))==0
                                   geo.foil(s,t,1)={data};%Naca xxxx profile
                                   ok=1;
                                   cd(settings.hdir)
                                else
                                    try
                                        %foo=str2num(data);
                                        %cd(settings.afdir)
                                             load(data)  %Testload to see that the file exists
                                        cd(settings.hdir)
                                        geo.foil(s,t,1)={data};
                                        data=0;
                                        ok=1;
                                    catch    
                                       cd(settings.hdir)
                                        disp(' ')
                                       disp(' + + + ')
                                        disp(' ')
                                        disp(' No such file! ')
                                       disp(' ')
                                    end
                                end  
                                
                            catch
                                cd(settings.hdir)
                                terror(4)
                            end
                        end
                     end 
                  
            	     case 6   
                     if t==1  %Only on first partition
                        data=input('Number of panels chord wise: ','s');
                        type=1;
                        if isinptok(data,type)==1
                           geo.nx(s,1)=str2num(data);
                        end
                        
                     end
    
    
                    case 7
                     if s==1	% First wing base chord is reference for twist 
          				    geo.TW(s,1,1)=0;
                     else
                            if t==1
                                data=input('Base chord twist [deg]: ','s');
                                type=1;
                                if isinptok(data,type)==1
                                    geo.TW(s,1,1)=str2num(data)*pi/180;
                                end                          
                            end
                            
                        
                     end
                     
                  case 8
                     %if t==1;
                     
    		             data=input('Partition dihedral [deg]: ','s');
                         type=1;
                          if isinptok(data,type)==1
                             geo.dihed(s,t)=str2num(data)*pi/180;
                          end
                     %end
                     
                  case 9
     				    data=input('Number of panels semi-span wise: ','s');
                        type=1;
                     if isinptok(data,type)==1
                        geo.ny(s,t)=str2num(data);
                     end
                     
                  case 10
                     data=input('Span of partition: ','s');
                     type=1;
                     if isinptok(data,type)==1
                        geo.b(s,t)=str2num(data);
                     end
                     
                  case 11
                     data=input('Taper ratio: ','s');
                     type=1;
                     if isinptok(data,type)==1
                        geo.T(s,t)=str2num(data);
                     end
                     
                  case 12
                             
                        ok=0;  
                        while ok==0
                            try
                                cd(settings.afdir)
                                disp(' ')
                                disp('____________________________ ')
                                disp(' ')
                                disp(' AVAILABLE AIRFOILS: ')
                                ls
                                disp('____________________________ ')
                                disp(' ')
                                disp('Enter profile filename from the list above (ex CLARKY.DAT) ')
                                disp('OR any NACA four digits series numer (ex: 2412)')
                                disp('0 (zero) for a flat plate. ')
                                disp(' ')
                                data=input('Tip chord airfoil: ','s');
                                
                                                                
                                if isempty(str2num((cell2mat({data}))))==0
                                   geo.foil(s,t,2)={data};%Naca xxxx profile
                                   ok=1;
                                   cd(settings.hdir)
                                else
                                    try
                                        %foo=str2num(data);
                                        %cd(settings.afdir)
                                             load(data)  %Testload to see that the file exists
                                        cd(settings.hdir)
                                        geo.foil(s,t,2)={data};
                                        data=0;
                                        ok=1;
                                    catch    
                                       cd(settings.hdir)
                                        disp(' ')
                                       disp(' + + + ')
                                        disp(' ')
                                        disp(' No such file! ')
                                       disp(' ')
                                    end
                                end
                                
                            catch
                                cd(settings.hdir)
                                terror(4)
                            end
                        end
                      
    
                  case 13
                     data=input('Quarter chord line sweep [deg]: ','s');
                     type=1;
                     if isinptok(data,type)==1
                        geo.SW(s,t)=str2num(data)*(pi/180);
                     end
                     
                  case 14   
					      data=input('Outboard twist [deg]: ','s');
                      type=1;
                      if isinptok(data,type)==1
                         geo.TW(s,t,2)=str2num(data)*(pi/180);
                      end
    
                  case 15
                      disp(' ')
                      disp(' _______________________ ')
                      disp(' Available mesh distribution types:')
                      disp('   [1] Linear')
                      disp('   [2] Spanwise half-cosine')
                      disp('   [3] Spanwise half-cosine, chordwise cosine')
                      disp('   [5] Spanwise cosine')
                      disp('   [6] Chordwise cosine')% (Added 22/08/2008 AT)
                      disp('   [7] 3:rd order centerpacking. (Not for wings)') 
                      disp(' ')
					      data=input('Mesh type: ','s');
                      type=1;
                      if isinptok(data,type)==1
                         geo.meshtype(s,t)=str2num(data);
                      end
                      
                   case 16
                      disp(' _______________________ ')
                      disp(' ')
                      data=input('Is partition flapped [1 0]:','s');
                      type=1;
                      if isinptok(data,type)==1
                         geo.flapped(s,t)=str2num(data);
                         
                         if geo.flapped==0
                             stepper=100;
                             geo.fc(s,t)=0;
                             geo.fnx(s,t)=0;
                             geo.fsym(s,t)=0;
                         end
                      end
                   case 17    
              	    if geo.flapped(s,t)==1
        				    data=input('Flap chord in fraction of local chord (0..1): ','s');
                            type=1;
                        if isinptok(data,type)==1
                           geo.fc(s,t)=str2num(data);
                        end
                        
                        if geo.fc(s,t)> 0.9
                            disp(' + + + Warning + + +  Are you sure? Consider an allmoving surface instead.')
                        end
                        
                    end
                                       
                  case 18
                     if geo.flapped(s,t)==1
                        data=input('Number of chord wise panels on flap: ','s');
                        type=1;
                        if isinptok(data,type)==1
                           geo.fnx(s,t)=str2num(data);
                        end
                        
                     end
    
                     case 19
                     if and(geo.symetric(s),geo.flapped(s,t))
          				    data=input('Do control surfaces deflect symmetrically [1 0]:','s');
                            type=1;
                          if isinptok(data,type)==1
                             geo.fsym(s,t)=str2num(data);
                          end
                          
                     else
           			    geo.fsym(s,t)=0;
                     end
                     
                  end %caseblock
                  
                  if isinptok(data,type)==1
                     stepper=stepper+1;
                  elseif isinptok(data,type)==-1
                     stepper=stepper-1;
                  elseif isinptok(data,type)==-2
                     stepper=99;
                     return
                  else
                  end
            end% while partitionblock
               
               if geo.flapped(s,t)==0
                  stepper=100;
           	    geo.fc(s,t)=0;
                geo.fnx(s,t)=0;
           	    geo.fsym(s,t)=0;
               end
     		    geo.nx(s,t)=geo.nx(s,1);
        	    %dihed(s,t)=dihed(s,1);
                if t~=1
        		    geo.TW(s,t,1)=geo.TW(s,t-1,2);	%continious twist
        		    geo.foil(s,t,1)=geo.foil(s,t-1,2); %continious camber
                end
             
  	    end%partitionblock
    end%wingblock
    
    geo.flap_vector=zeros(size(geo.flapped));
    
    model.geo = geo;

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
  
    [mesh, ref] = fLattice_setup(g, state, 1);

    model.parameters.S  = ref.S_ref;                %reference area;
    model.parameters.c  = ref.C_mac;                %mean aerodynamic choord
    model.parameters.b  = ref.b_ref;                %reference span
    model.parameters.XM = ref.mac_pos(1);
    model.parameters.YM = ref.mac_pos(2);
    model.parameters.ZM = ref.mac_pos(3);

    model.parameters.lB     = NaN;
    model.parameters.hB     = NaN;

    model.parameters.lR     = 0;
    model.parameters.hR     = 0;

    model.parameters.lT     = NaN;
    model.parameters.hT     = NaN;
    model.parameters.aT     = NaN;

    model.geo.mesh          = mesh;


    model.geo.patches = struct();

    model.geo.bodies  = struct();

    model.parameters.m      = NaN;
    model.parameters.Ixx    = NaN;
    model.parameters.Iyy    = NaN;
    model.parameters.Izz    = NaN;
    model.parameters.Ixz    = NaN;

    model = assign_geometryframe(model);

    model = assign_thrustframe(model);

    flap_matrix = model.geo.flapped;
    
    model.parameters.delta_e         = struct();
    model.parameters.delta_e.lmatrix = zeros(size(flap_matrix));
    model.parameters.delta_e.tau     = 0;
    model.parameters.delta_e.min     = NaN;
    model.parameters.delta_e.max     = NaN;
    model.parameters.delta_e.info    = NaN;


    model.parameters.delta_a         = struct();
    model.parameters.delta_a.lmatrix = zeros(size(flap_matrix));
    model.parameters.delta_a.tau     = 0;
    model.parameters.delta_a.min     = NaN;
    model.parameters.delta_a.max     = NaN;
    model.parameters.delta_a.info    = NaN;

    model.parameters.delta_r         = struct();
    model.parameters.delta_r.lmatrix = zeros(size(flap_matrix));
    model.parameters.delta_r.tau     = 0;
    model.parameters.delta_r.min     = NaN;
    model.parameters.delta_r.max     = NaN;
    model.parameters.delta_r.info    = NaN;

    model.parameters.propulsion        = struct();
    model.parameters.propulsion.tau    = 0;
    model.parameters.propulsion.Tmax   = NaN;
    model.parameters.propulsion.info   = NaN;
    model.parameters.propulsion.method = NaN;

catch
    emergency_save(model)

end


end
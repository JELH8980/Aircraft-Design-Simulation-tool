function [results]=coeff_create(results,lattice,state,ref,geo)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficient creator: Essential function for TORNADO						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes aerodynamic coefficients			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautical 
%                               and Vehicle Engineering	
%			copyright 2003											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidary function for TORNADO					
% Called by:	solverloop											
% Calls:			MATLAB standard functions																			
% Loads:																
% Saves: 												
% Input: 			
% Output: forces moments coefficients							
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=0.01;
q=0.5*state.rho*state.AS^2;			    %calculating dynamic pressure
										%for coefficient calculation		
[a b void]=size(results.F);


npan=cumsum(sum((geo.nx+geo.fnx).*geo.ny,2).*(geo.symetric+1)'); %Number of panels per wing

for s=1:a
	normal_force(s)=squeeze(results.F(s,1,:))'*lattice.N(s,:)';                                
end                                
panel_area=tarea(lattice.XYZ);
stat_press=normal_force./panel_area;	%Delta pressure, top/bottom
results.cp=((stat_press)./(q))';


CX=results.FORCE(:,:,1)/(q*ref.S_ref);
CY=results.FORCE(:,:,2)/(q*ref.S_ref);
CZ=results.FORCE(:,:,3)/(q*ref.S_ref);

B2WTransform=[cos(state.betha)*cos(state.alpha),        -sin(state.betha),          cos(state.betha)*sin(state.alpha) ;...
              cos(state.alpha)*sin(state.betha),         cos(state.betha),          sin(state.betha)*sin(state.alpha) ;...
                              -sin(state.alpha),                        0,                           cos(state.alpha)];
for i=1:b                          
    lemma(i,:)=B2WTransform*squeeze(results.FORCE(:,i,:));
end

D=lemma(:,1)';
C=lemma(:,2)';
L=lemma(:,3)';

CL=L/(q*ref.S_ref);
CD=D/(q*ref.S_ref);
CC=C/(q*ref.S_ref);

Cl=results.MOMENTS(1,:,1)/(q*ref.S_ref*ref.b_ref);
Cm=results.MOMENTS(1,:,2)/(q*ref.S_ref*ref.C_mac);
Cn=results.MOMENTS(1,:,3)/(q*ref.S_ref*ref.b_ref);



%% ------------ CL per wing computation
npan=cumsum(sum((geo.nx+geo.fnx).*geo.ny,2).*(geo.symetric+1)'); %Number of panels per wing

index1=1;

for i=1:geo.nwing
    index2=npan(i);
    
    lemma2=B2WTransform*(sum(squeeze(results.F(index1:index2,1,:))))';
    
    results.CLwing(i)=lemma2(3)/(q*ref.S_ref);
    results.CDwing(i)=lemma2(1)/(q*ref.S_ref);
    results.CYwing(i)=lemma2(2)/(q*ref.S_ref);
    
    index1=npan(i)+1;
end
%% ----------

%%Setting output
results.L=L(1);
results.D=D(1);
results.C=C(1);

results.CX=CX(:,1);
results.CY=CY(:,1);
results.CZ=CZ(:,1);
results.CL=CL(:,1);
results.CD=CD(:,1);
results.CC=CC(:,1);
results.Cl=Cl(:,1);
results.Cm=Cm(:,1);
results.Cn=Cn(:,1);

results.F=squeeze(results.F(:,1,:));
results.M=squeeze(results.M(:,1,:));

results.FORCE=squeeze(results.FORCE(:,1,:));
results.MOMENTS=squeeze(results.MOMENTS(:,1,:));


delta=config('delta');
fac1=ref.b_ref /(2*state.AS);
fac2=ref.C_mac /(2*state.AS);


%%Differentiating
dCX=(CX-CX(:,1))./delta;
dCY=(CY-CY(:,1))./delta;
dCZ=(CZ-CZ(:,1))./delta;

dCL=(CL-CL(:,1))./delta;
dCD=(CD-CD(:,1))./delta;
dCC=(CC-CC(:,1))./delta;

dCl=(Cl-Cl(:,1))./delta;
dCm=(Cm-Cm(:,1))./delta;
dCn=(Cn-Cn(:,1))./delta;




results.CL_a=dCL(2);
results.CD_a=dCD(2);
results.CC_a=dCC(2);
results.CX_a=dCX(2);
results.CY_a=dCY(2);
results.CZ_a=dCZ(2);
results.Cl_a=dCl(2);
results.Cm_a=dCm(2);
results.Cn_a=dCn(2);

results.CL_b=dCL(3);
results.CD_b=dCD(3);
results.CC_b=dCC(3);
results.CX_b=dCX(3);
results.CY_b=dCY(3);
results.CZ_b=dCZ(3);
results.Cl_b=dCl(3);
results.Cm_b=dCm(3);
results.Cn_b=dCn(3);

results.CL_P=dCL(4)/fac1;
results.CD_P=dCD(4)/fac1;
results.CC_P=dCC(4)/fac1;
results.CX_P=dCX(4)/fac1;
results.CY_P=dCY(4)/fac1;
results.CZ_P=dCZ(4)/fac1;
results.Cl_P=dCl(4)/fac1;
results.Cm_P=dCm(4)/fac1;
results.Cn_P=dCn(4)/fac1;

results.CL_Q=dCL(5)/fac2;
results.CD_Q=dCD(5)/fac2;
results.CC_Q=dCC(5)/fac2;
results.CX_Q=dCX(5)/fac2;
results.CY_Q=dCY(5)/fac2;
results.CZ_Q=dCZ(5)/fac2;
results.Cl_Q=dCl(5)/fac2;
results.Cm_Q=dCm(5)/fac2;
results.Cn_Q=dCn(5)/fac2;

results.CL_R=dCL(6)/fac1;
results.CD_R=dCD(6)/fac1;
results.CC_R=dCC(6)/fac1;
results.CX_R=dCX(6)/fac1;
results.CY_R=dCY(6)/fac1;
results.CZ_R=dCZ(6)/fac1;
results.Cl_R=dCl(6)/fac1;
results.Cm_R=dCm(6)/fac1;
results.Cn_R=dCn(6)/fac1;

try
    results.CL_d=dCL(7:end);
    results.CD_d=dCD(7:end);
    results.CC_d=dCC(7:end);
    results.CX_d=dCX(7:end);
    results.CY_d=dCY(7:end);
    results.CZ_d=dCZ(7:end);
    results.Cl_d=dCl(7:end);
    results.Cm_d=dCm(7:end);
    results.Cn_d=dCn(7:end);
end


   
[results]=fStripforce(geo,results,lattice,state,ref);




results.rundate=char(datetime('now'));



%[lemma]=fStripforce(geo,results,lattice,state,ref,vCfraction)

end%function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [results]=fStripforce(geo,results,lattice,state,ref)
%This lemma function computes the aerodynamic force on each strip.

dynp=0.5*state.rho*state.AS^2;
S=ref.S_ref;



%%Regenerating the lattice with whe wings at the centerline in order to
%%have the wings start at [0 0 0]

geo.startx=geo.startx*0;            
geo.starty=geo.starty*0;
geo.startz=geo.startz*0;
geo.dihed=geo.dihed.*0;             %No dihedral to avoid sign==0 at vertical tail



[lattice,ref]=fLattice_setup(geo,state,1);


try
    vCfraction=geo.vCfraction;      %catch cases with no spar position
catch
    vCfraction=0.25*ones(size(geo.c));   %defaulting to put spar shearcenter at c/4 
end



B2WTransform=[cos(state.betha)*cos(state.alpha),        -sin(state.betha),          cos(state.betha)*sin(state.alpha) ;...
                  cos(state.alpha)*sin(state.betha),         cos(state.betha),          sin(state.betha)*sin(state.alpha) ;...
                              -sin(state.alpha),                        0,                           cos(state.alpha)];


F=results.F;                            %Reassigning to save some space
%% Vortex points
[s1 s2 s3]=size(lattice.VORTEX);
if s2==8
    pV1=squeeze(lattice.VORTEX(:,4,:));
    pV2=squeeze(lattice.VORTEX(:,5,:));
elseif s2==4
    pV1=squeeze(lattice.VORTEX(:,2,:));
    pV2=squeeze(lattice.VORTEX(:,3,:));
end
pV=(pV1+pV2)/2;
%%    

[ai bi]=size(geo.nx);           %number of wings and panels
cnx=geo.nx+geo.fnx;             %corrected number of xpanels

cnx2=[];                        %mucking about with index variables
for i=1:ai
    for j=1:bi
        if cnx(i,j)>0
             cnx2=[cnx2 cnx(i,j)];
        end
    end
end



for i=1:geo.nwing;
    cny(i,:)=geo.ny(i,:).*(geo.symetric(i)+1); %corrected number of ypanels
end

stripsperwing=sum(cny,2);

%% Compute force action point and strip moment axis
m=0;
index2=0;
q=0;


for i=1:ai          %loop per wing
   
    for j=1:bi
        if cny(i,j)>=1
            q=q+1;%loop per partition
        end
        for k=1:cny(i,j)  
            %per strip loop
            index1=index2+1;
            index2=index2+cnx2(q);
            m=m+1;
                      
            %% Compute force action point and  strip moment axis
            cornerp=squeeze([lattice.XYZ(index1,1,:);
                              lattice.XYZ(index1,2,:);
                              lattice.XYZ(index2,3,:);
                              lattice.XYZ(index2,4,:)]);
            
           localC1=[(cornerp(1,:)+cornerp(2,:))/2];
           localC2=[(cornerp(3,:)+cornerp(4,:))/2];
           Mpoint=(1-vCfraction(i))*localC1+(vCfraction(i))*localC2;
           
           if sign(Mpoint(2))==0
               SGN=1;
           else
               SGN=sign(Mpoint(2));
           end
           
           yprimestation(m)=SGN*sqrt(Mpoint(2)^2+Mpoint(3)^2);
           
           %Local chord
           lemma1=localC1-localC2;
           lc=sqrt(sum(lemma1.^2));
           
           %local span
           lemma1=(-cornerp(1,:)+cornerp(2,:));
           lemma2=lemma1.*[0 1 1];%Disregarding x component
           ls(m)=sqrt(sum(lemma2.^2));
           
           %Strip Area
           la=ls(m)*lc;
           
            %%
            %Forces
            F0(m)=sum(sqrt(F(index1:index2,2).^2+F(index1:index2,3).^2)); %Only Z and Y component

            h(:,1)=Mpoint(1)-pV(index1:index2,1);
            h(:,2)=Mpoint(2)-pV(index1:index2,2);
            h(:,3)=Mpoint(3)-pV(index1:index2,3);             

            F3(m,:)=sum(F(index1:index2,:),1);
            M3(m,:)=sum(cross(F(index1:index2,:),h),1);
            
            F_w(m,:)=(B2WTransform*F3(m,:)')';              %strip forces in wind system
            Cx_w(m,:)=2*(F_w(m,:))./(state.rho*state.AS^2*la);   %strip coefficients in wind system.
      
            clear h
            %% Coefficients
            CZprime(m)=F0(m)/(dynp*la);
      
        end        
    end
    
    

return
    


end %function stripforce    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [panel_area]=tarea(XYZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tarea: Subsidary function for TORNADO					   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the area of each panel								
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautics	%
%				Copyright 2000											
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidaty function for TORNADO					
% Called by:	coeff_create
% 
% Calls:			MATLAB 5.2 std fcns								
% Loads:	none
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b c]=size(XYZ);
for i=1:a
   p1=[XYZ(i,1,1) XYZ(i,1,2) XYZ(i,1,3)];	%sets up the vectors 
   p2=[XYZ(i,2,1) XYZ(i,2,2) XYZ(i,2,3)];	%to the corners of the		
   p3=[XYZ(i,3,1) XYZ(i,3,2) XYZ(i,3,3)];	%panel.
   p4=[XYZ(i,4,1) XYZ(i,4,2) XYZ(i,4,3)];
   
   a=p2-p1;	%sets up the edge vectors
   b=p4-p1;
   c=p2-p3;
   d=p4-p3;
   
   ar1=norm(cross(b,a))/2;	%claculates the ctoss product of
   ar2=norm(cross(c,d))/2;	%two diagonal corners
   
 	panel_area(i)=ar1+ar2;	%Sums up the product to make the
end						    %Area
end% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[lc]=fLocal_chord2(geo,lattice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Geometry function 						 	%		 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Computes the Local chord at each collocation 
%  point row.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Tomas Melin, KTH, Department of% 
%	Aeronautics, copyright 2002				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Context: Auxillary function for TORNADO%
%	Called by: TORNADO spanload            %
%	Calls:	None									%
%	Loads:	None									%
%	Generates:	Local chord vector lc, same 
%  order as colloc, N, and the others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[indx1 indx2]=size(geo.b);

for s=1:indx1;	   		%Looping over wings
	CHORDS(s,1)=geo.c(s);		%calculating chords of first element
end

for s=1:indx1				%Looping over wings
	for t=1:indx2			%Looping over partitions
	%Chord loop, generating chords for wing partitions
            CHORDS(s,t+1)=CHORDS(s,t)*geo.T(s,t);	%calculating
      												%element root-chord
   end
end




lc=[];	%Local chord vector.


panelchords1=sqrt(sum((lattice.XYZ(:,1,:)-lattice.XYZ(:,4,:)).^2,3)); %inboard
panelchords2=sqrt(sum((lattice.XYZ(:,2,:)-lattice.XYZ(:,3,:)).^2,3)); %outboard
panelchords3=(panelchords1+panelchords2)/2; %Chord of each panel, CAUTION 
                                            %this is really camber line
                                            %length, so not really chord
                                            %for very cambered profiles

for i=1:indx1;			%Wing	
   for j=1:indx2;		%Partition
      lemma=[]; %local chord lemma vector.
      chordwisepanels=geo.nx(i,j)+geo.fnx(i,j); %number of panels chordwise on 
                                                %this partition 
      for k=1:geo.ny(i,j)                       %loop over panel strips.
          if geo.ny(i,j)~=0
              lemma=[lemma sum(panelchords3(1:chordwisepanels))];
              panelchords3=panelchords3((chordwisepanels+1):end);
              %size(panelchords3);
          end
      end  
      if geo.symetric(i)==1	%symmetric wings got two sides
         lc=[lc lemma lemma];
         panelchords3=panelchords3((chordwisepanels*geo.ny(i,j)+1):end);
      else
         lc=[lc lemma];
      end
          
   end
end
end

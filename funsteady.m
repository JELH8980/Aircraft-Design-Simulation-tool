function [output]=funsteady(lattice,geo,ref,Mach)



%% Transforming Tornado geometry to Tornagi
lattice2.colloc(:,1)=lattice.COLLOC(:,1);
lattice2.colloc(:,3)=lattice.COLLOC(:,2);  %Russian coordinate system y=z
lattice2.colloc(:,2)=lattice.COLLOC(:,3);

[a b c]=size(lattice.VORTEX);

try
if b>5

        lattice.VORTEX(:,2,:)=lattice.VORTEX(:,4,:); % JOp 0902
        lattice.VORTEX(:,3,:)=lattice.VORTEX(:,5,:); % JOp 0902
    
end
end

lattice2.vortex(:,1)=(lattice.VORTEX(:,2,1)+lattice.VORTEX(:,3,1))/2;
lattice2.vortex(:,3)=(lattice.VORTEX(:,2,2)+lattice.VORTEX(:,3,2))/2;
lattice2.vortex(:,2)=(lattice.VORTEX(:,2,3)+lattice.VORTEX(:,3,3))/2;

lattice2.normal(:,1)=-lattice.N(:,1);
lattice2.normal(:,3)=-lattice.N(:,2);
lattice2.normal(:,2)=-lattice.N(:,3);

x(1,:)=lattice.VORTEX(:,2,1);
x(2,:)=lattice.VORTEX(:,3,1);

y(1,:)=lattice.VORTEX(:,2,2);
y(2,:)=lattice.VORTEX(:,3,2);

z(1,:)=lattice.VORTEX(:,2,3);
z(2,:)=lattice.VORTEX(:,3,3);

dx=diff(x);
dy=diff(y);
dz=diff(z);


lattice2.semisp=sqrt(dz.^2+dy.^2)./2;
lattice2.vsweep=dx./lattice2.semisp;
cntr_pos=zeros(size(dx));

Sym=0;
Hground=40;


Ref.Swn=ref.S_ref;
Ref.MAC=ref.C_mac;
Ref.Spn=ref.b_ref;
Ref.Xcg=geo.CG(1);
Ref.Ycg=geo.CG(2);
Ref.Zcg=geo.CG(3);


%% Calling Tornagi functions
%% A=matr(lattice2,Sym,Hground,Mach,Ref);  
A = matrJO(lattice2,Sym,Hground,Mach,Ref);%Downwash influence matrix.
Am1=inv(A);                             %Inverting
task_id='alpha';
b=rhs(task_id,lattice2,Ref,cntr_pos);
g=A\b;

%% Adot=matr_dot(lattice2,Sym,Hground,Mach,Ref);
  Adot = matr_dotJO(lattice2,Sym,Hground,Mach,Ref);

gdot=-Am1*Adot*g/Ref.MAC;

[cy,cz,mx,my,mz]=aer_coef(lattice2,gdot,Ref,Sym);
% Russian notations
%output.cyda=cy;
%output.mzda=mz;
       
% American notations
%output.CZ_a_dot=-2*cy;
%output.Cm_a_dot=2*mz;

output(1)=-2*cy;
output(2)=2*mz;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%================ added 080729 JO: vectirized functions for speed
function A=matrJO(lattice,Sym,Hground,Mach,Ref)
% Calculation of matrix for steady task
% VECTORIZED
% Input:
% lattice - aircraft vortex lattice
% Sym - symmetry sign
% Hground - height above ground
% Mach - flow Mach number
% Ref - structure with aircraft reference data
%
% Output:
% A - resulting matrix
%

Hinf=5;
H=Hground*Ref.MAC;

n=length(lattice.semisp);
A=zeros(n,n);
switch Sym
    case 1 % Symmetrical solution
        if Hground >= Hinf
            % Far from ground
            for i=1:n
                % vectorize over j
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x =lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y =dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z =-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % symmetrical vortex
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)+lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=(vy+vy1).*lattice.normal(:,2)+(vz+vz1).*lattice.normal(:,3);
                %end
            end
        else
            % Account of ground effects
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % under ground vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % total from pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % symmetrical vortex
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)+lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                % under ground symmetrical vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_sta_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                % total from symmetrical pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % matrix coeeficient
                A(:,i)=(vy+vy1).*lattice.normal(:,2)+(vz+vz1).*lattice.normal(:,3);
                %end
            end
        end
    case 0  %  Solution without symmetry
        if Hground >= Hinf
            % Far from ground
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=vy.*lattice.normal(:,2)+vz.*lattice.normal(:,3);
                %end
            end
        else
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_sta_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % under ground vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_sta_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % total from pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=vy.*lattice.normal(:,2)+vz.*lattice.normal(:,3);
                %end
            end
        end
    otherwise
        disp('Unknown symmetry switch')
end
end %Function
%=============================================================
function A=matr_dotJO(lattice,Sym,Hground,Mach,Ref)
% Calculation of matrix for "dot" task
% VECTORIZED
% Input:
% lattice - aircraft vortex lattice
% Sym - symmetry sign
% Hground - height above ground
% Mach - flow Mach number
% Ref - structure with aircraft reference data
%
% Output:
% A - resulting matrix
%

Hinf=5;
H=Hground*Ref.MAC;

n=length(lattice.semisp);
A=zeros(n,n);
switch Sym
    case 1 % Symmetrical solution
        if Hground >= Hinf
            % Far from ground
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % symmetrical vortex
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)+lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=(vy+vy1).*lattice.normal(:,2)+(vz+vz1).*lattice.normal(:,3);
                %end
            end
        else
            % Account of ground effects
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % under ground vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % total from pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % symmetrical vortex
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)+lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                % under ground symmetrical vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_dyn_loop(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                % total from symmetrical pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                % matrix coeeficient
                A(:,i)=(vy+vy1).*lattice.normal(:,2)+(vz+vz1).*lattice.normal(:,3);
                %end
            end
        end
    case -1 %  Antisymmetrical solution
        disp('Antisymmetrical solution')
    case 0  %  Solution without symmetry
        if Hground >= Hinf
            % Far from graund
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x =lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=vy.*lattice.normal(:,2)+vz.*lattice.normal(:,3);
                %end
            end
        else
            for i=1:n
                %for j=1:n
                % j - number of collocation point
                % i - numer of vortex
                x=lattice.colloc(:,1)-lattice.vortex(i,1);
                dy=lattice.colloc(:,2)-lattice.vortex(i,2);
                dz=lattice.colloc(:,3)-lattice.vortex(i,3);
                y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy0,vz0]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % under ground vortex
                dy=dy+2*H;
                y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                [vy2,vz2]=sks_dyn_loop(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                % total from pair
                vy0=vy0-vy2; vz0=vz0-vz2;
                vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                A(:,i)=vy.*lattice.normal(:,2)+vz.*lattice.normal(:,3);
                %end
            end
        end
    otherwise
        disp('Unknown symmetry switch')
end
end%Function





%-===================================================================
function [vyvec,vzvec]=sks_sta_loop(xvec,yvec,zvec,l_2,tg,M)
% Calculation of steady downwash from horseshoe vortex
%
% Input:
% x - longitudinal coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% y - vertical coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% z - side coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% l_2 - semispan of horseshoe vortex
% tg - tangent of horseshoe vortex sweep
% M - flow Mach number
%
% Output:
% vy - vertical downwash
% vz - side downwash
%
vyvec = xvec;
vzvec = xvec;
xvec = xvec(:);
yvec = yvec(:);
zvec = zvec(:);

bet2=1-M*M;
bet2tg2=bet2+tg*tg;
for k = 1:length(xvec)
    x = xvec(k);
    y = yvec(k);
    z = zvec(k);
    %------
    %l_2
    %z
    z1=-l_2-z;
    z2= l_2-z;
    x1=x-tg*z;
    y2=y*y;
    r1=x1*x1+bet2tg2*y2;
    r3=y2+z1*z1;
    r4=y2+z2*z2;
    a1=x1-z1*tg;
    a2=x1-z2*tg;
    r5=sqrt(a1*a1+bet2*r3);
    r6=sqrt(a2*a2+bet2*r4);
    if r3==0
        disp('asas')
    end
    vy=-z2/r4*((x1-z2*tg)/r6+1)+z1/r3*((x1-z1*tg)/r5+1);
    vy=vy-x1/r1*((bet2tg2*z2-x1*tg)/r6-(bet2tg2*z1-x1*tg)/r5);
    vy=.5*vy;
    vz=1/r3-1/r4-(x1-z2*tg)/r6/r4+(x1-z1*tg)/r5/r3-tg/r1*((bet2tg2*z2-x1*tg)/r6-(bet2tg2*z1-x1*tg)/r5);
    vz=.5*y*vz;
    vyvec(k) = vy;
    vzvec(k)=vz;
end
end
%==============================================================
function [VYvec,VZvec]=sks_dyn_loop(Xvec,Yvec,Zvec,AL2,TG,AM)
% Calculation of dynamic downwash from horseshoe vortex
%
% Input:
% X - longitudinal coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% Y - vertical coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% Z - side coordinate of the collocation point with respect to the
%     center of horseshoe vortex
% AL2 - semispan of horseshoe vortex
% TG - tangent of horseshoe vortex sweep
% AM - flow Mach number
%
% Output:
% VY - vertical downwash
% VZ - side downwash
%
VYvec = Xvec;
VZvec = Xvec;
Xvec = Xvec(:);
Yvec = Yvec(:);
Zvec = Zvec(:);

EPS=.1*AL2*AL2;
B=1.-AM*AM;
BETSQ=B+TG*TG;
BET=sqrt(BETSQ);
for k = 1:length(Xvec);
    X = Xvec(k);
    Y = Yvec(k);
    Z = Zvec(k);
    Z1=-AL2-Z;
    Z2=AL2-Z;
    X1=X-Z*TG;
    if   (Y < 0 || Y > 0)
        YSQ=Y*Y;
        R1=Z1*Z1+YSQ;
        R2=Z2*Z2+YSQ;
        R3=X1-Z1*TG;
        R4=X1-Z2*TG;
        R0=X1*X1+YSQ*BETSQ;
        R5=sqrt(R3*R3+B*R1);
        R6=sqrt(R4*R4+B*R2);
        R7=B*Z1-R3*TG;
        R8=B*Z2-R4*TG;
        if (R4 <= 0)
            R10=R6-R4;
        else
            R10=B*R2/(R6+R4);
        end
        if (R3 <= 0)
            R10=R10/(R5-R3);
        else
            R10=R10*(R5+R3)/(B*R1);
        end
        VY=(1.+R4/R6)*(R4*Z2+R2*TG)/R2-(1.+R3/R5)*(R3*Z1+R1*TG)/R1...
            +(R0+AM*AM*YSQ)/R0*(R8/R6-R7/R5)-(1.+TG*TG)/BET*log((BET*R6...
            +R8)/(BET*R5+R7))+TG*log(R10);
        VZ=Y*((R6+R4)/R2-(R5+R3)/R1+AM*AM*X1/R0*(R4/R6-R3/R5))+...
            TG*(atan(Z2/Y)-atan(Z1/Y)+atan((R4*Z2-R2*TG)/(Y*R6))-...
            atan((R3*Z1-R1*TG)/(Y*R5)));
        VY=.5*VY;
        VZ=.5*VZ;
    else
        R1=Z1*Z1*B;
        R2=Z2*Z2*B;
        R3=X1-Z1*TG;
        R4=X1-Z2*TG;
        R5=sqrt(R3*R3+R1);
        R6=sqrt(R4*R4+R2);
        R7=B*Z1-R3*TG;
        R8=B*Z2-R4*TG;
        if (R4 <= 0)
            R10=R6-R4;
        else
            R10=R2/(R6+R4);
        end
        if (R3 <= 0)
            R10=R10/(R5-R3);
        else
            R10=R10*(R5+R3)/R1;
        end
        VY=(X1+R6)/Z2-(X1+R5)/Z1+TG*log(R10);
        R0=sqrt(B)*abs(X1);
        if (R0 < EPS)
            VY=VY-sign(Z2)*(1.+TG*TG)/BET*log(abs(Z2/Z1));
        else
            VY=VY-(1.+TG*TG)/BET*(sign(R8)*log((BET*R6+abs(R8))/...
                R0)-sign(R7)*log((BET*R5+abs(R7))/R0));
        end
        VY=.5*VY;
        VZ=0.;
    end
    VYvec(k) = VY;
    VZvec(k) = VZ;
end






end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=matr(lattice,Sym,Hground,Mach,Ref)
% Calculation of matrix for steady task
%
% Input:
% lattice - aircraft vortex lattice
% Sym - symmetry sign
% Hground - height above ground
% Mach - flow Mach number
% Ref - structure with aircraft reference data
%
% Output:
% A - resulting matrix
%

Hinf=5;
H=Hground*Ref.MAC;

n=length(lattice.semisp);
A=zeros(n,n);
switch Sym
    case 1 % Symmetrical solution
       if Hground >= Hinf
           % Far from ground
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % symmetrical vortex
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)+lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                    vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=(vy+vy1)*lattice.normal(j,2)+(vz+vz1)*lattice.normal(j,3);
                end
            end
       else
            % Account of ground effects
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % under ground vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % total from pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % symmetrical vortex
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)+lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    % under ground symmetrical vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_sta(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    % total from symmetrical pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                    vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % matrix coeeficient
                    A(j,i)=(vy+vy1)*lattice.normal(j,2)+(vz+vz1)*lattice.normal(j,3);
                end
            end
       end           
    case 0  %  Solution without symmetry
       if Hground >= Hinf
           % Far from ground
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=vy*lattice.normal(j,2)+vz*lattice.normal(j,3);
                end
            end
       else
             for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % under ground vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_sta(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % total from pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=vy*lattice.normal(j,2)+vz*lattice.normal(j,3);
                end
            end
       end
    otherwise
        disp('Unknown symmetry switch')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=rhs(task_id,lattice,Ref,cntr_pos)
% Calculation of right hand sides for main system of linear equations
% 
% Input:
% task_id - identificator of solving task
% lattice - aircraft vortex lattice 
% Ref - structure with aircraft reference data
% cntr_pos - 
%
% Output:
% b - vector of right hand sides
%
N=length(lattice.semisp);
b=zeros(N,1);
switch task_id
    case 'zero'
        b(:)=-2*pi*lattice.normal(:,1);
    case 'alpha'
        b(:)=-2*pi*lattice.normal(:,2);
    case 'beta'
        b(:)=2*pi*lattice.normal(:,3);
    case 'roll'
        b(:)=(lattice.colloc(:,2)-Ref.Ycg).*lattice.normal(:,3)-...
             (lattice.colloc(:,3)-Ref.Zcg).*lattice.normal(:,2);
        b=4*pi*b/Ref.Spn;
    case 'yaw'
        b(:)=4*pi*(lattice.colloc(:,1)-Ref.Xcg).*lattice.normal(:,3)/Ref.Spn;
    case 'pitch'
        b(:)=-2*pi*(lattice.colloc(:,1)-Ref.Xcg).*lattice.normal(:,2)/Ref.MAC;
    case 'control'
        b(:)=-2*pi*cntr_pos(:);     
    otherwise
        disp(strcat('RHS: unknown task:_',task_id))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=matr_dot(lattice,Sym,Hground,Mach,Ref)
% Calculation of matrix for "dot" task
%
% Input:
% lattice - aircraft vortex lattice
% Sym - symmetry sign
% Hground - height above ground
% Mach - flow Mach number
% Ref - structure with aircraft reference data
%
% Output:
% A - resulting matrix
%

Hinf=5;
H=Hground*Ref.MAC;

n=length(lattice.semisp);
A=zeros(n,n);
switch Sym
    case 1 % Symmetrical solution
       if Hground >= Hinf
           % Far from ground
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % symmetrical vortex
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)+lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                    vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=(vy+vy1)*lattice.normal(j,2)+(vz+vz1)*lattice.normal(j,3);
                end
            end
       else
            % Account of ground effects
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % under ground vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % total from pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % symmetrical vortex
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)+lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    % under ground symmetrical vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_dyn(x,y,z,lattice.semisp(i),-lattice.vsweep(i),Mach);
                    % total from symmetrical pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy1=vy0*lattice.normal(i,2)+vz0*lattice.normal(i,3);
                    vz1=-vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    % matrix coeeficient
                    A(j,i)=(vy+vy1)*lattice.normal(j,2)+(vz+vz1)*lattice.normal(j,3);
                end
            end
       end           
    case -1 %  Antisymmetrical solution
        disp('Antisymmetrical solution')
    case 0  %  Solution without symmetry
       if Hground >= Hinf
           % Far from graund
            for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=vy*lattice.normal(j,2)+vz*lattice.normal(j,3);
                end
            end
       else
             for i=1:n
                for j=1:n
                    % j - number of collocation point
                    % i - numer of vortex
                    x=lattice.colloc(j,1)-lattice.vortex(i,1);
                    dy=lattice.colloc(j,2)-lattice.vortex(i,2);
                    dz=lattice.colloc(j,3)-lattice.vortex(i,3);
                    y=dy*lattice.normal(i,2)+dz*lattice.normal(i,3);
                    z=-dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy0,vz0]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % under ground vortex
                    dy=dy+2*H;
                    y=dy*lattice.normal(i,2)-dz*lattice.normal(i,3);
                    z=dy*lattice.normal(i,3)+dz*lattice.normal(i,2);
                    [vy2,vz2]=sks_dyn(x,y,z,lattice.semisp(i),lattice.vsweep(i),Mach);
                    % total from pair
                    vy0=vy0-vy2; vz0=vz0-vz2;
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    vy=vy0*lattice.normal(i,2)-vz0*lattice.normal(i,3);
                    vz=vy0*lattice.normal(i,3)+vz0*lattice.normal(i,2);
                    A(j,i)=vy*lattice.normal(j,2)+vz*lattice.normal(j,3);
                end
            end
       end
    otherwise
        disp('Unknown symmetry switch')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cy,cz,mx,my,mz]=aer_coef(lattice,circ,Ref,Sym)
% calculation of aerodynamic loads
%
% Output:
% cy - vertical force coefficient
% cz - side force coefficient
% mx - rolling moment coefficient
% my - yawing moment coefficient
% mz - pitching moment coefficient
%
% Input:
% lattice - vortex lattice structure
% circ - circulations of vortex lattice
% Ref - structure with referensial data
% Sym - symmetry sign

N=length(lattice.vsweep);
switch Sym
    case 1
        cy=0; cz=0; mx=0; my=0; mz=0;   
        for i=1:N
           dcy=2*circ(i)*lattice.semisp(i)*lattice.normal(i,2);
           cy=cy+dcy;
           mz=mz-dcy*(lattice.vortex(i,1)-Ref.Xcg);
        end
        cy=4*cy/Ref.Swn;
        mz=4*mz/Ref.Swn/Ref.MAC;
    case -1
        cy=0; cz=0; mx=0; my=0; mz=0;  
    case 0
        cy=0; cz=0; mx=0; my=0; mz=0;  
        for i=1:N
           dx=lattice.vortex(i,1)-Ref.Xcg;
           dy=lattice.vortex(i,2)-Ref.Ycg;
           dz=lattice.vortex(i,3)-Ref.Zcg;
           dp=2*circ(i)*lattice.semisp(i);
           dcy=dp*lattice.normal(i,2);
           dcz=dp*lattice.normal(i,3);
           cy=cy+dcy;
           cz=cz+dcz;
           mx=mx+dcz*dy-dcy*dz;
           my=my+dcz*dx;
           mz=mz-dcy*dx;
        end
        cy=2*cy/Ref.Swn;
        cz=2*cz/Ref.Swn;
        mx=2*mx/Ref.Swn/Ref.Spn;
        my=2*my/Ref.Swn/Ref.Spn;
        mz=2*mz/Ref.Swn/Ref.MAC;
    otherwise
        disp(strcat('AER_COEF: unknown Sym:_',num2str(Sym)))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vy,vz]=sks_sta(x,y,z,l_2,tg,M)
% Calculation of steady downwash from horseshoe vortex
%
% Input:
% x - longitudinal coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% y - vertical coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% z - side coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% l_2 - semispan of horseshoe vortex
% tg - tangent of horseshoe vortex sweep
% M - flow Mach number 
%
% Output:
% vy - vertical downwash
% vz - side downwash
%

bet2=1-M*M;
bet2tg2=bet2+tg*tg;
%l_2
%z
z1=-l_2-z;
z2= l_2-z;
x1=x-tg*z;
y2=y*y;
r1=x1*x1+bet2tg2*y2;
r3=y2+z1*z1;
r4=y2+z2*z2;
a1=x1-z1*tg;
a2=x1-z2*tg;
r5=sqrt(a1*a1+bet2*r3);
r6=sqrt(a2*a2+bet2*r4);
if r3==0
    disp('asas')
end
vy=-z2/r4*((x1-z2*tg)/r6+1)+z1/r3*((x1-z1*tg)/r5+1);



vy=vy-x1/r1*((bet2tg2*z2-x1*tg)/r6-(bet2tg2*z1-x1*tg)/r5);
vy=.5*vy;
vz=1/r3-1/r4-(x1-z2*tg)/r6/r4+(x1-z1*tg)/r5/r3-tg/r1*((bet2tg2*z2-x1*tg)/r6-(bet2tg2*z1-x1*tg)/r5);
vz=.5*y*vz;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VY,VZ]=sks_dyn(X,Y,Z,AL2,TG,AM)
% Calculation of dynamic downwash from horseshoe vortex
%
% Input:
% X - longitudinal coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% Y - vertical coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% Z - side coordinate of the collocation point with respect to the 
%     center of horseshoe vortex
% AL2 - semispan of horseshoe vortex
% TG - tangent of horseshoe vortex sweep
% AM - flow Mach number 
%
% Output:
% VY - vertical downwash
% VZ - side downwash
%
       EPS=.1*AL2*AL2;
       B=1.-AM*AM;
       BETSQ=B+TG*TG;
       BET=sqrt(BETSQ);
       Z1=-AL2-Z;
       Z2=AL2-Z;
       X1=X-Z*TG;
if   (Y < 0 || Y > 0) 
       YSQ=Y*Y;
       R1=Z1*Z1+YSQ;
       R2=Z2*Z2+YSQ;
       R3=X1-Z1*TG;
       R4=X1-Z2*TG;
       R0=X1*X1+YSQ*BETSQ;
       R5=sqrt(R3*R3+B*R1);
       R6=sqrt(R4*R4+B*R2);
       R7=B*Z1-R3*TG;
       R8=B*Z2-R4*TG;
       if (R4 <= 0)
           R10=R6-R4;
       else
           R10=B*R2/(R6+R4);
       end
       if (R3 <= 0)
           R10=R10/(R5-R3);
       else
           R10=R10*(R5+R3)/(B*R1);
       end
       VY=(1.+R4/R6)*(R4*Z2+R2*TG)/R2-(1.+R3/R5)*(R3*Z1+R1*TG)/R1...
           +(R0+AM*AM*YSQ)/R0*(R8/R6-R7/R5)-(1.+TG*TG)/BET*log((BET*R6...
           +R8)/(BET*R5+R7))+TG*log(R10);
       VZ=Y*((R6+R4)/R2-(R5+R3)/R1+AM*AM*X1/R0*(R4/R6-R3/R5))+...
           TG*(atan(Z2/Y)-atan(Z1/Y)+atan((R4*Z2-R2*TG)/(Y*R6))-...
           atan((R3*Z1-R1*TG)/(Y*R5)));
        VY=.5*VY;
        VZ=.5*VZ;
 else
       R1=Z1*Z1*B;
       R2=Z2*Z2*B;
       R3=X1-Z1*TG;
       R4=X1-Z2*TG;
       R5=sqrt(R3*R3+R1);
       R6=sqrt(R4*R4+R2);
       R7=B*Z1-R3*TG;
       R8=B*Z2-R4*TG;
       if (R4 <= 0)
           R10=R6-R4;
       else
           R10=R2/(R6+R4);
       end
       if (R3 <= 0)
           R10=R10/(R5-R3);
       else
           R10=R10*(R5+R3)/R1;
       end
       VY=(X1+R6)/Z2-(X1+R5)/Z1+TG*log(R10);
       R0=sqrt(B)*abs(X1);
       if (R0 < EPS)
           VY=VY-sign(Z2)*(1.+TG*TG)/BET*log(abs(Z2/Z1));
       else
           VY=VY-(1.+TG*TG)/BET*(sign(R8)*log((BET*R6+abs(R8))/...
               R0)-sign(R7)*log((BET*R5+abs(R7))/R0));
       end
       VY=.5*VY;
       VZ=0.;
end
end%function



















 
            

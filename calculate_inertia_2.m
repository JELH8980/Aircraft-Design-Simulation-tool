function [model] = calculate_inertia_2(model)
% CALCULATE_INERTIA_2 - Computes the mass properties of an aircraft model.
%
% INPUT:
%   model - Struct containing geometric definitions of the aircraft.
%
% OUTPUT:
%   model - Updated struct with computed mass and inertia properties.
%
% USER INPUT:
%   BEW  - Basic empty weight of the aircraft [kg].
%
% PROCESS:
% - Iterates over defined bodies in the model.
% - Checks for required geometric fields.
% - Computes total volume and estimates uniform density.
% - Computes mass moments and inertia tensor elements.
% - Determines the aircraft's center of gravity (CG).
% - Adjusts inertia values to reference the body frame.
%
% VARIABLES:
%   m           - Total aircraft mass [kg]
%   V           - Total volume of all bodies [m^3]
%   density     - Estimated uniform density [kg/m^3]
%   XB, YB, ZB  - CG coordinates in the body frame [m]
%   Ixx, Iyy, Izz, Ixz - Inertia tensor elements [kg*m^2]
%
% CALCULATIONS:
% - Volume estimation: Sum of body volumes with symmetry scaling.
% - Mass distribution: Elemental volume * density.
% - Inertia: Summed contributions from elements.
% - CG shift correction: Adjusts inertia relative to body reference frame.
%
% OUTPUT STRUCTURE:
% - model.parameters.m      : Mass [kg]
% - model.parameters.lB     : CG x-coordinate [m]
% - model.parameters.hB     : CG z-coordinate [m]
% - model.parameters.Ixx    : Inertia about x-axis [kg*m^2]
% - model.parameters.Iyy    : Inertia about y-axis [kg*m^2]
% - model.parameters.Izz    : Inertia about z-axis [kg*m^2]
% - model.parameters.Ixz    : Product of inertia [kg*m^2]
% - model.geo.CG           : CG position relative to reference [m]
%
% METHOD:
% - This function is marked under model.parameters.method.inertia.


bodies_names = fields(model.geo.bodies);
no_bodies = numel(bodies_names);


m =  str2double(input('    Enter BEW: ', 's'));

% Initialize total moment accumulators
XB_moment = 0;
YB_moment = 0;
ZB_moment = 0;

m_acc = 0;

Ixx = 0;

Iyy = 0;

Izz = 0;

Ixz = 0;

V = 0;

for i = 1:no_bodies

    body_struct = model.geo.bodies.(bodies_names{i});


    if any([~isfield(body_struct, 'elements'), ...
            ~isfield(body_struct, 'volume'), ...
            ~isfield(body_struct, 'points'), ...
            ~isfield(body_struct, 'convexhull')])
             disp(append(bodies_names{i}, ' skipped (not properly defined).'))

    else

        sym_factor = 1 + body_struct.symmetry;

        V = V + body_struct.volume*sym_factor;

    end

end

density = m/V;

for i = 1:no_bodies   

    body_struct = model.geo.bodies.(bodies_names{i});
    
    if any([~isfield(body_struct, 'elements'), ...
            ~isfield(body_struct, 'volume'), ...
            ~isfield(body_struct, 'points'), ...
            ~isfield(body_struct, 'convexhull')])
             disp(append(bodies_names{i}, ' skipped (not properly defined).'))

    else

        % Calculating Moment of Inertia and center of gravity
    
        no_elements = height(body_struct.elements.centroids);
    
        dv = body_struct.elements.volume;
    
        dm = dv*density;
    
        for j = 1:no_elements
            centroid = body_struct.elements.centroids(j, :); % Corrected indexing
            
            % Compute relative coordinates
            XBel = centroid(1);
            YBel = centroid(2);
            ZBel = centroid(3);
            
    
            % Moment of Inertia calculations
            Ixx = Ixx +  dm * (YBel^2 + ZBel^2) * sym_factor;
            Iyy = Iyy +  dm * (XBel^2 + ZBel^2) * sym_factor;
            Izz = Izz +  dm * (XBel^2 + YBel^2) * sym_factor;
            Ixz = Ixz +  dm * (XBel * ZBel) * sym_factor;

    
            % Accumulate weighted moments
            XB_moment = XB_moment + dm * XBel * sym_factor;
            YB_moment = YB_moment + dm * YBel * (2 - sym_factor); 
            ZB_moment = ZB_moment + dm * ZBel * sym_factor;

            % Accumulate total mass
            m_acc = m_acc + dm * sym_factor; 

        end

    end
    
end

% Compute final center of mass positions
if m_acc > 0
    XB = XB_moment / m_acc;
    YB = YB_moment / m_acc;
    ZB = ZB_moment / m_acc;
else
    XB = NaN;
    YB = NaN;
    ZB = NaN;
end

% Translating inertia to B frame



model.parameters.m = m;

model.parameters.lB = XB;
model.parameters.hB = ZB;
lR = model.parameters.lR;
hR = model.parameters.hR;

model.parameters.Ixx = Ixx;
model.parameters.Iyy = Iyy - m*((lR-XB)^2+(ZB-hR)^2);
model.parameters.Izz = Izz - m*(XB-lR)^2;
model.parameters.Ixz = Ixz - m*(XB-lR)*(ZB-lR);


model.parameters.method.inertia = 'calculate_inertia_2';

model.geo.CG(1) = model.parameters.lB - model.parameters.lR;
model.geo.CG(2) = 0;
model.geo.CG(3) = model.parameters.hB - model.parameters.hR;

end
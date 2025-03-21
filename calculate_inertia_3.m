function [model] = calculate_inertia_3(model)
% CALCULATE_INERTIA_3 - Computes the mass properties of an aircraft model 
% considering material properties.
%
% INPUT:
%   model - Struct containing geometric definitions and material properties.
%
% OUTPUT:
%   model - Updated struct with computed mass and inertia properties.
%
% PROCESS:
% - Iterates over defined bodies in the model.
% - Checks for required geometric and material fields.
% - Computes total mass based on material density fractions.
% - Computes mass moments and inertia tensor elements.
% - Determines the aircraft's center of gravity (CG).
% - Adjusts inertia values relative to the body reference frame.
%
% VARIABLES:
%   m           - Total aircraft mass [kg]
%   V_body      - Volume of an individual body [m^3]
%   density     - Density of each material [kg/m^3]
%   fraction    - Material fraction in body composition [-]
%   weighted_sum- Weighted sum of material densities [-]
%   XB, YB, ZB  - CG coordinates in the body frame [m]
%   Ixx, Iyy, Izz, Ixz - Inertia tensor elements [kg*m^2]
%
% CALCULATIONS:
% - Volume estimation: Sum of body volumes with symmetry scaling.
% - Mass estimation: Summed contributions from material densities.
% - Inertia: Elemental volume and mass properties determine moments.
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


% Calculating mass

m = 0;

XB_moment = 0;
YB_moment = 0;
ZB_moment = 0;

m_acc = 0;

Ixx = 0;

Iyy = 0;

Izz = 0;

Ixz = 0;

for i = 1:no_bodies

    body_struct = model.geo.bodies.(bodies_names{i});



    if any([~isfield(body_struct, 'elements'), ...
            ~isfield(body_struct, 'volume'), ...
            ~isfield(body_struct, 'points'), ...
            ~isfield(body_struct, 'convexhull'), ...
            ~isfield(body_struct, 'materials')])
             disp(append(bodies_names{i}, ' skipped (not properly defined).'))

    else
        
        
        % Calculating mass


        sym_factor = 1 + body_struct.symmetry;

        V_body = body_struct.volume;

        material_names = fields(body_struct.materials);

        weighted_sum = 0;

        for j=1:numel(material_names)
            
           
            material_name = material_names{j};

            density = str2double(body_struct.materials.(material_name).density);

            if ~isfield(body_struct.materials.(material_name), 'fraction')
                disp(append(bodies_names{i}, ' skipped (not properly defined).'))

            else

                fraction = body_struct.materials.(material_name).fraction;

                weighted_sum = weighted_sum + density*fraction;


            end

   

        end

        m_body = V_body*weighted_sum;

        m = m + m_body;

        % Calculating Moment of Inertia and center of gravity

        XR = model.parameters.lR;
        YR = 0;
        ZR = model.parameters.hR;

        no_elements = height(body_struct.elements.centroids);

        dv = body_struct.elements.volume;

        dx = body_struct.elements.dx;

        dy = body_struct.elements.dy;

        dz = body_struct.elements.dz;

        dm = dv*weighted_sum;

        for j = 1:no_elements
            
            centroid = body_struct.elements.centroids(j, :); % Corrected indexing
            
            % Compute relative coordinates
            XBel = centroid(1) + XR;
            YBel = centroid(2) + YR;
            ZBel = centroid(3) + ZR;
    
            % Moment of Inertia calculations
            Ixx = Ixx + (1/12 * dm * (dy^2 + dz^2) + dm * (YBel^2 + ZBel^2)) * sym_factor;
            Iyy = Iyy + (1/12 * dm * (dx^2 + dz^2) + dm * (XBel^2 + ZBel^2));
            Izz = Izz + (1/12 * dm * (dx^2 + dy^2) + dm * (XBel^2 + YBel^2)) * sym_factor;
            Ixz = Ixz + dm * (XBel * ZBel) * sym_factor;
    
            % Accumulate weighted moments
            XB_moment = XB_moment + dm * XBel * sym_factor;
            YB_moment = YB_moment + dm * YBel * (2- sym_factor); 
            ZB_moment = ZB_moment + dm * ZBel * sym_factor;

            % Accumulate total mass
            m_acc = m_acc + dm * sym_factor; 

        end
        
    end


end

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

model.parameters.Ixx = Ixx + m*(ZB^2+YB^2);
model.parameters.Iyy = Iyy + m*(XB^2+YB^2);
model.parameters.Izz = Izz + m*(XB^2+YB^2);
model.parameters.Ixz = Ixz + m*XB*ZB;

model.parameters.lB = XB;
model.parameters.hB = ZB;

model.parameters.method.inertia = 'calculate_inertia_3';

model.geo.CG(1) = model.parameters.lB - model.parameters.lR;
model.geo.CG(2) = 0;
model.geo.CG(3) = model.parameters.hB - model.parameters.hR;

end
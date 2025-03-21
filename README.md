# Aircraft Design Simulation Tool



The **Aircraft Design Simulation Tool** is a MATLAB-based software package developed to support the design, analysis, and simulation of aerospace vehicles. It provides an infrastructure for creating geometry, generating aerodynamic data, solving for equilibirium conditions, studying the small disturbance models and also has simulating the model.

More specifically the toolkit enables the following capabilities:

- **Design Geometry for VLM**: Define aircraft geometry suitable for Vortex Lattice Method (VLM) aerodynamic analysis (Using Redhammer Tornado's solver).
  
- **Specify Thrustline, Engine Model, and Control Surfaces**: Configure propulsion system parameters (thrustline and engine model) and control surfaces (e.g., elevator, aileron, rudder).
  
- **Define Batch Calculations**: Calculate stability and control derivatives across a large multidimensional grid (`ndgrid`) of input states.
  
- **Solve for Straight Flight Conditions Using Optimization**: Determine trimmed straight flight conditions through optimization techniques.
  
- **Organize Results**: Store derivatives and straight flight condition data within a filesystem, referenced to the aircraft model.
  
- **Perform Modal Analysis**: Conduct modal analysis at any identified straight flight condition, completed within seconds.
  
- **Begin Control Design**: Initiate control system design using a Linear Time-Invariant (LTI) state-space formulation, followed by testing on a nonlinear plant.
  
- **Simulate with Simulink and FlightGear**: Perform simulations using Simulink, with optional integration of FlightGear and Logitech Extreme 3D Pro joystick support (requires Simulink configuration).
  
- **Animate and Visualize**: Animate flight dynamics in real or accelerated time.


# Software Infrastructure
![image](https://github.com/user-attachments/assets/3f70803a-c80f-44ef-8bac-027fbd1041d3)

# Geometry Parametrization

![image](https://github.com/user-attachments/assets/629743c4-7256-4e61-ac21-7e818b7c4cc6)

# Solving for SF conditions

![image](https://github.com/user-attachments/assets/263b1d50-e513-4771-984e-066e4899d3eb)


# Modal Analysis

![image](https://github.com/user-attachments/assets/46a2ce38-4094-455b-8945-515451ebe7cb)


# Visualization using Flight Gear
![image](https://github.com/user-attachments/assets/fcda63f0-0031-4193-9238-da0601250adc)

# Simulation output from predefined reference

![image](https://github.com/user-attachments/assets/4b1eb3ba-2eee-4677-8f79-f9367884cc43)

# State Trajectory Visualization

![image](https://github.com/user-attachments/assets/26c4b895-5116-4519-8ef9-835f9458285d)




## Future work includes:
- **Optimization-based solver** for stationary turn, pull-up and pull-down.
- **Measurement models** for typical fixed wing drone hardware configurations to allow for observer design and development.
- **3D-dubins** trajectory generation constrained to solved for trim conditions.
- **Refined infrastructure** for seamless testing and implementation of controllers in Simulink.

### Comments:
- **Ground effect** This software does not provide any explicit tool for including ground effect in your model, this is due to the limited CFD knowledge of the author to this code, if Redhammer Tornado comes with an update including such capabilities then they will be implemented. 
- **Old syntax** Due to some old syntax in the Tornado code errors may be raised in future relesaes of MatLab, the code-lines in question trigger warnings in V2024a/b while no warnings occures in earlier versions. The code is made in Simulink for V2024 so to exploit all capabiltiies in earlier versions of MatLab the code must be exported to an older version.

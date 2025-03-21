# Aerospace Vehicle Analysis Toolkit

The **Aerospace Vehicle Analysis Toolkit** is a MATLAB-based software package developed to support the design, analysis, and simulation of aerospace vehicles. It provides a structured set of tools for aerodynamic modeling, flight dynamics, and control system development.

The toolkit enables the following capabilities:

- **Design Geometry for VLM**: Define aircraft geometry suitable for Vortex Lattice Method (VLM) aerodynamic analysis.
- **Specify Thrustline, Engine Model, and Control Surfaces**: Configure propulsion system parameters (thrustline and engine model) and control surfaces (e.g., elevator, aileron, rudder).
- **Define Batch Calculations**: Calculate stability and control derivatives across a large multidimensional grid (`ndgrid`) of input states.
- **Solve for Straight Flight Conditions Using Optimization**: Determine trimmed straight flight conditions through optimization techniques.
- **Organize Results**: Store derivatives and straight flight condition data within a filesystem, referenced to the aircraft model.
- **Perform Modal Analysis**: Conduct modal analysis at any identified straight flight condition, completed within seconds.
- **Begin Control Design**: Initiate control system design using a Linear Time-Invariant (LTI) state-space formulation, followed by testing on a nonlinear plant.
- **Simulate with Simulink and FlightGear**: Perform simulations using Simulink, with optional integration of FlightGear and Logitech Extreme 3D Pro joystick support (requires Simulink configuration).
- **Animate and Visualize**: Animate flight dynamics in real or accelerated time and display trajectories in geodetic coordinates.

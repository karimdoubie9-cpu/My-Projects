# Adaptive Control of a Hexacopter for Package Delivery

Final project for the **Adaptive and Autonomous Aerospace Systems** course (M.Sc. Space Engineering, Politecnico di Milano).

Team: Antonio Fraccastoro, Karim Ahmed Mohamed Doubie, Mohamed Ghoraba.

## About

Design, adaptive augmentation and mission testing of a controller for a hexacopter that delivers packages. A baseline cascaded controller is built first. A predictor-based model reference adaptive controller (PBMRAC) is then added to handle the unknown package mass, rotor faults and wind. Finally the vehicle flies a full delivery mission through an urban environment. Everything is implemented and simulated in MATLAB / Simulink.

## Contents

| Task | Topic | Main points |
|---|---|---|
| 1 | Model implementation | Simulink model of the hexacopter (actuator dynamics, rigid-body dynamics, position and attitude kinematics), hover trim and rotor-speed allocation |
| 2 | Baseline controller | Cascaded architecture (PID position loop, geometric attitude control on SO(3)), loop-shaping tuning with 2-DOF loops, successive loop closure, setpoint and helical trajectory tracking, reduced propeller effectiveness, force steps and wind gust |
| 3 | Adaptive augmentation | Package modelled as a rigid body, predictor-based MRAC with e-modification, projection and saturation hedging, gain selection, effectiveness estimation, robustness to rotor fault and gust, Monte Carlo validation on 500 sampled vehicles |
| 4 | Delivery mission | Occupancy map with obstacle inflation, RRT* route planning with restarts and path shortcutting, minimum-snap trajectory with time optimization and speed limit, yaw generation, disturbed 2-D and 3-D flights comparing baseline and PBMRAC |

## Tools used

MATLAB / Simulink (`sisotool` for loop shaping), RRT* and minimum-snap planning.

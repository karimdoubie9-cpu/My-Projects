# LADEE Reverse Engineering Study

Final report for the **Space System Engineering and Operations** course (M.Sc. Space Engineering, Politecnico di Milano, A.Y. 2024-2025).

Team 15: Arturo Calderazzi, Riccardo Giovenzana, Matteo Candita, Karim Ahmed Mohamed Doubie, Riccardo Cislaghi, Mohamed Ghoraba.

Course staff: Prof. M. Lavagna, L. Bianchi, M. Bussolino.

## About

A reverse engineering study of NASA's **Lunar Atmosphere and Dust Environment Explorer (LADEE)**, a robotic mission that orbited the Moon from September 2013 to April 2014 to study its tenuous atmosphere and dust exosphere. LADEE also demonstrated a laser communication system (LLCD) and the low-cost NASA Ames Modular Common Spacecraft Bus (MCSB).

Starting from public literature, the team reconstructed the mission and spacecraft design and then re-sized each subsystem from the mission requirements. The results were compared with real LADEE values.

## Contents

| Section | Topic |
|---|---|
| 1 | Mission introduction: goals, drivers, functional analysis, requirements, phases, ConOps, payload |
| 2 | Mission analysis: trajectory, orbit segments, delta-V budget |
| 3 | Propulsion: bipropellant MMH/NTO system, reverse sizing of tanks and pressurant |
| 4 | TMTC: S-band architecture, ground segment, link budget |
| 5 | AOCS: sensors, actuators, modes, pointing budget, disturbance torques, actuator sizing |
| 6 | Thermal control: passive/active design, temperature ranges, single-node and 18-node models |
| 7 | EPS: body-mounted solar arrays, Li-ion battery, power budget, sizing |
| 8 | OBDH: RAD750-based avionics, data bus, onboard data budget |
| 9 | Overall comments: structure, mechanisms, launcher interface (Minotaur V) |
| 10 | Mass budget: sizing results vs. literature |

## Tools used

MATLAB / Simulink, LaTeX (report).

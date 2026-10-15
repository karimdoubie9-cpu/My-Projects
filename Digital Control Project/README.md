# Digital Control of a Flexible Transmission

Project for the **Digital Control Technology for Aeronautics** course (M.Sc. Space Engineering, Politecnico di Milano, A.Y. 2025-26).

Authors: Riccardo Cislaghi, Karim Ahmed Mohamed Doubie.

## About

Design and validation of a digital control system for a flexible aeronautical actuator: a motor and a load connected through a gearbox and an elastic transmission, which gives a fourth-order plant with a lightly damped flexible mode. The work has two parts. Part I designs the controllers in the Z-domain assuming an ideal digital implementation. Part II adds the real signal chain (sensors, amplification, filtering, ADC and DAC) in a hybrid Simulink simulation.

## Contents

| Part | Section | Topic |
|---|---|---|
| I | 1 | Methodology and design criteria: θ1/ωm cascade architecture (outer position loop, inner speed loop) and performance requirements |
| I | 2 | System analysis: modelling, controllability and observability, structural dynamics (resonant poles and anti-resonant zeros) |
| I | 3 | Choice of the sampling time |
| I | 4 | Controller design: zero-order-hold discretization, inner speed loop and outer position loop controllers |
| I | 5 | Stability margins of the resulting system |
| I | 6 | Performance verification and robustness analysis (sensitivity function) |
| II | 7 | Signal acquisition: tachometer and encoder sensors, amplifier, anti-aliasing filter, ADC and DAC sizing |
| II | 8 | Complete Simulink hybrid simulation |
| II | 9 | Performance evaluation against the design requirements |
| II | 10 | Conclusions |

## Tools used

MATLAB / Simulink (`c2d`, `sisotool`), LaTeX (report).

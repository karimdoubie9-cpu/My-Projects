# SGN Assignment #1: Spacecraft Guidance and Navigation

Report for the **Spacecraft Guidance and Navigation (SGN)** course (M.Sc. Space Engineering, Politecnico di Milano, A.Y. 2025-26).

Author: Karim Ahmed Mohamed Doubie.

Course staff: Prof. F. Topputo and P. Di Lizia; TA: A. Morselli and S. Bonaccorsi.

## About

The assignment covers three guidance problems, from trajectory design in the three-body problem to low-thrust optimal control and multiple-shooting transfers in a full n-body model. Each exercise gives the theory, the numerical implementation in MATLAB, and the results with plots and tables.

## Contents

| Exercise | Topic | Main points |
|---|---|---|
| 1 | Lambert problem in the Earth–Moon CR3BP | Single-shooting differential corrector with the State Transition Matrix, numerical continuation in time of flight (0.9 to 0.3), minimum injection maneuver into LEO |
| 2 | Continuous guidance around asteroid Eros | Low-thrust orbit lowering that minimizes dust exposure, nondimensionalization, Pontryagin's Minimum Principle and TPBVP, planar solution, continuation to a 3.5° inclined orbit, n-body vs. Keplerian propagation |
| 3 | Impulsive guidance for the LUMOS mission | Transfer from an Earth orbit to a Distant Prograde Orbit (DPO), n-body propagation, multiple-shooting problem statement with N = 5 nodes, optimization with `fmincon`, alternative full-state formulation |

## Tools used

MATLAB (`ode113`, `fsolve`, `fmincon`), SPICE kernels for ephemerides, LaTeX (report).

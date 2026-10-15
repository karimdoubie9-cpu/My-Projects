# SGN Assignment #2: Spacecraft Guidance and Navigation

Report for the **Spacecraft Guidance and Navigation (SGN)** course (M.Sc. Space Engineering, Politecnico di Milano, A.Y. 2025-26).

Author: Karim Ahmed Mohamed Doubie.

Course staff: Prof. F. Topputo and P. Di Lizia; TA: A. Morselli and S. Bonaccorsi.

## About

The assignment covers three navigation problems: uncertainty propagation, batch orbit determination and sequential filtering. Each exercise gives the theory, the numerical implementation in MATLAB, and the results with plots and tables.

## Contents

| Exercise | Topic | Main points |
|---|---|---|
| 1 | Uncertainty propagation for the INTEGRAL mission | Apsis event detection in Keplerian motion, linearized covariance (LinCov) vs. Unscented Transform (UT), Monte Carlo validation and minimum sample size, uncertainty projected on the orbital plane, time evolution of the 3σ bounds |
| 2 | Batch filters for tracking INTEGRAL | Visibility windows for the TROLL and KOUROU ground stations, simulated azimuth/elevation/range measurements from SGP4 and a TLE, weighted batch least squares (`lsqnonlin`) with Keplerian and J2 dynamics, single-station and two-station solutions, comparison with the SGP4 ground truth over six orbital periods |
| 3 | Sequential filters for the FUTURE CubeSat | Ground track and AI image-processing acquisition windows, simulated GNSS measurements, Unscented Kalman Filter (UKF) with GNSS and with AI-based position measurements, comparison of error and 3σ bounds inside and outside the acquisition windows |

## Tools used

MATLAB (`ode113`, `lsqnonlin`), SPICE kernels and SGP4, LaTeX (report).

# System Identification of a Multirotor UAV: Longitudinal Dynamics

Project for the **Estimation and Learning in Aerospace** course (M.Sc. Space Engineering, Politecnico di Milano). Professor: Marco Lovera.

Authors: Mohamed Ghoraba, Karim Ahmed Mohamed Doubie, Riccardo Cislaghi.

## About

Identification of the longitudinal dynamics of a multirotor UAV from simulated input/output data. The input is the normalised pitching moment (excitation signal) and the outputs are the pitch rate and the longitudinal acceleration. The model is a linearized grey-box model, identified with two methods, tested for sensitivity to sensor noise and finally assessed statistically with a Monte Carlo analysis.

## Contents

| Task | Topic | Main points |
|---|---|---|
| 1 | Grey-box identification | Data prefiltering (4th-order Butterworth band-pass around the PSD peak), non-parametric FRF estimation (direct FFT and Welch method) with coherence check, `greyest` on the FRF, parametric Output-Error (OE) identification in frequency and time domain, asymptotic variance of the parameters, validation with a 3211 sequence |
| 2 | Sensitivity analysis | Identification repeated for different sensor noise levels on pitch rate and longitudinal acceleration, comparison of fit and asymptotic variance, selection of the best model |
| 3 | Monte Carlo analysis | Different noise realizations from the Simulink model, number of iterations from detection guarantees and convergence of the mean, distribution of the parameters, dispersion of poles, zeros and frequency response functions |

## Tools used

MATLAB / Simulink (System Identification Toolbox: `greyest`, OE models).

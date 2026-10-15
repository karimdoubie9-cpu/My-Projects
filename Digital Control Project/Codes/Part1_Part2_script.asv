% =========================================================================
% PROJECT: Digital Controller Design for a Flexible Actuator System
% PART 1: Modeling, Design, and Analytical Verification (Ideal Case)
% =========================================================================
% Objective: To design and verify a discrete-time cascaded control system
% for an aeronautical actuator characterized by flexible transmission.
%
% Note: Part 1 considers an ideal digital control implementation. Sensor, 
% ADC, and actuator non-idealities are introduced in Part 2.
% =========================================================================

clear; clc; close all;

%% STAGE 1: System Identification and Structural Analysis
% Goal: Fully characterize the plant's dynamic behavior and structural 
% properties (Controllability, Observability, and Resonance).

% --> Step 1.1: Define System Parameters and State-Space Model
% -------------------------------------------------------------------------
Ra = 10;                % Armature resistance [Ohm]
La = 0.18e-3;           % Motor armature inductance [H]
Kt = 7.68e-3;           % Torque constant [Nm/A]
Km = 7.68e-3;           % Back-EMF constant [V/(rad/s)]
Kg = 70;                % Total gear ratio
eta_m = 0.85;           % Motor efficiency
eta_g = 0.90;           % Gearbox efficiency
Jeq = 0.0022;           % Equivalent motor inertia [kg*m^2]
Beq = 0.015;            % Equivalent motor damping [Nms/rad]
K1 = 1.5;               % Flexible coupling stiffness [Nm/rad]
J1 = 0.025;             % Load moment of inertia [kg*m^2]
B1 = 0.05;              % Load viscous damping [Nms/rad]

% The state vector is defined as x(t) = [theta_m; omega_m; theta_1; omega_1]
% and the input vector is u(t) = [Va (Control); TL (Disturbance)].

A_full = [   0,                           1,                         0,       0;
          -K1/Jeq,   -Beq/Jeq - eta_m*eta_g*Kg^2*Kt*Km/Ra/Jeq,    K1/Jeq,     0;
             0,                           0,                         0,       1;
           K1/J1,                         0,                      -K1/J1,   -B1/J1];

B_full = [           0,               0;
          eta_m*eta_g*Kg*Kt/Ra/Jeq,   0;
                     0,               0;
                     0,              1/J1];

% Output Selection: Select only speeds for Stage 1.1 analysis
C_full = [0, 1, 0, 0;   % omega_m (Output 1)
          0, 0, 0, 1];  % omega_1 (Output 2)
D_full = zeros(2,2);

sys_full = ss(A_full, B_full, C_full, D_full);

% --> Step 1.2: Structural and Spectral Analysis
% -------------------------------------------------------------------------
% Verification of structural properties ensures a valid control problem.
Co = ctrb(A_full, B_full);
Ob = obsv(A_full, eye(4)); % Verify visibility of all 4 states
isControllable = rank(Co) == size(A_full,1);
isObservable   = rank(Ob) == size(A_full,1);

fprintf('--- Structural Properties ---\n');
fprintf('Controllable: %d, Observable: %d\n', isControllable, isObservable);

% Extract control-relevant Transfer Functions (Input 1: Va)
G_full_tf = tf(sys_full);
G_motor = minreal(G_full_tf(1, 1));      % G_motor = omega_m / Va

      %    18.69 s^2 + 37.39 s + 1122
      % --------------------------------
      % s^3 + 18.87 s^2 + 775.6 s + 2376
G_load_speed = minreal(G_full_tf(2, 1)); % G_load_speed = omega_1 / Va
      %               1122
      % --------------------------------
      % s^3 + 18.87 s^2 + 775.6 s + 2376
 
% SPECTRAL ANALYSIS: Direct Result - Resonance vs Anti-Resonance
plant_poles = pole(G_motor);
motor_zeros = zero(G_motor);
load_zeros  = zero(G_load_speed);

fprintf('\n--- Spectral Analysis ---\n');
disp('System Poles (from G_motor):'); disp(plant_poles);
disp('Motor Speed TF Zeros (from G_motor):'); disp(motor_zeros);

% Custom Pole-Zero Map Plot
figure('Name', 'Pole-Zero Map: Motor vs. Load Dynamics');
hold on;
plot(real(plant_poles), imag(plant_poles), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(motor_zeros), imag(motor_zeros), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(load_zeros), imag(load_zeros), 'ms', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
xlabel('Real Axis [s^{-1}]'); ylabel('Imaginary Axis [s^{-1}]');
title('Pole-Zero Map: Motor vs. Load Dynamics');
legend('System Poles (Resonant)', 'Motor Zeros (Anti-Resonant)', 'Load Zero (No Anti-Resonance)');
axis equal; hold off;

%% STAGE 2: Requirements Translation and Sampling
% Objective: Define performance bandwidth and determine sampling time.

% Rule of Thumb: Bandwidth must be safely below the anti-resonance (wz).
wspeed_loop = 5; % Target Speed Bandwidth [rad/s]

% Rationale for Sampling: Ts derived from speed requirement (fs = 20 * f_loop)
Ts = (2*pi) / (20 * wspeed_loop); 
fs = 1/Ts;
fprintf('\n--- Sampling Parameters ---\n');
fprintf('Chosen Sampling Frequency: %.2f Hz\n', fs);
fprintf('Selected Sampling Time (Ts): %.4f s\n', Ts);

%% STAGE 3: Direct Digital Control Design

% --> Step 3.1: Plant Discretization
Gz_motor = c2d(G_motor, Ts, 'zoh');

% --> Step 3.2: Inner Speed Loop Synthesis (Ragazzini Method)
s = tf('s');
F_s = 1 / (s/wspeed_loop + 1)^2;
F_star = c2d(F_s, Ts, 'zoh');

% The synthesis formula: C(z) = (1/G) * (F / (1-F))
C_speed = minreal( (1/Gz_motor) * (F_star / (1 - F_star)) );

% Inner loop verification
Lspeed=C_speed * Gz_motor;
G_inner_closed = feedback(Lspeed, 1);
      % G_inner_closed = ωₘ(z) / ω_ref(z)
figure('Name', 'Inner Speed Loop Verification');
step(G_inner_closed); title('Step Response of Closed Inner Speed Loop (ω_m / ω_ref)'); grid on;

figure('Name', 'Inner Speed Loop Bode Plot');
margin(minreal(Lspeed)); grid on;

% --> Step 3.3: Outer Position Loop Design
% Goal: Control Output Position (theta_1).
    % The project requirement is to control the OUTPUT shaft position θ1,
    % not the motor shaft position θm.
    %
    % Therefore, the "outer-loop plant" must be:
    %   G_pos_plant(z) = θ1(z) / ω_ref(z)
    % not θm(z)/ω_ref(z).

% G_wm = G_motor;        % ωm/Va
% G_w1 = G_load_speed;   % ω1/Va
G_w1_wm = minreal(G_load_speed / G_motor);% % ω1/ωm
Gz_w1_wm = c2d(G_w1_wm, Ts, 'zoh');
z = tf('z', Ts);
% G_pos_plant(z) = theta_1(z) / w_ref(z)
G_pos_plant = minreal( G_inner_closed * Gz_w1_wm * (Ts/(z - 1)) );
        % G_pos_plant(z) = θ₁(z) / ω_ref(z)=output/input
           % G_pos_plant(z) represents the complete open-loop transfer function for the outer position controller
            % t describes the relationship from the desired angular rate of the motor shaft (ω_ref),
            %  which is the output of the position controller, to the actual angular position of the final load (θ₁).

        % Ts/ (z-1)    : his is a discrete-time integrator that converts the load's angular rate (ω₁) into the load's 
                        % angular position (θ₁). The gain Ts is included to correctly scale the integration in the discrete domain.
        % G_inner_closed=  ωₘ(z) / ω_ref(z): Tthe closed-loop transfer function of the inner speed control loop
                                    % It describes the relationship from the desired angular 
                                    % rate of the motor shaft (ω_ref) to the actual angular 
                                    % rate of the motor shaft (ωₘ). 

% Controller Design: Simple Proportional Gain
omega_bw_pos = 0.9; % initial Bandwidth choice
[mag, ~] = bode(G_pos_plant, omega_bw_pos);
 % Choose Kp so that |Kp * G_pos_plant(j*omega_bw_pos)| = 1
Kp_initial = 1 / squeeze(mag);

% Final Engineering Choice: Manually tuned for optimal time-domain response
Kp = 0.77; 
C_pos = Kp;
L_pos = C_pos * G_pos_plant;

% Search for actual achieved bandwidth
w_vec = logspace(-2, 2, 5000); 
[mag_vec, ~] = bode(L_pos, w_vec);
[~, idx_c] = min(abs(squeeze(mag_vec) - 1));
omega_c_pos_actual = w_vec(idx_c);
fprintf('\nActual Achieved Cut-off Frequency: %.4f rad/s\n', omega_c_pos_actual);

% Verification of the Final Cascaded Position System
G_outer_closed = feedback(L_pos, 1);
        % G_final_closed(z) = θ₁(z) / θ_ref(z)   
S_pos = feedback(1, L_pos); % Sensitivity Function

figure('Name', 'Final Position Loop Verification');
step(G_outer_closed); title('Step Response of Final Closed Position Loop (θ_1 / θ_ref)'); grid on;
S_info = stepinfo(G_outer_closed);
fprintf('\n--- Final Position Loop Metrics --\n');
fprintf('Settling Time: %.2f s\n', S_info.SettlingTime);
fprintf('Overshoot: %.2f %%\n', S_info.Overshoot);

figure('Name', 'Final Position Loop Bode Plot');
margin(L_pos); title('Open-Loop Frequency Response (L_{pos})');

figure('Name', 'Outer Loop Sensitivity Function Bode Plot');
bode(S_pos); title('Outert loop Sensitivity Function (S_{pos}(z))'); grid on;

% Frequency-domain TFs and Bode plots: Va -> speeds
figure('Name','Bode Diagram: Va to Angular Speeds');
bode(G_motor, G_load_speed); grid on;
title('Bode Diagram: $V_a \rightarrow \omega_m,\ \omega_1$','Interpreter','latex');
legend('$\omega_m/V_a$','$\omega_1/V_a$','Interpreter','latex','Location','best');

%% STAGE 4: Full Multi-Loop Disturbance Rejection Analysis
% Objective: Verify rejection of Torque Disturbance TL on theta_1 
% using the total closed-loop MIMO interconnection.


% 0) Discretize the FULL plant (MIMO) once
% Inputs:  u = [V_a ; T_L]
% States:  x = [theta_m; omega_m; theta_1; omega_1]
% Outputs: y = [omega_m; omega_1; theta_m; theta_1]
C_all = [0 1 0 0;   % omega_m
         0 0 0 1;   % omega_1
         1 0 0 0;   % theta_m
         0 0 1 0];  % theta_1
D_all = zeros(4,2);

P = ss(A_full, B_full, C_all, D_all);
Pd = c2d(P, Ts, 'zoh');

% Name I/Os for clean block-diagram interconnection
Pd.InputName  = {'Va','TL'};
Pd.OutputName = {'wm','w1','thm','th1'};

% 1) Ensure controllers have I/O names (must be discrete-time TFs/SS)
% Inner speed controller: input = speed error e_w, output = Va
C_speed.InputName  = {'e_w'};
C_speed.OutputName = {'Va'};

% Outer position controller: input = position error e_th, output = omega_ref
% (Your C_pos is scalar Kp, but keep it as a discrete-time TF or double.)
C_pos_tf = C_pos;                 % if C_pos is numeric (e.g., 0.78), it is OK
C_pos_tf = tf(C_pos_tf, 1, Ts);   % enforce discrete-time TF (safe)
C_pos_tf.InputName  = {'e_th'};
C_pos_tf.OutputName = {'wref'};

% 2) Define summing junctions (cascade structure)
% Inner loop: e_w = wref - wm
Sum_w = sumblk('e_w = wref - wm');

% Outer loop: e_th = thref - th1
Sum_th = sumblk('e_th = thref - th1');

% 3) Build FULL closed-loop interconnection
% Exogenous inputs:  thref (position reference), TL (disturbance torque)
% Outputs of interest: thm, th1 (motor and load positions)
CL_all = connect(Pd, C_speed, C_pos_tf, Sum_w, Sum_th, ...
                 {'thref','TL'}, {'thm','th1','wm','w1'});

% Extract disturbance -> positions (set thref = 0 by selecting input #2)
G_cl_TL_to_thm = minreal(CL_all(1,2));   % thm / TL
G_cl_TL_to_th1 = minreal(CL_all(2,2));   % th1 / TL

% 4) Frequency-domain: Closed-loop Bode from T_L to positions
figure('Name','Closed-Loop Disturbance Rejection (FULL): Bode');
bode(G_cl_TL_to_thm, G_cl_TL_to_th1);
grid on;
title('Closed-Loop Bode: $T_L \rightarrow \theta_m,\theta_1$ (full cascaded loop)','Interpreter','latex');
legend('$\theta_m/T_L$','$\theta_1/T_L$','Interpreter','latex','Location','best');

% 5) Time-domain: Closed-loop step response to a unit step torque disturbance
figure('Name','Closed-Loop Disturbance Rejection (FULL): Step');
step(G_cl_TL_to_thm, G_cl_TL_to_th1);
grid on;
title('Closed-Loop Step: Deviation due to Unit Step $T_L$ (full cascaded loop)','Interpreter','latex');
xlabel('Time (s)');
ylabel('Angular deviation (rad)');
legend('$\theta_m$','$\theta_1$','Interpreter','latex','Location','best');

disp('STAGE 4 Complete: FULL cascaded disturbance rejection computed (TL -> theta_m, theta_1).');

%% STAGE 5: RUN THE SIMULATION OF TASK1
% --- Run the Simulation ---
 out = sim('Part1_2024b.slx');
% --- Create the Plots ---

figure('Name', 'Simulation Results', 'Color', 'w');

% Subplot 1: Step Response (Position Tracking)
subplot(3, 1, 1);
% Plot theta1 against its OWN time vector
plot(out.theta1.Time, out.theta1.Data, 'b-', 'LineWidth', 2);
hold on;
% Plot the step reference against its OWN time vector
plot(out.step.Time, out.step.Data, 'r--', 'LineWidth', 1.5);
hold off;
grid on;
title('System Step Response: Position Tracking');
ylabel('Position \theta_1 (rad)');
xlabel('Time (s)');
legend('Output Position', 'Reference', 'Location', 'SouthEast');

% Subplot 2: Actuator Voltage
subplot(3, 1, 2);
% Plot Va against its OWN time vector
plot(out.Va.Time, out.Va.Data, 'g-', 'LineWidth', 2);
grid on;
title('Actuator Control Input');
ylabel('Voltage V_a (V)');
xlabel('Time (s)');

% Subplot 3: Load Angular Speed
subplot(3, 1, 3);
% Plot w1 against its OWN time vector
plot(out.w1.Time, out.w1.Data, 'm-', 'LineWidth', 2);
% Add a vertical line at t = 0.7854 on the current subplot
yline(0.7854, 'k--', 'LineWidth', 1.5);
grid on;
title('Load Angular Speed');
ylabel('Speed \omega_1 (rad/s)');
xlabel('Time (s)');

% Add an overall title to the figure window
sgtitle('Performance Analysis of the Controlled System');

%% =========================================================================
% PART 2: Hybrid Simulation and Hardware Implementation Verification
% =========================================================================
% Objective: To evaluate the robustness and real-world performance of the
% designed cascaded control law by explicitly modeling the physical
% non-idealities of the data acquisition (ADC) and actuation (DAC) signal chains.
%
% Note: This phase transitions from the ideal mathematical model to a
% high-fidelity hybrid simulation, incorporating effects such as sampling
% jitter, signal quantization, and physical actuator saturation.
% =========================================================================

% parameters needed for simulink

wfilter_position = 20 * omega_c_pos_actual;   % outer loop filter
wfilter_speed = 27;                         % inner loop filter

nbit = 8;   % number of bit for ADC

% --- Run the Simulation ---
 out = sim('Part2_2024b.slx');
% --- Create the Plots ---

figure('Name', 'Simulation Results', 'Color', 'w');

% Subplot 1: Step Response (Position Tracking)
subplot(3, 1, 1);
hold on;
% Plot the step reference against its OWN time vector
plot(out.step.Time, out.step.Data, 'r--', 'LineWidth', 1.5);
% Plot thetam against its OWN time vector
plot(out.thetam.Time, out.thetam.Data, 'LineWidth', 2);
% Plot theta1 against its OWN time vector
plot(out.theta1.Time, out.theta1.Data, 'b-', 'LineWidth', 2);
hold off;
grid on;
title('System Step Response: Position Tracking');
ylabel('Position \theta (deg)');
xlabel('Time (s)');
legend('Reference', '\theta_m', '\theta_1', 'Location', 'SouthEast');

% Subplot 2: Actuator Voltage
subplot(3, 1, 2);
% Plot Va against its OWN time vector
plot(out.Va.Time, out.Va.Data, 'g-', 'LineWidth', 2);
grid on;
title('Actuator Control Input');
ylabel('Voltage V_a (V)');
xlabel('Time (s)');

% Subplot 3: Load Angular Speed
subplot(3, 1, 3);
hold on
% Plot wm against its OWN time vector
plot(out.wm.Time, out.wm.Data, 'LineWidth', 2);
% Plot w1 against its OWN time vector
plot(out.w1.Time, out.w1.Data, 'm-', 'LineWidth', 2);
% Add a vertical line at t = 0.7854 on the current subplot
yline(45, 'k--', 'LineWidth', 1.5);
grid on;
title('Load Angular Speed');
ylabel('Speed \omega (deg/s)');
xlabel('Time (s)');
legend('\omega_m', '\omega_1', '45 deg/s')
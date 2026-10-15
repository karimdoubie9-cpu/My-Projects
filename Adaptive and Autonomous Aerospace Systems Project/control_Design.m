viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

s = tf('s');
Gm = 20/(s+20); 
Gd = exp(-0.01*s);

J = 0.01;

UAV.D_fv       = diag([0.26 0.26 0]);                % linear velocity damping
UAV.D_tauomega = diag([3.6e-4 3.6e-4 1.20e-2]);      % angular velocity damping

D = UAV.D_tauomega(1,1);
G = Gm*Gd/(J*s^2 + D*s);

%% use sisotool with G and get KR and Kw
KR = 2;
Kw = 0.4;

C2 = Kw*s; 
C1=KR;
Hatt = 1 + C2/C1;

Tatt = feedback(G*C1,Hatt); % closed loop of attitude loop

dcgain(Tatt)

%% position loop


Dv = UAV.D_fv(1,1);

Ptrans = 1/(1.2*s^2 + Dv*s);

Pxy = minreal(Tatt * Ptrans);
%%
controlSystemDesigner(Pxy)

%%
s=tf('s');
Gm = 20/(s+20);
Gd = exp(-0.01*s);

Pz =Gm*Gd* 1/(1.2*s^2);
controlSystemDesigner(Pz)
%%
gains = position_gains_from_CSD(C1,C2)

%% Bode: actuator G_m*G_d, with vs. without the delay
%  G_m = k_m/(s+k_m) is the motor lag, G_d = exp(-tau*s) the 10 ms delay.
%  A pure delay has |G_d| = 1 at every frequency, so the two magnitude curves
%  coincide: the delay's whole cost is extra phase lag, -omega*tau, which grows
%  linearly with frequency and is what erodes the attitude loop's phase margin.
s  = tf('s');
Gm = 20/(s+20);
Gd = exp(-0.01*s);
w  = logspace(-1, 3, 800);                         % [rad/s]

[m0, p0] = bode(Gm,    w);
[m1, p1] = bode(Gm*Gd, w);
m0 = squeeze(m0);  p0 = rad2deg(unwrap(deg2rad(squeeze(p0))));
m1 = squeeze(m1);  p1 = rad2deg(unwrap(deg2rad(squeeze(p1))));

figure('Name','Bode G_m G_d with vs without delay','Color','w','Position',[140 140 760 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
title(tl, 'G_m G_d: motor lag with and without the 10 ms delay', 'FontName', get(groot, 'defaultAxesFontName'));

ax1 = nexttile;
semilogx(w, 20*log10(m0), 'LineWidth', 2.1); hold on
semilogx(w, 20*log10(m1), '--', 'LineWidth', 2.1); grid on
ylabel('|G|  [dB]');
legend('G_m  (no delay)', 'G_m G_d  (10 ms delay)', 'Location', 'southwest');

ax2 = nexttile;
semilogx(w, p0, 'LineWidth', 2.1); hold on
semilogx(w, p1, '--', 'LineWidth', 2.1); grid on
xline(20, ':', 'motor pole, 20 rad/s', 'HandleVisibility', 'off');
ylabel('\angle G  [deg]');  xlabel('\omega  [rad/s]');

linkaxes([ax1 ax2], 'x');  xlim([w(1) w(end)]);

fprintf('Extra phase lag from the delay: %.1f deg at 20 rad/s, %.1f deg at 100 rad/s\n', ...
        -rad2deg(20*0.01), -rad2deg(100*0.01));

%% Closed-loop attitude T_att: Bode and 1 rad step response
%  Uses Tatt from the first section, so it reflects whatever KR and Kw are set
%  there. The model works in radians; the step response is plotted in degrees
%  (1 rad = 57.3 deg) because attitude reads more naturally that way.
if ~exist('Tatt', 'var')
    error('control_Design:noTatt', 'Run the first section first: it builds Tatt.');
end

w = logspace(-1, 3, 800);                         % [rad/s]
[mT, pT] = bode(Tatt, w);
mT_dB = 20*log10(squeeze(mT));
pT    = rad2deg(unwrap(deg2rad(squeeze(pT))));

% -3 dB bandwidth (first crossing, the standard definition) and the resonance.
% With these gains the magnitude dips below -3 dB and then climbs back to about
% 0 dB at a lightly damped mode before rolling off, so the resonance is looked
% for AFTER the first crossing; it is the same mode that shows up as the ripple
% in the step response.
dcT   = dcgain(Tatt);
i_bw  = find(mT_dB < 20*log10(dcT) - 3, 1);
w_bw  = NaN;  if ~isempty(i_bw), w_bw = w(i_bw); end
if isfinite(w_bw)
    beyond = find(w > w_bw);
    [Mr_dB, k] = max(mT_dB(beyond));  i_r = beyond(k);
else
    [Mr_dB, i_r] = max(mT_dB);
end

figure('Name','Closed-loop attitude Bode','Color','w','Position',[140 140 760 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
title(tl, sprintf('Closed-loop attitude T_{att}   (K_R = %g, K_\\omega = %g)', KR, Kw), 'FontName', get(groot, 'defaultAxesFontName'));

ax1 = nexttile;
semilogx(w, mT_dB, 'LineWidth', 2.1); grid on; hold on
yline(20*log10(dcT) - 3, ':', '-3 dB', 'HandleVisibility', 'off');
if isfinite(w_bw)
    xline(w_bw, '--', sprintf('\\omega_{BW} = %.1f rad/s', w_bw), ...
          'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end
plot(w(i_r), Mr_dB, 'o', 'MarkerSize', 7, 'LineWidth', 1.8, 'HandleVisibility', 'off');
text(w(i_r)*1.15, Mr_dB - 6, sprintf('resonance %.1f dB at %.1f rad/s', Mr_dB, w(i_r)));
ylabel('|T_{att}|  [dB]');

ax2 = nexttile;
semilogx(w, pT, 'LineWidth', 2.1); grid on
if isfinite(w_bw), xline(w_bw, '--', 'HandleVisibility', 'off'); end
ylabel('\angle T_{att}  [deg]');  xlabel('\omega  [rad/s]');

linkaxes([ax1 ax2], 'x');  xlim([w(1) w(end)]);

% ---- 1 rad step, shown in degrees -----------------------------------
step_rad = 1;                                     % [rad] commanded attitude step
[y, t]   = step(step_rad * Tatt);
S        = stepinfo(y, t, step_rad * dcT);        % 2% settling band by default

figure('Name','Closed-loop attitude step 1 rad','Color','w','Position',[160 160 760 420]);
plot(t, rad2deg(y), 'LineWidth', 2.1); grid on; hold on
yline(rad2deg(step_rad), '--', sprintf('command: %.1f deg (%g rad)', rad2deg(step_rad), step_rad), ...
      'HandleVisibility', 'off');
xlabel('t  [s]');  ylabel('attitude  [deg]');
title('Closed-loop attitude response to a 1 rad step');
text(0.62, 0.30, sprintf('rise time  %.3f s\novershoot  %.1f %%\nsettling (2%%)  %.3f s', ...
     S.RiseTime, S.Overshoot, S.SettlingTime), 'Units', 'normalized', ...
     'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);

fprintf('\nClosed-loop attitude T_att: DC gain %.3f | -3 dB bandwidth %.2f rad/s | resonance %.2f dB at %.2f rad/s\n', ...
        dcT, w_bw, Mr_dB, w(i_r));
fprintf('1 rad step: rise %.3f s | overshoot %.1f %% | settling (2%%) %.3f s | peak %.1f deg\n', ...
        S.RiseTime, S.Overshoot, S.SettlingTime, rad2deg(S.Peak));

%% Closed-loop xy position loop T_xy: Bode and 1 m step response
%  Same 2-DOF structure as the attitude loop, one level up:
%
%     u = C1*(r - y) - C2*y,    C1 = K_p + K_i/s,    C2 = K_v*s
%
%  on the design plant Pxy = Tatt*Ptrans (the closed attitude loop acting as
%  the force actuator). From reference to position:
%
%     T_xy = Pxy*C1 / (1 + Pxy*(C1 + C2))
%
%  C2 = K_v*s is improper, and MATLAB cannot combine an improper model with
%  the delay inside Tatt. So the plant is given two proper outputs, position
%  AND velocity (velocity = s*position), and K_v acts on the velocity output.
%  The controller is realised as one system with a single integrator, so the
%  closed loop has no spurious pole-zero cancellation at s = 0.
%  Gains are the ones the Simulink model uses (design.position).
if ~exist('Tatt', 'var')
    error('control_Design:noTatt', 'Run the first section first: it builds Tatt.');
end
s   = tf('s');
m   = 1.2;                                        % [kg] design mass, no package
Dv  = UAV.D_fv(1,1);
pos = design.position(UAV);
Kp  = pos.K_p;  Kv = pos.K_v;  Ki = pos.K_i;

Pv  = Tatt * (1/(m*s + Dv));                      % force -> velocity
Pyv = [1/s; tf(1)] * Pv;                          % force -> [position; velocity]
Pyv.InputName  = 'u';
Pyv.OutputName = {'y','v'};

C1p = Kp + Ki/s;
Kc  = minreal(ss([C1p, -C1p, -Kv]), [], false);   % u = C1*(r - y) - K_v*v
Kc.InputName  = {'r','y','v'};
Kc.OutputName = 'u';

Txy = connect(Kc, Pyv, 'r', 'y');                 % reference -> position

w = logspace(-2, 3, 800);                         % [rad/s]
[mX, pX] = bode(Txy, w);
mX_dB = 20*log10(squeeze(mX));
pX    = rad2deg(unwrap(deg2rad(squeeze(pX))));

dcX    = dcgain(Txy);
i_bwX  = find(mX_dB < 20*log10(dcX) - 3, 1);
w_bwX  = NaN;  if ~isempty(i_bwX), w_bwX = w(i_bwX); end
[MrX_dB, i_rX] = max(mX_dB);                      % resonant peak, if any

% Stability margins come from an OPEN loop, not from T_xy. As in Control
% System Designer (LoopTransfer_C1), the loop is broken at C1's output with the
% state feedback C2 (velocity) still closed:
%
%     L_C1 = C1 * Pxy/(1 + Pxy*C2)
%
% The inner C2 loop is closed on the two-output plant (K_v acting on the
% velocity output), which keeps every piece proper despite the delay in Tatt.
Pin  = feedback(Pyv, [0 Kv]);                     % u -> [y; v], with u = w - K_v*v
L_C1 = C1p * Pin('y', 'u');
[~, PMxy, ~, wc_xy] = margin(L_C1);               % [Gm, Pm, Wcg, Wcp]
fprintf('\nxy position loop (LoopTransfer_C1): phase margin %.1f deg at crossover %.3f rad/s\n', ...
        PMxy, wc_xy);

figure('Name','Closed-loop xy position Bode','Color','w','Position',[140 140 760 560]);
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
title(tl, 'Closed-loop xy position T_{xy}', 'FontName', get(groot, 'defaultAxesFontName'));

ax1 = nexttile;
semilogx(w, mX_dB, 'LineWidth', 2.1); grid on
ylabel('|T_{xy}|  [dB]');

ax2 = nexttile;
semilogx(w, pX, 'LineWidth', 2.1); grid on
ylabel('\angle T_{xy}  [deg]');  xlabel('\omega  [rad/s]');

linkaxes([ax1 ax2], 'x');  xlim([w(1) w(end)]);

% ---- 1 m step -------------------------------------------------------
step_m = 1;                                       % [m] commanded position step
[yX, tX] = step(step_m * Txy);
SX       = stepinfo(yX, tX, step_m * dcX);        % 2% settling band by default

figure('Name','Closed-loop xy position step 1 m','Color','w','Position',[160 160 760 420]);
plot(tX, yX, 'LineWidth', 2.1); grid on; hold on
yline(step_m, '--', sprintf('command: %g m', step_m), 'HandleVisibility', 'off');
xlabel('t  [s]');  ylabel('position  [m]');
title('Closed-loop xy position response to a 1 m step');
text(0.62, 0.30, sprintf('rise time  %.2f s\novershoot  %.1f %%\nsettling (2%%)  %.2f s', ...
     SX.RiseTime, SX.Overshoot, SX.SettlingTime), 'Units', 'normalized', ...
     'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);

fprintf('\nClosed-loop xy position T_xy: DC gain %.3f | -3 dB bandwidth %.3f rad/s | peak %.2f dB at %.3f rad/s\n', ...
        dcX, w_bwX, MrX_dB, w(i_rX));
fprintf('1 m step: rise %.2f s | overshoot %.1f %% | settling (2%%) %.2f s | peak %.3f m\n', ...
        SX.RiseTime, SX.Overshoot, SX.SettlingTime, SX.Peak);

%%
function gains = position_gains_from_CSD(C1,C2)
%POSITION_GAINS_FROM_CSD
% Extract project position gains from the 2-DOF architecture:
%
%   u = C1*(r-y) - C2*y
%
% expected:
%   C1 = Kp + Ki/s
%   C2 = Kv*s
%
% Returns:
%   gains.K_p
%   gains.K_v
%   gains.K_i

    s = tf('s');

    % Clean numerical pole/zero cancellations
    C1 = minreal(tf(C1));
    C2 = minreal(tf(C2));

    % ----- C1 = Kp + Ki/s -----
    % Ki = lim(s->0) s*C1(s)
    Ki = dcgain(minreal(s*C1));

    % Kp = high-frequency gain of C1
    Kp = evalfr(C1,1e8);

    % ----- C2 = Kv*s -----
    % Kv = C2(s)/s
    Kv = dcgain(minreal(C2/s));

    gains.K_p = real(Kp);
    gains.K_v = real(Kv);
    gains.K_i = real(Ki);

    fprintf('\nPosition controller gains\n');
    fprintf('K_p = %.6g\n',gains.K_p);
    fprintf('K_v = %.6g\n',gains.K_v);
    fprintf('K_i = %.6g\n',gains.K_i);

end
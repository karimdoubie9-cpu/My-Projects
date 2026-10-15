function att = attitude(UAV, mode, opts)
%ATTITUDE  Design and analyse the geometric attitude-loop gains.
%
%   att = design.attitude(UAV)
%   att = design.attitude(UAV,'designed')
%   att = design.attitude(UAV,'as_flown')
%   att = design.attitude(UAV,'legacy_wn')
%
%   The default 'designed' mode is a frequency-domain design based on a
%   requested open-loop crossover frequency and phase margin.  Unlike the old
%   wn/zeta placement, it explicitly includes:
%       - rigid-body rotational inertia;
%       - optional aerodynamic rotational damping;
%       - the first-order actuator pole;
%       - a pure sensing/computation delay used for tuning.
%
%   This matches the exam requirement to tune on a simplified/linearised model
%   while accounting for actuator/sensor bandwidth and delay.  Saturation is
%   an amplitude constraint and must still be verified after control allocation
%   in the nonlinear model; this function reports an optional torque-demand
%   check when limits are supplied.
%
%   CONTROL LAW IMPLEMENTED IN Final.slx
%   ------------------------------------
%   The geometric error block evaluates
%
%       e_R = k_R1*(Re'*e1 x e1) + k_R2*(Re'*e2 x e2) ...
%           + k_R3*(Re'*e3 x e3),      Re = R_d' * R
%
%   and the rate term uses ctrl.Komega.  For a small attitude error
%   Re ~= I + skew(delta),
%
%       e_R ~= diag(k_R2+k_R3, k_R1+k_R3, k_R1+k_R2) * delta.
%
%   Therefore the effective per-axis proportional gain is
%
%       Keff_R = [k_R2+k_R3, k_R1+k_R3, k_R1+k_R2],
%
%   not k_R itself.  This function always designs in terms of Keff_R and then
%   inverts that map to obtain the three geometric k_R coefficients actually
%   required by the Simulink block.
%
%   SMALL-SIGNAL AXIS MODEL USED FOR TUNING
%   ----------------------------------------
%       P_i(s) = [wm/(s+wm)] * 1/(J_i*s^2 + D_i*s) * exp(-tau_d*s)
%       C_i(s) = K_R,i + K_omega,i*s
%
%   At the requested crossover wc_i we impose
%
%       |C_i(j wc_i) P_i(j wc_i)| = 1
%       angle(C_i P_i) = -180 deg + PM_i.
%
%   Since C_i(jw) = K_R,i + j*w*K_omega,i, these two equations give the
%   gains directly.  This is preferable to treating wn as if it were the
%   actual closed-loop bandwidth once actuator dynamics and delay are present.
%
%   DEFAULT DESIGN
%   --------------
%       roll/pitch crossover : 6 rad/s
%       yaw crossover        : 3 rad/s
%       target phase margin  : 55 deg
%       sensor/compute delay : 10 ms (explicit design assumption)
%       actuator pole        : read from UAV if possible, otherwise 20 rad/s
%
%   For the nominal project inertia J = diag([0.010 0.010 0.018]) kg m^2,
%   wm = 20 rad/s, negligible D and tau_d = 0.01 s, this gives approximately
%
%       Keff_R    = [0.0964  0.0964  0.0686]
%       Keff_w    = [0.0605  0.0605  0.0496]
%       k_R       = [0.0343  0.0343  0.0621].
%
%   Name-value examples:
%       att = design.attitude(UAV,'designed', ...
%               wc_rollpitch=7, wc_yaw=3, ...
%               phase_margin_rollpitch=55, phase_margin_yaw=55, ...
%               sensor_delay=0.010);
%
%       % Optional amplitude sanity check (NOT a replacement for allocation):
%       att = design.attitude(UAV,'designed', ...
%               tau_limit=[1.2;1.2;1.0], ...
%               design_angle_deg=[20;20;25], ...
%               design_rate_deg_s=[120;120;90]);
%
%   IMPORTANT
%   ---------
%   This file tunes the gains of the controller already implemented in
%   Final.slx.  It does not add desired-angular-rate / angular-acceleration
%   feedforward terms.  Such terms are required for exact tracking of a general
%   time-varying attitude trajectory; if Final.slx neglects them, the resulting
%   tracking error must be assessed in Task 2.3.

arguments
    UAV  struct
    mode (1,:) char {mustBeMember(mode,{'designed','as_flown','legacy_wn'})} = 'designed'

    % Frequency-domain design targets
    opts.wc_rollpitch              (1,1) double {mustBePositive}    = 6      % [rad/s]
    opts.wc_yaw                    (1,1) double {mustBePositive}    = 3      % [rad/s]
    opts.phase_margin_rollpitch    (1,1) double {mustBePositive}    = 55     % [deg]
    opts.phase_margin_yaw          (1,1) double {mustBePositive}    = 55     % [deg]
    opts.sensor_delay              (1,1) double {mustBeNonnegative} = 0.010  % [s], design assumption

    % Set NaN to auto-detect from UAV; falls back to 20 rad/s (= 0.05 s).
    opts.motor_pole                (1,1) double = NaN                        % [rad/s]

    % Rotational aerodynamic damping diag(D_tauomega), [N m s/rad].
    % NaNs request auto-detection from the UAV struct; zero is used if absent.
    opts.Domega                    (3,1) double = [NaN;NaN;NaN]

    % Legacy second-order design, retained only for comparison.
    opts.legacy_zeta               (1,1) double {mustBePositive} = 0.8
    opts.legacy_wn_rollpitch       (1,1) double {mustBePositive} = 8          % [rad/s]
    opts.legacy_wn_yaw             (1,1) double {mustBePositive} = 3          % [rad/s]

    % Optional amplitude check.  Leave tau_limit as NaN unless you have an
    % allocation-aware torque limit at the operating collective thrust.
    opts.tau_limit                 (3,1) double = [NaN;NaN;NaN]               % [N m]
    opts.design_angle_deg          (3,1) double {mustBeNonnegative} = [20;20;25]
    opts.design_rate_deg_s         (3,1) double {mustBeNonnegative} = [120;120;90]
end

att.mode = mode;
att.design_method = 'frequency-domain crossover / phase-margin placement';

% -------------------------------------------------------------------------
% Plant data used by the simplified tuning model
% -------------------------------------------------------------------------
Jd = diag(UAV.J)';
if numel(Jd) ~= 3 || any(~isfinite(Jd)) || any(Jd <= 0)
    error('design:attitude:BadInertia', ...
        'UAV.J must contain three positive principal inertias.');
end

[wm, wm_source] = resolveMotorPole(UAV, opts.motor_pole);
[Domega, D_source] = resolveRotDamping(UAV, opts.Domega);

att.motor_pole = wm;
att.motor_tau  = 1/wm;
att.motor_pole_source = wm_source;
att.sensor_delay = opts.sensor_delay;
att.Domega = Domega;
att.Domega_source = D_source;
att.Jdiag = Jd;

% Requested targets (stored for plots/reporting regardless of mode).
att.wc_target = [opts.wc_rollpitch opts.wc_rollpitch opts.wc_yaw];
att.phase_margin_target = [opts.phase_margin_rollpitch ...
                           opts.phase_margin_rollpitch ...
                           opts.phase_margin_yaw];

% -------------------------------------------------------------------------
% Gain calculation
% -------------------------------------------------------------------------
switch mode
    case 'as_flown'
        % Old gains kept only as a reproducibility/reference option.
        att.k_R    = [2 2 2];
        att.Komega = eye(3);

    case 'legacy_wn'
        % Old textbook second-order placement.  This intentionally ignores
        % actuator dynamics and delay and is therefore NOT the recommended
        % final exam design.
        wn = [opts.legacy_wn_rollpitch ...
              opts.legacy_wn_rollpitch ...
              opts.legacy_wn_yaw];

        K_R_eff = Jd .* wn.^2;
        K_w_eff = 2*opts.legacy_zeta .* wn .* Jd;
        
        att.k_R    = effectiveToGeometric(K_R_eff);
        att.Komega = diag(K_w_eff);

    case 'designed'
        wc = att.wc_target;
        PM = att.phase_margin_target;

        K_R_eff = zeros(1,3);
        K_w_eff = zeros(1,3);
        axis_design = [];

        for i = 1:3
            [K_R_eff(i), K_w_eff(i), info_i] = ...
                designAxisAtCrossover(Jd(i), Domega(i), wm, ...
                                      opts.sensor_delay, wc(i), PM(i));
            if i == 1
                axis_design = repmat(info_i,1,3);
            else
                axis_design(i) = info_i;
            end
        end
        K_R_eff=2*ones(3,1);
        K_w_eff=0.4*ones(3,1);
        % K_R_eff=1*ones(3,1);
        % K_w_eff=0.2 *ones(3,1);
        
        att.k_R    = effectiveToGeometric(K_R_eff);
        att.Komega = diag(K_w_eff);
        att.axis_design = axis_design;
end

% -------------------------------------------------------------------------
% Linear-equivalent gains ACTUALLY seen by the small-error loop
% -------------------------------------------------------------------------
att.Keff_R = [att.k_R(2)+att.k_R(3), ...
              att.k_R(1)+att.k_R(3), ...
              att.k_R(1)+att.k_R(2)];
att.Keff_w = diag(att.Komega)';

% Backward-compatible names used by existing plotting scripts.
att.K_R_target     = att.Keff_R;
att.K_omega_target = att.Keff_w;

% Equivalent 2nd-order quantities if actuator/delay are ignored.  These are
% diagnostics only; they are NOT the bandwidth/crossover of the real design.
att.wn_equiv = sqrt(att.Keff_R ./ Jd);
att.zeta_equiv = (Domega + att.Keff_w) ./ ...
                 (2*sqrt(Jd .* att.Keff_R));
att.wn_rollpitch = mean(att.wn_equiv(1:2));
att.wn_yaw       = att.wn_equiv(3);
att.zeta         = mean(att.zeta_equiv(1:2));

% -------------------------------------------------------------------------
% Frequency-domain verification on the SAME model used for design
% -------------------------------------------------------------------------
att.freq = [];
for i = 1:3
    info_i = analyseAxis(Jd(i), Domega(i), att.Keff_R(i), ...
                         att.Keff_w(i), wm, opts.sensor_delay, ...
                         att.wc_target(i));
    if i == 1
        att.freq = repmat(info_i,1,3);
    else
        att.freq(i) = info_i;
    end
end

att.wc_actual = [att.freq.wc];
att.phase_margin = [att.freq.phase_margin_deg];
att.bandwidth_3db = [att.freq.bandwidth_3db];
att.gain_margin_db = [att.freq.gain_margin_db];
att.sensitivity_peak = [att.freq.Ms];
att.complementary_sensitivity_peak = [att.freq.Mt];
att.additional_delay_margin = [att.freq.additional_delay_margin_s];

% A simple design-health summary.
att.actuator_separation = wm ./ att.wc_actual;
att.design_ok = all(att.phase_margin >= 45) && ...
                all(att.actuator_separation >= 2.5);

% -------------------------------------------------------------------------
% Optional amplitude sanity check
% -------------------------------------------------------------------------
% The true feasibility constraint is on individual rotor thrusts after
% allocation.  This only checks whether a representative small-error PD demand
% is below a user-supplied torque envelope.
ang = deg2rad(opts.design_angle_deg(:))';
rate = deg2rad(opts.design_rate_deg_s(:))';
att.design_tau_demand = att.Keff_R .* ang + att.Keff_w .* rate;
att.tau_limit = opts.tau_limit(:)';
att.torque_check_available = all(isfinite(att.tau_limit));
if att.torque_check_available
    att.torque_utilisation = att.design_tau_demand ./ att.tau_limit;
    att.torque_check_ok = all(att.torque_utilisation <= 1);
else
    att.torque_utilisation = [NaN NaN NaN];
    att.torque_check_ok = NaN;
end

end

% =========================================================================
% Local helpers
% =========================================================================
function [KR, Kw, info] = designAxisAtCrossover(J, D, wm, tau, wc, PMdeg)
%DESIGNAXISATCROSSOVER Solve |L|=1 and phase(L)=-180+PM at wc.

% Rigid-body plant 1/(J s^2 + D s)
rb_den_mag = hypot(J*wc^2, D*wc);
rb_mag = 1/rb_den_mag;
rb_phase = -atan2(D*wc, -J*wc^2);      % rad, roughly -pi

% First-order actuator wm/(s+wm)
act_mag   = 1/sqrt(1 + (wc/wm)^2);
act_phase = -atan(wc/wm);

% Sensor/computation delay exp(-tau*s)
delay_phase = -wc*tau;

Pmag   = rb_mag * act_mag;
Pphase = rb_phase + act_phase + delay_phase;

required_loop_phase = -pi + deg2rad(PMdeg);
controller_phase = required_loop_phase - Pphase;

% The implemented PD/rate-feedback controller can only supply 0...+90 deg.
if ~(controller_phase > 0 && controller_phase < pi/2)
    error('design:attitude:InfeasibleTarget', ...
        ['Requested wc = %.3g rad/s and PM = %.1f deg require %.1f deg ' ...
         'of controller phase lead.  The implemented C(s)=KR+Kw*s can only ' ...
         'provide between 0 and 90 deg.  Reduce wc/PM or reduce the assumed delay.'], ...
         wc, PMdeg, rad2deg(controller_phase));
end

Cmag = 1/Pmag;
KR = Cmag*cos(controller_phase);
Kw = Cmag*sin(controller_phase)/wc;

if KR <= 0 || Kw <= 0
    error('design:attitude:NonPositiveGain', ...
        'Crossover design produced non-positive attitude gains.');
end

info.wc = wc;
info.phase_margin_target_deg = PMdeg;
info.rigid_body_phase_deg = rad2deg(rb_phase);
info.actuator_phase_deg = rad2deg(act_phase);
info.delay_phase_deg = rad2deg(delay_phase);
info.plant_phase_deg = rad2deg(Pphase);
info.controller_phase_deg = rad2deg(controller_phase);
info.actuator_mag_db = 20*log10(act_mag);
info.KR = KR;
info.Kw = Kw;
end

function kR = effectiveToGeometric(Keff)
%EFFECTIVETOGEOMETRIC Invert Keff=[k2+k3,k1+k3,k1+k2].
A = ones(3) - eye(3);
kR = (A \ Keff(:))';

if any(kR <= 0)
    error('design:attitude:NegativeGeometricGain', ...
        ['Requested effective stiffness cannot be represented with all-positive ' ...
         'geometric k_R coefficients.  Effective gains were [%g %g %g].'], ...
        Keff(1),Keff(2),Keff(3));
end
end

function out = analyseAxis(J, D, KR, Kw, wm, tau, wc_hint)
%ANALYSEAXIS Numerical frequency check without requiring Control System Toolbox.

wmin = max(1e-3, wc_hint/200);
wmax = max([20*wc_hint, 10*wm, 200]);
w = logspace(log10(wmin), log10(wmax), 8000);
s = 1i*w;

P = (wm ./ (s + wm)) .* exp(-s*tau) ./ (J*s.^2 + D*s);
C = KR + Kw*s;
L = C .* P;
T = L ./ (1 + L);
S = 1 ./ (1 + L);

magdb = 20*log10(abs(L));
phase = unwrap(angle(L));
phase_deg = rad2deg(phase);

% Gain crossover nearest the intended crossover.
idx_cross = find(magdb(1:end-1).*magdb(2:end) <= 0);
if isempty(idx_cross)
    wc = NaN;
    pm = NaN;
else
    cand = zeros(size(idx_cross));
    for k = 1:numel(idx_cross)
        ii = idx_cross(k);
        cand(k) = interpLogX(w(ii),w(ii+1),magdb(ii),magdb(ii+1),0);
    end
    [~,kbest] = min(abs(log(cand/wc_hint)));
    wc = cand(kbest);
    ph_wc = interp1(log(w), phase_deg, log(wc), 'linear');
    pm = 180 + ph_wc;
end

% Closed-loop -3 dB bandwidth (T(0) ~= 1 for this type-2 attitude loop).
Tdb = 20*log10(abs(T));
idx_bw = find(Tdb <= -3,1,'first');
if isempty(idx_bw) || idx_bw == 1
    bw = NaN;
else
    bw = interpLogX(w(idx_bw-1),w(idx_bw),Tdb(idx_bw-1),Tdb(idx_bw),-3);
end

% High-frequency -180 deg crossing after gain crossover -> gain margin.
gm_db = Inf;
wpc = NaN;
if isfinite(wc)
    idx0 = find(w >= wc,1,'first');
    ph_rel = phase_deg + 180;
    idx_pc_rel = find(ph_rel(idx0:end-1) >= 0 & ph_rel(idx0+1:end) < 0,1,'first');
    if ~isempty(idx_pc_rel)
        ii = idx0 + idx_pc_rel - 1;
        wpc = interpLogX(w(ii),w(ii+1),phase_deg(ii),phase_deg(ii+1),-180);
        mag_pc = interp1(log(w),magdb,log(wpc),'linear');
        gm_db = -mag_pc;
    end
end

out.wc = wc;
out.phase_margin_deg = pm;
out.bandwidth_3db = bw;
out.gain_margin_db = gm_db;
out.phase_crossover = wpc;
out.Ms = max(abs(S));
out.Mt = max(abs(T));
if isfinite(wc) && isfinite(pm)
    out.additional_delay_margin_s = deg2rad(pm)/wc;
else
    out.additional_delay_margin_s = NaN;
end
end

function x = interpLogX(x1,x2,y1,y2,yt)
%INTERPLOGX Linear interpolation in y against log(x).
if y2 == y1
    x = sqrt(x1*x2);
    return
end
a = (yt-y1)/(y2-y1);
x = exp(log(x1) + a*(log(x2)-log(x1)));
end

function [wm, source] = resolveMotorPole(UAV, user_wm)
%RESOLVEMOTORPOLE Read actuator pole from common project field names.
if isfinite(user_wm)
    if user_wm <= 0
        error('design:attitude:BadMotorPole','motor_pole must be positive.');
    end
    wm = user_wm;
    source = 'name-value motor_pole';
    return
end

if isfield(UAV,'k_m') && isscalar(UAV.k_m) && isfinite(UAV.k_m) && UAV.k_m > 0
    wm = UAV.k_m;
    source = 'UAV.k_m';
elseif isfield(UAV,'km') && isscalar(UAV.km) && isfinite(UAV.km) && UAV.km > 0
    wm = UAV.km;
    source = 'UAV.km';
elseif isfield(UAV,'tau_m') && isscalar(UAV.tau_m) && isfinite(UAV.tau_m) && UAV.tau_m > 0
    wm = 1/UAV.tau_m;
    source = '1/UAV.tau_m';
elseif isfield(UAV,'motor_tau') && isscalar(UAV.motor_tau) && isfinite(UAV.motor_tau) && UAV.motor_tau > 0
    wm = 1/UAV.motor_tau;
    source = '1/UAV.motor_tau';
else
    wm = 20;
    source = 'fallback 20 rad/s (tau_m = 0.05 s)';
end
end

function [D, source] = resolveRotDamping(UAV, user_D)
%RESOLVEROTDAMPING Resolve diagonal rotational aerodynamic damping.
if all(isfinite(user_D))
    D = user_D(:)';
    if any(D < 0)
        error('design:attitude:BadDamping','Domega must be nonnegative.');
    end
    source = 'name-value Domega';
    return
end

names = {'D_tauomega','D_tau_omega','Dtauomega','D_tau','Dtau','Domega','D_omega','Drot','D_w'};
for k = 1:numel(names)
    if isfield(UAV,names{k})
        val = UAV.(names{k});
        if isequal(size(val),[3 3])
            val = diag(val);
        end
        if numel(val) == 3 && all(isfinite(val(:))) && all(val(:) >= 0)
            D = val(:)';
            source = ['UAV.' names{k}];
            return
        end
    end
end

D = [0 0 0];
source = 'zero (no rotational damping field detected)';
end

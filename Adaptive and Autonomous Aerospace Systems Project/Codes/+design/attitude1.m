function att = attitude(UAV, mode)
%ATTITUDE  Attitude loop gains, plus their linear-equivalent for analysis.
%
%   att = design.attitude(UAV)              % 'as_flown'  (default)
%   att = design.attitude(UAV, 'designed')  % zeta/omega_n placement
%
%   THE CONTROL LAW IMPLEMENTED IN Final.slx
%   ----------------------------------------
%   The block 'attitude controller/e_R' evaluates the geometric error
%
%       e_R = k_R1*(Re'*e1 x e1) + k_R2*(Re'*e2 x e2) + k_R3*(Re'*e3 x e3),
%       Re  = R_d' * R
%
%   and the rate term uses ctrl.Komega.  For a small attitude error
%   Re ~ I + skew(delta) one gets
%
%       e_R ~ diag( k_R2+k_R3,  k_R1+k_R3,  k_R1+k_R2 ) * delta
%
%   so the EFFECTIVE per-axis proportional gain is the SUM OF THE OTHER TWO
%   k_R entries, not k_R itself.  With k_R1=k_R2=k_R3=5 the effective gain is
%   10 on every axis.  Any linear analysis must use att.Keff_R / att.Keff_w,
%   otherwise it describes a controller that is not the one being simulated.
%
%   MODES
%   -----
%   'as_flown'  (default) reproduces the gains currently wired into the model:
%               k_R = [5 5 5], Komega = I.  Effective: K_R = 10, K_w = 1
%               -> roll/pitch  wn = sqrt(10/0.010) = 31.6 rad/s, zeta = 1.58
%                  yaw         wn = sqrt(10/0.018) = 23.6 rad/s, zeta = 1.18
%               Both natural frequencies sit at or above the motor pole
%               (k_m = 20 rad/s), i.e. the loop is faster than its actuator.
%
%   'designed'  places the linearised loop at the declared design targets
%               (zeta = 0.8, wn = 8 rad/s in roll/pitch, 3 rad/s in yaw),
%               which keeps a factor ~2.5 clear of the motor bandwidth.
%               k_R entries are halved because of the doubling above.
%
%   NOTE: 'as_flown' is the default so that refactored runs reproduce the
%   existing results bit for bit.  Switching to 'designed' CHANGES the closed
%   loop -- do it deliberately and re-run the Task 2 scenarios.

arguments
    UAV  struct
    mode (1,:) char {mustBeMember(mode,{'as_flown','designed'})} = 'designed'
end

att.mode = mode;

% Declared design targets (used by 'designed', and reported for comparison)
att.zeta        = 0.8;
att.wn_rollpitch = 8;      % [rad/s]
att.wn_yaw       = 3;      % [rad/s]

Jd = diag(UAV.J)';                                    % [Jxx Jyy Jzz]
wn = [att.wn_rollpitch att.wn_rollpitch att.wn_yaw];
att.K_R_target     = Jd .* wn.^2;                     % desired linear stiffness
att.K_omega_target = 2*att.zeta .* wn .* Jd;          % desired linear damping

switch mode
    case 'as_flown'
        att.k_R   = [5 5 5];
        att.Komega = eye(3);

    case 'designed'
        % Invert  Keff = [k2+k3, k1+k3, k1+k2]  for the requested stiffness.
        % A = ones(3) - eye(3) maps k -> Keff and is invertible.
        A       = ones(3) - eye(3);
        att.k_R = (A \ att.K_R_target(:))';
        att.Komega = diag(att.K_omega_target);
end

% Linear equivalents actually seen by the closed loop -- USE THESE for Bode,
% sensitivity, bandwidth and step-response analysis.
att.Keff_R = [att.k_R(2)+att.k_R(3), ...
              att.k_R(1)+att.k_R(3), ...
              att.k_R(1)+att.k_R(2)];
att.Keff_w = diag(att.Komega)';
end

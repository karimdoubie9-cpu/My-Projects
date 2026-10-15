function R = unpack(so)
%UNPACK  Pull the logged signals out of a SimulationOutput into plain arrays.
%
%   R = harness.unpack(so)
%
%   Returns a struct with time in R.t and each logged signal as an N-by-k
%   array (time down the rows), regardless of whether the model logged it as
%   k-by-1-by-N or N-by-k.  Optional signals that the model does not log are
%   simply absent from R, so callers can test with isfield.
%
%   Caveat worth knowing: Final.slx logs some signals under the SAME names as
%   its From-Workspace inputs (traj_ref, wind_vel, wind_force, lambda_sig).
%   They stay inside the SimulationOutput object, so they do not clobber
%   anything -- but do not treat so.wind_vel as "the wind I asked for".

CORE     = {'p','p_d','v','v_d'};
% Everything else the model logs. Absent signals are simply skipped, so the
% same unpack works whether or not a given To-Workspace block is enabled.
%   f_bl, f_ad, f_total  baseline, adaptive and summed force in the SAME run
%                       (N-by-3 arrays, with individual *_t timestamps)
%   f_c    commanded force, inertial      [N]   (outer-loop output)
%   T_c    commanded thrust                     (allocation input)
%   T_out  thrust actually produced             (allocation output, post-rotor)
%   w_a    control wrench applied to the body
%   euler  attitude, ZYX                  [rad]
%   omegab body angular velocity          [rad/s]
%   omega_c commanded body angular velocity     (inner-loop reference)
%   Omega_cmd, Omega_limited, Omega_actual      rotor request, clipped request,
%                                               and motor-integrator state
%   ground_contact                              selected plant at the floor
OPTIONAL = {'e_p','euler','omegab','vb','omega_c','lambda_hat', ...
            'yaw_d','e_hat_pred','theta_a_hat_pred', ...
            'f_c','f_bl','f_ad','f_total', ...
            'T_c','T_out','w_a','pb_active','pb_flight_phase', ...
            'ground_contact','Omega_cmd','Omega_limited','Omega_actual'};

for k = 1:numel(CORE)
    n = CORE{k};
    if ~hasSignal(so, n)
        error('harness:unpack:missingSignal', ...
              ['Signal ''%s'' is not in the simulation output. Check that the ' ...
               'To-Workspace block is enabled and named correctly.'], n);
    end
    R.(n) = grab(so, n);
end
R.t = timeOf(so, 'p');

for k = 1:numel(OPTIONAL)
    n = OPTIONAL{k};
    if hasSignal(so, n)
        R.(n)          = grab(so, n);
        R.([n '_t'])   = timeOf(so, n);
    end
end
end

% ---------------------------------------------------------------------
function tf = hasSignal(so, name)
try
    tf = isprop(so, name) || (isstruct(so) && isfield(so, name));
catch
    tf = false;
end
end

function t = timeOf(so, name)
t = so.(name).Time;
t = t(:);
end

function M = grab(so, name)
ts = so.(name);
D  = squeeze(ts.Data);
n  = numel(ts.Time);
if isvector(D)
    M = D(:);
elseif size(D,1) ~= n && size(D,2) == n
    M = D.';
else
    M = D;
end
end

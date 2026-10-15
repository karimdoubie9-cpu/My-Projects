function ref = pack(t, p, v, a, yaw)
%PACK  Assemble a reference into the From-Workspace struct the model expects.
%
%   ref = traj.pack(t, p, v, a, yaw)
%
%   t    1-by-N or N-by-1 time vector (strictly increasing)
%   p,v,a  3-by-N position / velocity / acceleration
%   yaw  1-by-N heading [rad], expected already unwrapped
%
%   The From-Workspace block in Final.slx reads a 10-channel signal laid out as
%   [p(3) v(3) a(3) yaw(1)], so the struct carries ONLY the fields the block
%   understands.  Use traj.pos/vel/acc/yaw to read the channels back out
%   instead of re-slicing the raw matrix at every call site.

arguments
    t   (:,1) double
    p   (3,:) double
    v   (3,:) double
    a   (3,:) double
    yaw (1,:) double
end

N = numel(t);
if size(p,2) ~= N || size(v,2) ~= N || size(a,2) ~= N || numel(yaw) ~= N
    error('traj:pack:sizeMismatch', ...
          'p, v, a and yaw must all have %d samples to match t.', N);
end
if any(diff(t) <= 0)
    error('traj:pack:nonMonotonicTime', ...
          'Time vector must be strictly increasing (From-Workspace requirement).');
end

ref.time                = t;
ref.signals.values      = [p.', v.', a.', yaw(:)];
ref.signals.dimensions  = 10;
end

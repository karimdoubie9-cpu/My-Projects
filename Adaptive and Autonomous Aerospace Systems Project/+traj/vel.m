function V = vel(ref, k)
%VEL  Velocity channels of a reference struct.
%   V = traj.vel(ref)     -> 3-by-N
%   V = traj.vel(ref, k)  -> 3-by-1 sample k
V = ref.signals.values(:, 4:6).';
if nargin > 1, V = V(:, k); end
end

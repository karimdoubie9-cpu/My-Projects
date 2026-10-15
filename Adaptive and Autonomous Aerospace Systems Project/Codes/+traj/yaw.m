function Y = yaw(ref, k)
%YAW  Heading channel of a reference struct.
%   Y = traj.yaw(ref)     -> 1-by-N
%   Y = traj.yaw(ref, k)  -> scalar sample k
Y = ref.signals.values(:, 10).';
if nargin > 1, Y = Y(k); end
end

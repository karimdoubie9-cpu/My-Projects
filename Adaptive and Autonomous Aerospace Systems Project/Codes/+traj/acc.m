function A = acc(ref, k)
%ACC  Acceleration channels of a reference struct.
%   A = traj.acc(ref)     -> 3-by-N
%   A = traj.acc(ref, k)  -> 3-by-1 sample k
A = ref.signals.values(:, 7:9).';
if nargin > 1, A = A(:, k); end
end

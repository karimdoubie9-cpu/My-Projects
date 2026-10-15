function P = pos(ref, k)
%POS  Position channels of a reference struct.
%   P = traj.pos(ref)     -> 3-by-N
%   P = traj.pos(ref, k)  -> 3-by-1 sample k
P = ref.signals.values(:, 1:3).';
if nargin > 1, P = P(:, k); end
end

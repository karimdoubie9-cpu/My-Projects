function S = crossmat(x)
%CROSSMAT  Skew-symmetric matrix such that crossmat(a)*b == cross(a,b).
%   Accepts x as either a 3x1 or 1x3 vector.

S = [    0  -x(3)   x(2);
      x(3)      0  -x(1);
     -x(2)   x(1)      0];
end

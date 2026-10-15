function [F1, F2, F2_pinv] = allocation(UAV)
%ALLOCATION  Control allocation matrices for the hexacopter.
%
%   [F1,F2,F2_pinv] = uav.allocation(UAV)
%
%   F2 (4 x Nr) maps the individual rotor thrusts T_r to the collective
%   wrench acting on the body:
%
%       [ T ; tau_x ; tau_y ; tau_z ] = F2 * [T_1 ... T_6]'
%
%   row 1 : total thrust                       ->   1
%   row 2 : roll  torque from arm geometry     ->   b*sin(gamma_r)
%   row 3 : pitch torque from arm geometry     -> - b*cos(gamma_r)
%   row 4 : yaw   torque from rotor drag       -> - xi_r*sigma
%
%   F1 (Nr x 4) lifts the 4 controlled wrench components into the 6-vector
%   body wrench used by the plant (the two horizontal force components are
%   not directly actuated by a coplanar rotor set).
%
%   F2_pinv is the Moore-Penrose pseudoinverse used by the allocation block
%   (P2CS1); with Nr = 6 > 4 it selects the minimum-norm thrust set.

g     = UAV.gamma(:)';
b     = UAV.b;
sigma = UAV.sigma;
xi    = UAV.xi(:)';

F2 = [ ones(size(g));
        b*sin(g);
       -b*cos(g);
       -xi*sigma ];

F1 = [zeros(2,4);
      eye(4)];

F2_pinv = pinv(F2);
end

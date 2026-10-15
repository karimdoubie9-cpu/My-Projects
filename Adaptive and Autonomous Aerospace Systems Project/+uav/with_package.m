function [UAV_pkg, package] = with_package(UAV, m_p, opts)
%WITH_PACKAGE  Augment the UAV rigid body with a rigidly attached package.
%
%   [UAV_pkg, package] = uav.with_package(UAV, m_p)
%   [UAV_pkg, package] = uav.with_package(UAV, m_p, r_c=..., L=...)
%
%   The package is modelled as a uniform cube of mass m_p and side L, rigidly
%   attached at r_c (body frame).  The combined body is described about the
%   ORIGINAL body-frame origin, so the static moment S becomes non-zero and
%   the mass matrix M is no longer block diagonal -- this coupling is exactly
%   what the baseline controller does not know about (Task 3.1).
%
%   Parallel-axis contribution of a point-like package about the body origin:
%       J_pkg_o = J_pkg_c + m_p*( (r_c'*r_c)*I - r_c*r_c' )
%
%   Aerodynamic damping is left unchanged: the package is small compared with
%   the rotor disc area, so its contribution to D is negligible.
%
%   See also uav.params, harness.config.

arguments
    UAV   struct
    m_p   (1,1) double {mustBeNonnegative}
    opts.r_c (3,1) double = [0.02; 0.02; -0.15]   % [m] attachment point
    opts.L   (1,1) double = 0.2                    % [m] cube side
end

package.m   = m_p;
package.r_c = opts.r_c;
package.L   = opts.L;
package.J   = (m_p * opts.L^2 / 6) * eye(3);       % uniform cube about its CoM

r_c     = opts.r_c;
m_total = UAV.m + m_p;

% Combined centre of mass, expressed in the body frame (row, as in UAV.r_bg)
r_bg = (UAV.m*UAV.r_bg(:) + m_p*r_c) / m_total;

% Package inertia transported to the body origin
J_pkg_o = package.J + m_p*((r_c'*r_c)*eye(3) - r_c*r_c');

UAV_pkg      = UAV;                 % inherit geometry, rotors, allocation
UAV_pkg.m    = m_total;
UAV_pkg.r_bg = r_bg';
UAV_pkg.J    = UAV.J + J_pkg_o;
UAV_pkg.S    = m_total * uav.crossmat(UAV_pkg.r_bg);
UAV_pkg.M    = [m_total*eye(3),  UAV_pkg.S';
                UAV_pkg.S,       UAV_pkg.J];
UAV_pkg.Minv = inv(UAV_pkg.M);
UAV_pkg.D    = UAV.D;
UAV_pkg.package = package;
end

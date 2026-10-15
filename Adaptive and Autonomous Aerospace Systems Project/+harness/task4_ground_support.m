function [w_ground, supported] = task4_ground_support(p_i, v_i, R, w_free, ground_z)
%TASK4_GROUND_SUPPORT Ideal sticking support on a flat landing surface.
%#codegen
% Used WITH the limited position state and the velocity impact reset added
% by harness.build_single_model. Wrench cancellation alone cannot stop a descent.
% R is body-to-inertial; w_free is the body wrench before ground reaction.
% This ideal support cancels force AND moment while resting on the floor.
% It has no compliance, bounce, finite friction cone, or landing-gear geometry.
% Tolerances (as in Antonio's original ground_support): without them,
% round-off in v_i = R*v_b makes v_i(3) flicker just above zero while the
% vehicle rests, the support drops for an instant, gravity adds a little
% downward velocity, and the support (which cancels force, not velocity)
% then freezes it: the vehicle "sinks" at a constant speed while its height
% is held at the floor.
z_tol = 1e-6;
v_tol = 1e-6;
w_ground = zeros(6,1);
supported = 0;
F_free_i = R*w_free(1:3);
if p_i(3) <= ground_z + z_tol && v_i(3) <= v_tol && F_free_i(3) <= 0
    w_ground = -w_free;
    supported = 1;
end
end

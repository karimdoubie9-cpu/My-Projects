function [pb_active, pb_phase] = task4_flight_gate(p_i, v_i, request, t, cfg, contact)
%TASK4_FLIGHT_GATE Measured-state supervisor for one delivery flight.
%#codegen
% Run this function at a DISCRETE sample time in Simulink, outside PBMRAC.
% cfg = [ground_z; enable_height; enable_vz; confirm_time].
% Phase: 0 waiting for lift-off, 1 airborne, 2 touchdown detected (latched).
% contact is the selected plant's floor-contact signal, not a height proxy.
% A separate continuous contact inhibit also disables the subsystem at impact,
% without waiting up to one supervisor sample for this latch to update.

persistent phase candidate_since
if isempty(phase)
    phase = 0;
    candidate_since = -1;
end

if request <= 0
    phase = 0;
    candidate_since = -1;
else
    height = p_i(3) - cfg(1);
    if phase == 0
        if contact <= 0.5 && height >= cfg(2) && v_i(3) >= cfg(3)
            if candidate_since < 0
                candidate_since = t;
            end
            if t - candidate_since >= cfg(4) - 1e-10
                phase = 1;
            end
        else
            candidate_since = -1;
        end
    elseif phase == 1
        % Stay active throughout descent, including the last millimetre.
        if contact > 0.5
            phase = 2;
        end
    end
end

pb_phase = phase;
pb_active = double(request > 0 && phase == 1);
end

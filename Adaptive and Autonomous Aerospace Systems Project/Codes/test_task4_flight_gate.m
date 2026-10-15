function test_task4_flight_gate()
%TEST_TASK4_FLIGHT_GATE Check lift-off dwell and actual-contact latching.
% Logic test only; it does not compile or simulate the rigid-body plant.
% Reset through request=0, which explicitly resets BOTH persistent states.
% Do not rely on clearing a packaged function or an anonymous call handle.
cfg = [0; 0.15; 0.02; 0.10];
cleanup_guard = onCleanup(@() reset_gate(cfg)); %#ok<NASGU>
fprintf('Testing flight gate: %s\n',which('harness.task4_flight_gate'));

% Repeat in the same MATLAB function instance. Each pass ends airborne,
% deliberately exercising reset from an existing state on the next pass.
for pass = 1:2
    % cfg, pass, label, z, vz, request, t, contact, active, phase
    check(cfg,pass,'explicit reset',          0,     0,   0, 0,    1,0,0);
    check(cfg,pass,'on ground',               0,     0,   1, 0,    1,0,0);
    check(cfg,pass,'short clearance starts',  0.2,   0.1, 1, 0.01, 0,0,0);
    check(cfg,pass,'clearance lost',          0.1,   0.1, 1, 0.06, 0,0,0);
    check(cfg,pass,'new candidate',           0.2,   0.1, 1, 0.10, 0,0,0);
    check(cfg,pass,'confirmation incomplete', 0.2,   0.1, 1, 0.19, 0,0,0);
    check(cfg,pass,'confirmed lift-off',      0.2,   0.1, 1, 0.20, 0,1,1);
    check(cfg,pass,'hover',                  10,     0,   1, 1,    0,1,1);
    check(cfg,pass,'descent',                 0.05, -0.1, 1, 2,    0,1,1);
    check(cfg,pass,'last millimetre airborne',0.001,-0.1, 1, 3,    0,1,1);
    % Impact has already reset velocity: vz need not be negative here.
    check(cfg,pass,'actual contact',          0,     0,   1, 3.01, 1,0,2);
    check(cfg,pass,'bounce remains disabled', 0.2,   0.1, 1, 4,    0,0,2);

    % Reset from the landed latch before starting another time sequence.
    check(cfg,pass,'reset after touchdown',   0,     0,   0, 0,    1,0,0);
    check(cfg,pass,'contact blocks candidate',0.2,   0.1, 1, 0.01, 1,0,0);
    check(cfg,pass,'next run candidate',      0.2,   0.1, 1, 0.10, 0,0,0);
    check(cfg,pass,'next run lift-off',       0.2,   0.1, 1, 0.20, 0,1,1);
    check(cfg,pass,'off request in flight',   1,     0.2, 0, 0.21, 0,0,0);
    check(cfg,pass,'candidate after reset',   0.2,   0.1, 1, 0.30, 0,0,0);
    check(cfg,pass,'airborne before repeat',  0.2,   0.1, 1, 0.40, 0,1,1);
end
fprintf('Task 4 flight-gate logic checks passed (2 consecutive passes).\n');
end

function check(cfg,pass,label,z,vz,request,t,contact,expected_active,expected_phase)
[active,phase] = harness.task4_flight_gate([0;0;z],[0;0;vz],request,t,cfg,contact);
assert(isscalar(active) && isscalar(phase),'task4:gateOutputSize', ...
    'Flight gate must return scalar active and phase outputs.');
assert(active==expected_active && phase==expected_phase,'task4:gateTest', ...
    ['Pass %d, %s at t=%.3f s: expected active=%g, phase=%g; ' ...
     'received active=%g, phase=%g. Gate file: %s'], ...
    pass,label,t,expected_active,expected_phase,active,phase, ...
    which('harness.task4_flight_gate'));
end

function reset_gate(cfg)
% Also runs after an assertion failure, leaving the standalone gate reset.
harness.task4_flight_gate(zeros(3,1),zeros(3,1),0,0,cfg,1);
end

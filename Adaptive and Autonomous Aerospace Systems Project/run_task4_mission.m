%% ========================================================================
%  TASK 4 - DELIVERY MISSION IN AN URBAN ENVIRONMENT
%  ------------------------------------------------------------------------

%     2-D mission  traj.mission    route searched at the 20 m flight level,
%                                  around buildings that are all too tall to
%                                  overfly      -> Figures/Task4_mission
%     3-D mission  traj.mission3d  buildings of varying heights, altitude chosen
%                                  by the planner: it climbs over the low B1
%                                  and passes north of the 40 m B4
%                                               -> Figures/Task4_mission_3d

%
%  ------------------------------------------------------------------------
%  STAGE 1 - PLANNING  (traj.mission, three steps)
%
%  1a. ENVIRONMENT. 
%  1b. ROUTE. plannerRRTStar
%  1c. TRAJECTORY. 
%

%  ------------------------------------------------------------------------
%  STAGE 2 - THE MISSION. Split into takeoff / cruise (2-D) or transit (3-D) /
%  landing
%
%  STAGE 3 - UNFORESEEN EVENTS, both injected during the middle phase, where
%  the vehicle has the least margin and the most speed:
%     - a 1-cosine wind gust at ~35% into it;
%     - a step loss of rotor effectiveness at ~70%, held to the end, so the
%       vehicle must complete the delivery AND land degraded.

%  ========================================================================

clear; close all; clc;

viz.defaults();   % Times New Roman, heavier lines (see +viz/defaults.m)

MISSIONS = {'2d', '3d'};     % fly both; remove one to run only the other

% Rotor start-up per mission (harness.config startup): 'stopped' spins the
% rotors up from zero, 'nominal_hover' starts them at the bare-UAV hover speed.
STARTUP = struct('m2d', 'nominal_hover', 'm3d', 'nominal_hover');


S.model        = 'Final';
S.hold         = 10;       % [s] stationary reference after touchdown
S.touchdown    = 0.05;     % [m] in the hold, ask to sit 5 cm below the pad ...
S.touchdown_t  = 1.0;      % [s] ... over 1 s, so every controller makes contact
S.gate_Ts      = 0.01;     % [s] gate sample time (predictor stays continuous)
S.enable_h     = 0.15;     % [m] measured height needed to enable PBMRAC
S.enable_vz    = 0.02;     % [m/s] measured climb rate needed to enable it
S.confirm      = 0.10;     % [s] both conditions must hold this long
S.hover_h      = 1.0;      % [m] takeoff: rise from the pad to this height ...
S.rise_t       = 2.0;      % [s] ... in this time (quintic, rest to rest) ...
S.hover_t      = 5.0;      % [s] ... and hover there before the mission starts.
                           %     PBMRAC switches on during the rise (gate above) and
                           %     its initial transient settles in the hover; the mission starts
                           %     from the hover point and lands as before.

test_task4_flight_gate();
test_task4_contact();      % small Simulink drop / rest / lift-off / re-landing

results = struct();
for i = 1:numel(MISSIONS)
    S.startup = STARTUP.(['m' MISSIONS{i}]);
    results.(['m' MISSIONS{i}]) = fly_mission(MISSIONS{i}, S);
end

%% ---- Side-by-side summary ----------------------------------------------
fprintf('\n================ Task 4 - both missions ================\n');
fprintf('%-6s %8s %8s %9s | %-8s %8s %8s %8s\n', ...
        'case', 'route', 'flight', 'clear.', 'ctrl', 'RMS', 'peak', 'final');
for i = 1:numel(MISSIONS)
    Rs = results.(['m' MISSIONS{i}]);
    for k = 1:2
        if k == 1
            head = sprintf('%-6s %6.1f m %6.1f s %7.2f m', upper(MISSIONS{i}), ...
                           Rs.M.plan.route_length, Rs.M.T_end, Rs.M.clearance);
        else
            head = repmat(' ', 1, 35);
        end
        fprintf('%s | %-8s %6.3f m %6.3f m %6.3f m\n', head, Rs.tags{k}, ...
                Rs.metrics(k).rms, Rs.metrics(k).peak, Rs.metrics(k).final);
    end
end
fprintf('\nTask 4 complete: figures in Figures/Task4_mission and Figures/Task4_mission_3d.\n');


%% ========================================================================
%  One mission: plan, fly both controllers, plot, label, export
%  ========================================================================
function out = fly_mission(MISSION, S)

tag = upper(MISSION);                    % '2D' or '3D'
fprintf('\n\n############ %s MISSION ############\n', tag);
close all;

%% ---- 1. Plan the mission -------------------------------------------
switch MISSION
    case '2d'
        [ref, M]   = traj.mission(cruise_alt = 20, speed = 3.0, seed = 1, takeoff_alt = S.hover_h);
        mid_phase  = 'Cruise';
        export_dir = 'Task4_mission';
    case '3d'
        [ref, M]   = traj.mission3d(speed = 3.0, seed = 1, max_speed = 4.5, takeoff_alt = S.hover_h);
        mid_phase  = 'Transit';
        export_dir = 'Task4_mission_3d';
end
phase_names = {'Ascent', mid_phase, 'Descent'};

% Planning pipeline, step by step (footprints, RRT*, shortcut, densification,
% min-snap, time allocation, dilation, reference interface).
viz.planning(ref, M);

viz.kinematics(ref, phases = M.phases, name = 'Mission reference');

% Planning-side figures: the map the planner saw, the profiles it produced,
% and the phase structure of the result.
viz.map2d(M, ref = ref, name = 'Occupancy map, 2D top view with inflation');
viz.minsnap(ref, phases = M.phases, limit = M.max_speed, ...
            name = 'Minimum-snap reference profiles');
viz.phase_shaded(ref, M.phases, names = phase_names, ...
                 name = 'Planned trajectory with phase shading');

% Takeoff: rise from the pad to the hover point and hover before the mission.
% All mission times in M move by the same offset.
[ref, M] = prepend_takeoff_hover(ref, M, S.hover_h, S.rise_t, S.hover_t);
fprintf('Takeoff: pad -> %.1f m in %.1f s, hover %.1f s; mission starts at %.1f s\n', ...
        S.hover_h, S.rise_t, S.hover_t, M.t_start);

% Hold the pad position for S.hold seconds after the planned touchdown, so
% the landing and the PBMRAC switch-off are observed on the ground.
M.t_land = ref.time(end);
ref = append_landing_hold(ref, S.hold, 0.01, S.touchdown, S.touchdown_t);

%% ---- 2. Unforeseen events, placed inside the middle phase -----------
%  Placed as FRACTIONS of the middle phase rather than at fixed times, so they
%  stay in the right part of the mission if the map or the speed cap changes
%  the total duration.
mid_span = M.t_descend - M.t_climb;
t_gust   = M.t_climb + 0.35*mid_span;
gust_dur = min(30, 0.5*mid_span);        % a 30 s gust would outlast a short 3-D transit
t_fail   = M.t_climb + 0.70*mid_span;

wind = disturb.gust(ref.time, t0 = t_gust, duration = gust_dur, peak = [-12;0;0]);
lam  = disturb.fault(ref.time, t0 = t_fail, lambda = [0.8;0.9;1;0.8;0.8;1]);
events     = [t_gust, t_gust+gust_dur; t_fail, ref.time(end)];
eventNames = {'Wind gust','Actuator degradation'};

fprintf('\nEvents: gust at %.1f s (%.0f s long), degradation from %.1f s\n', ...
        t_gust, gust_dur, t_fail);

%% ---- 3. Fly it, baseline vs adaptive --------------------------------
% Rotor start-up from S.startup; floor height and flight gate switched on.
P = harness.config(package_mass = 0.4, startup = S.startup);
ground_z         = M.p_takeoff(3);
P.pb_gate_Ts     = S.gate_Ts;
P.pb_gate_cfg    = [ground_z; S.enable_h; S.enable_vz; S.confirm];
P.task4_ground_z = ground_z;
P.pb_flight_gate = 1;
fprintf('Rotor start-up: %s (%.1f rad/s) | PBMRAC enabled after %.2f m and %.2f m/s held for %.2f s\n', ...
        S.startup, P.Omega_prop_0(1), S.enable_h, S.enable_vz, S.confirm);

R_bl = harness.scenario(ref, config=P, model=S.model, package=1, mrac=0, ...
                        wind=wind, lambda=lam, label='Baseline');
R_ad = harness.scenario(ref, config=P, model=S.model, package=1, mrac=1, robust=1, ...
                        wind=wind, lambda=lam, label='PBMRAC');

runs = {R_bl, R_ad};
tags = {'Baseline','PBMRAC'};

% Lift-off / touchdown / PBMRAC switching, landing contact, rotor states.
harness.task4_flight_report(runs, tags, ref, M, P.pb_gate_cfg);
harness.task4_motor_report(runs, tags, P.UAV);

%% ---- 4. Figures -----------------------------------------------------
lims = [0 M.map_size(1); 0 M.map_size(2); 0 M.map_size(3)];

viz.path3d(runs, tags, buildings = M.buildings, limits = lims, ...
           view = [40 25; 90 12], ...
           markers = struct('Takeoff', M.p_takeoff, 'Landing', M.p_land, ...
                            'RouteStart', M.p_climb, 'RouteEnd', M.p_descend), ...
           name = 'Delivery mission, 3D flight path');

% Position AND velocity per axis
viz.tracking(runs, tags, events = events, show = 'both', ...
             name = 'Mission tracking');

% Speed profile: magnitude against the 4-5 m/s limit, plus the vertical
% velocity, which shows the phase structure.
viz.speed(runs, tags, limit = M.max_speed, phases = M.phases, ...
          events = events, eventNames = eventNames, ...
          name = 'Mission speed profile');

viz.tracking(runs, tags, overlay = wind, overlayLabel = 'V_{gust}  [m/s]', ...
             name = 'Mission tracking with gust overlay');

viz.axis_errors(runs, tags, events = events, show3d = false, ...
                name = 'Mission per-axis errors');

viz.errors(runs, tags, events = events, eventNames = eventNames, ...
           name = 'Mission error norms');

% Actuator margin across the mission. The degradation from t_fail onward is
% exactly the case where margin can run out: the controller compensates by
% driving the healthy rotors harder, and this figure says whether it still
% had authority left when it reached the landing pad.
viz.actuators(runs, tags, P.UAV, events = events, ...
              name = 'Mission actuator usage');

viz.attitude(runs, tags, events = events, name = 'Mission attitude response');

% Commanded-force diagnostics from the adaptive run, zoomed on the onset of
% the rotor fault. f_bl is the baseline contribution WITHIN R_ad, not the
% force from the separate R_bl flight. harness.unpack returns N-by-3 arrays,
% so take norms across columns.
force_signals = {'f_bl','f_ad','f_total'};
for k = 1:numel(force_signals)
    assert(isfield(R_ad, force_signals{k}), ...
        'run_task4_mission:missingForceSignal', ...
        'Required logged signal %s is missing from the adaptive run.', force_signals{k});
end
norm_fc_bl = vecnorm(R_ad.f_bl, 2, 2);
norm_fc_ad = vecnorm(R_ad.f_ad, 2, 2);
norm_fc    = vecnorm(R_ad.f_total, 2, 2);
fault_win  = [max(0, t_fail - 5), min(ref.time(end), t_fail + 7)];

figure('Name','Mission force magnitudes PBMRAC','Color','w');
plot(R_ad.f_bl_t, norm_fc_bl, 'LineWidth',1.8);
hold on
plot(R_ad.f_ad_t, norm_fc_ad, 'LineWidth',1.8);
plot(R_ad.f_total_t, norm_fc, 'LineWidth',1.8);
xline(t_fail, '--', 'fault', 'HandleVisibility','off');
grid on
legend('|f_{bl}|','|f_{ad}|','|f_{total}|','Location','best');
xlabel('t  [s]'); ylabel('force magnitude  [N]');
title('Commanded force around the rotor fault (PBMRAC)');
xlim(fault_win);

figure('Name','Mission total force components PBMRAC','Color','w');
plot(R_ad.f_total_t, R_ad.f_total(:,3), 'LineWidth',1.8);
hold on
plot(R_ad.f_total_t, R_ad.f_total(:,1), 'LineWidth',1.8);
plot(R_ad.f_total_t, R_ad.f_total(:,2), 'LineWidth',1.8);
xline(t_fail, '--', 'fault', 'HandleVisibility','off');
grid on
legend('f_z','f_x','f_y','Location','best');
xlabel('t  [s]'); ylabel('force  [N]');
title('Total commanded force components around the rotor fault (PBMRAC)');
xlim(fault_win);

% Top-down flight path with the flown paths overlaid
viz.map2d(M, runs = runs, tags = tags, ref = ref, ...
          name = '2D flight path');

p_ref = traj.pos(ref);
s = viz.style();

% Altitude profile against the phase boundaries. In the 3-D mission this is
% where the planner's climbs over the low buildings show up.
figure('Name','Altitude profile','Color','w','Position',[140 140 820 420]);
hold on; grid on
plot(ref.time, p_ref(3,:), '--k', 'LineWidth',1.8, 'DisplayName','Planned');
for k = 1:2
    plot(runs{k}.t, runs{k}.p(:,3), 'LineWidth',2.1, ...
         'Color', s.colors{k}, 'DisplayName', tags{k});
end
xline(M.phases, '--', lower(phase_names(2:3)), 'HandleVisibility','off');
xlabel('t  [s]'); ylabel('altitude  [m]');
title('Altitude across the mission phases'); legend('Location','best');

%% ---- 5. Numbers -----------------------------------------------------
%  Reported per PHASE as well as overall. The phases are not equivalent:
%  takeoff and landing are near-hover with the vehicle at its heaviest and
%  slowest, the middle phase is where the disturbances hit. A single
%  whole-mission RMS would average those very different regimes.
harness.report(runs, tags, window = [t_gust, t_gust+gust_dur], ...
               title = sprintf('Task 4 - %s delivery mission', tag));

fprintf('\n%s mission - per-phase RMS |e_p| [m]:\n', tag);
fprintf('%-10s %10s %10s %10s %10s\n', 'case', 'takeoff', lower(mid_phase), 'landing', 'end hold');
bounds = [0 M.t_climb; M.t_climb M.t_descend; M.t_descend M.t_land; M.t_land ref.time(end)];
met = struct('rms', {}, 'peak', {}, 'final', {}, 'phase', {});
for k = 1:2
    m = harness.metrics(runs{k});
    v = zeros(1,4);
    for ph = 1:4
        idx = m.t >= bounds(ph,1) & m.t <= bounds(ph,2);
        v(ph) = sqrt(mean(m.nep(idx).^2));
    end
    fprintf('%-10s %10.4f %10.4f %10.4f %10.4f\n', tags{k}, v);
    % 'final' = where it ended up vs the pad (not vs the hold reference,
    % which is S.touchdown below the floor on purpose).
    met(k) = struct('rms', m.rms_p, 'peak', m.peak_p, ...
                    'final', norm(runs{k}.p(end,:) - M.p_land(:).'), 'phase', v);
end

%% ---- Label and save the figures -------------------------------------
%  Each figure gets the mission in its window name (hence in the exported file
%  name, "2D - ..." / "3D - ...") and a tag printed in its top-left corner, so a
%  figure is unambiguous even once copied out of its folder.
label_figures(tag);
viz.export(export_dir);
close all;

out = struct('mission', MISSION, 'ref', ref, 'M', M, 'R_bl', R_bl, 'R_ad', R_ad, ...
             'events', events, 'eventNames', {eventNames}, 'tags', {tags}, ...
             'metrics', met);
end

function [ref, M] = prepend_takeoff_hover(ref, M, h, t_rise, t_hover)
%PREPEND_TAKEOFF_HOVER  Pad -> hover point (quintic, rest to rest), hover, then the mission.
% The mission reference must start at rest at the hover point (traj.mission takeoff_alt).
dt = 0.01;
p = traj.pos(ref); v = traj.vel(ref); a = traj.acc(ref); yaw = traj.yaw(ref); t = ref.time(:).';
assert(abs(p(3,1) - (M.p_takeoff(3))) < 1e-9 && h > 0, 'task4:hover', 'Mission must start at the hover point.');
tr = 0:dt:t_rise; if tr(end) < t_rise, tr(end+1) = t_rise; end
s = tr/t_rise; e3 = [0;0;1];
pr = p(:,1) - h*e3 + h*e3*(10*s.^3 - 15*s.^4 + 6*s.^5);
vr = h*e3*(30*s.^2 - 60*s.^3 + 30*s.^4)/t_rise;
ar = h*e3*(60*s - 180*s.^2 + 120*s.^3)/t_rise^2;
th = dt:dt:t_hover; if isempty(th) || th(end) < t_hover, th(end+1) = t_hover; end
nh = numel(th); t0 = t_rise + t_hover;
ref = traj.pack([tr, t_rise + th, t0 + t(2:end)].', ...
    [pr, repmat(p(:,1),1,nh), p(:,2:end)], [vr, zeros(3,nh), v(:,2:end)], ...
    [ar, zeros(3,nh), a(:,2:end)], [repmat(yaw(1),1,numel(tr)+nh), yaw(2:end)]);
M.p_takeoff = M.p_takeoff - h*e3.';          % the vehicle takes off from the pad
M.t_start   = t0;
M.t_climb   = M.t_climb + t0;
M.t_descend = M.t_descend + t0;
M.phases    = M.phases + t0;
M.T_end     = M.T_end + t0;
end

function ref = append_landing_hold(ref, duration, dt, push, t_push)
%APPEND_LANDING_HOLD  Hold the terminal pose after the rest-to-rest reference.
% During the hold the height reference is lowered smoothly by PUSH metres
% (below the floor) over T_PUSH seconds, i.e. the controller is asked to
% "sit down" on the pad. Without it, a controller that tracks the pad height
% exactly (PBMRAC) hovers millimetres above the floor, never touches it, and
% is never switched off. The floor stops the vehicle; the reference does not
% move it anywhere else.
assert(duration > 0 && dt > 0 && t_push > 0 && t_push <= duration, ...
       'task4:hold', 'Invalid hold duration, time step or push time.');
t   = ref.time(:);
tau = (dt:dt:duration).';
if isempty(tau) || tau(end) < duration - 1e-10
    tau(end+1,1) = duration;
end
p = traj.pos(ref); v = traj.vel(ref); a = traj.acc(ref); yaw = traj.yaw(ref);
nh = numel(tau);
% Quintic rest-to-rest ramp 0 -> -push over t_push, then constant.
s  = min(tau.'/t_push, 1);
dz = -push*(10*s.^3 - 15*s.^4 + 6*s.^5);
vz = -push*(30*s.^2 - 60*s.^3 + 30*s.^4)/t_push .* (s < 1);
az = -push*(60*s - 180*s.^2 + 120*s.^3)/t_push^2 .* (s < 1);
ph = repmat(p(:,end),1,nh); ph(3,:) = ph(3,:) + dz;
vh = zeros(3,nh);           vh(3,:) = vz;
ah = zeros(3,nh);           ah(3,:) = az;
ref = traj.pack([t; t(end)+tau], [p, ph], [v, vh], [a, ah], ...
    [yaw(:).', repmat(yaw(end),1,nh)]);
end

function label_figures(tag)
%LABEL_FIGURES  Prefix every open figure's name with the mission and stamp it.
figs = findall(groot, 'Type', 'figure');
for f = reshape(figs, 1, [])
    if ~startsWith(f.Name, [tag ' - '])
        f.Name = sprintf('%s - %s', tag, f.Name);
    end
    annotation(f, 'textbox', [0.005 0.955 0.2 0.04], ...
               'String', sprintf('%s mission', tag), ...
               'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 11, ...
               'Color', [0.09 0.21 0.36], 'FitBoxToText', 'on', ...
               'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end
end

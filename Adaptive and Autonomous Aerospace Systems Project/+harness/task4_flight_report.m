function task4_flight_report(runs, tags, ref, M, cfg)
%TASK4_FLIGHT_REPORT Ground startup, gate events, and flown mission metrics.
% RMS is time weighted on each run's logged time grid. Building clearance is
% the minimum over logged positions; it is not a continuous collision proof.

fprintf('\n=== TASK 4 FLIGHT AND ACTIVATION CHECK ===\n');
for k = 1:numel(runs)
    R = runs{k};
    [t_core,i_core] = unique(R.t(:),'last');
    p_core = R.p(i_core,:); v_core = R.v(i_core,:);
    assert(isfield(R,'pb_active') && isfield(R,'pb_flight_phase'), ...
        'harness:task4:missingGateLog', ...
        'Missing gate logs: use the prepared model and supplied harness.unpack.');
    on = R.pb_active(:,1) > 0.5;
    assert(~on(1),'harness:task4:gateAtStart','PBMRAC must be disabled at grounded startup.');
    fprintf('\n%s\n',tags{k});
    fprintf('  logged interval: [%.12f %.12f] s; reference end %.12f s\n', ...
        t_core(1),t_core(end),ref.time(end));
    i_on = find(on,1,'first');
    if isempty(i_on)
        fprintf('  PBMRAC activation: none\n');
    else
        t_on = R.pb_active_t(i_on);
        p_on = interp1(t_core,p_core,t_on,'linear');
        v_on = interp1(t_core,v_core,t_on,'linear');
        fprintf('  PBMRAC on at %.4f s: height %.4f m, vz %.4f m/s\n', ...
            t_on,p_on(3)-cfg(1),v_on(3));
        if isfield(R,'e_hat_pred')
            ie = find(R.e_hat_pred_t >= t_on,1,'first');
            if ~isempty(ie) && size(R.e_hat_pred,2)==6
                fprintf('  first predictor sample after enable: |ep| %.3e m, |ev| %.3e m/s\n', ...
                    norm(R.e_hat_pred(ie,1:3)),norm(R.e_hat_pred(ie,4:6)));
            end
        end
    end
    if isfield(R,'ground_contact')
        [ct,ci] = unique(R.ground_contact_t(:),'last');
        contact = double(R.ground_contact(ci,1)>0.5);
        % Ignore initial floor support; report landing-phase actual contact.
        i_touch = find(ct>=M.t_descend & contact>0.5,1,'first');
        if isempty(i_touch)
            fprintf('  Actual landing contact not reached before simulation end.\n');
        else
            t_touch = ct(i_touch);
            fprintf('  actual floor contact during landing at %.6f s\n',t_touch);
            before = find(t_core<t_touch,1,'last');
            if ~isempty(before)
                fprintf('  last pre-contact sample: height %+.6f m, vz %+.6f m/s\n', ...
                    p_core(before,3)-cfg(1),v_core(before,3));
            end
            fprintf('  observation after contact: %.4f s\n',t_core(end)-t_touch);
        end
        contact_at_gate = interp1(ct,contact,R.pb_active_t(:),'previous');
        assert(~any(on & contact_at_gate>0.5),'harness:task4:activeInContact', ...
            'PBMRAC must be disabled while the plant is in floor contact.');
    else
        fprintf('  No actual-contact log: these results predate the v3 contact update.\n');
    end
    i_off = [];
    if ~isempty(i_on)
        i_off = find(~on & (1:numel(on)).'>i_on,1,'first');
    end
    if ~isempty(i_off)
        t_off = R.pb_active_t(i_off);
        p_off = interp1(t_core,p_core,t_off,'linear');
        v_off = interp1(t_core,v_core,t_off,'linear');
        fprintf('  PBMRAC off at %.6f s\n',t_off);
        fprintf('  height at off: %+.6f m; vz at off: %+.6f m/s\n', ...
            p_off(3)-cfg(1),v_off(3));
        fprintf('  observation after PB off: %.4f s\n',max(0,t_core(end)-t_off));
    elseif ~isempty(i_on)
        fprintf('  PBMRAC remains enabled at simulation end.\n');
    end
    [h_min,i_min] = min(p_core(:,3)-cfg(1));
    fprintf('  minimum height above floor: %.6f m at t=%.6f s\n',h_min,t_core(i_min));
    ep = vecnorm(R.p(i_core,:)-R.p_d(i_core,:),2,2);
    [ep_peak,i_peak] = max(ep);
    fprintf('  whole-run peak position error: %.6f m at t=%.6f s\n', ...
        ep_peak,t_core(i_peak));
    flown_clearance = box_clearance(R.p,M.buildings);
    fprintf('  sampled building clearance: %.3f m (requested %.3f m)\n', ...
        flown_clearance,M.uav_radius);
    fprintf('  final position error to landing pad: %.4f m\n',norm(R.p(end,:)-M.p_land));
    fprintf('  final speed: %.4f m/s; final vertical velocity: %+.4f m/s\n', ...
        norm(R.v(end,:)),R.v(end,3));
    if min(R.p(:,3))-cfg(1) < -1e-6
        warning('harness:task4:groundPenetration', ...
            '%s penetrates the floor. Inspect the contact model before interpreting landing results.',tags{k});
    end
    if flown_clearance < M.uav_radius
        warning('harness:task4:flownClearance', ...
            '%s enters the requested building-clearance margin.',tags{k});
    end
end
assert(~any(runs{1}.pb_active(:)>0.5),'harness:task4:baselineGate', ...
    'The baseline run unexpectedly enabled PBMRAC.');
if ~any(runs{2}.pb_active(:)>0.5)
    warning('harness:task4:noEnable', ...
        'The adaptive run never enabled PBMRAC; it is not a valid controller comparison.');
end

bounds = [ref.time(1) M.t_climb; M.t_climb M.t_descend; ...
          M.t_descend M.t_land; M.t_land ref.time(end)];
phase_names = {'takeoff','cruise','landing','end hold'};
fprintf('\nTime-weighted per-phase RMS |e_p| [m]:\n');
fprintf('%-12s %10s %10s %10s %10s\n','case','takeoff','cruise','landing','end hold');
for k = 1:numel(runs)
    R = runs{k};
    [t_core,i_core] = unique(R.t(:),'last');
    e2 = sum((R.p(i_core,:)-R.p_d(i_core,:)).^2,2);
    values = zeros(1,size(bounds,1));
    for j = 1:size(bounds,1)
        values(j) = window_rms(t_core,e2,bounds(j,1),bounds(j,2), ...
            tags{k},phase_names{j});
    end
    fprintf('%-12s %10.4f %10.4f %10.4f %10.4f\n',tags{k},values);
end

figure('Name','Task 4 PBMRAC activation','Color','w');
tiledlayout(3,1);
nexttile; hold on; grid on
for k = 1:numel(runs)
    plot(runs{k}.t,runs{k}.p(:,3)-cfg(1),'LineWidth',1.3,'DisplayName',tags{k});
end
yline(cfg(2),'--','enable height','HandleVisibility','off');
ylabel('height above floor [m]'); legend('Location','best');
nexttile; hold on; grid on
for k = 1:numel(runs)
    stairs(runs{k}.pb_active_t,runs{k}.pb_active(:,1),'LineWidth',1.3,'DisplayName',tags{k});
end
ylim([-0.1 1.1]); yticks([0 1]); ylabel('PBMRAC active'); xlabel('t [s]');
nexttile; hold on; grid on
t_zoom = max(ref.time(1),M.t_land-2);
for k = 1:numel(runs)
    R = runs{k};
    keep = R.t>=t_zoom & R.t<=ref.time(end);
    plot(R.t(keep),R.p(keep,3)-cfg(1),'LineWidth',1.3,'DisplayName',tags{k});
    on = R.pb_active(:,1)>0.5;
    i_on = find(on,1,'first'); i_off = [];
    if ~isempty(i_on), i_off=find(~on & (1:numel(on)).'>i_on,1,'first'); end
    if ~isempty(i_off)
        t_off = R.pb_active_t(i_off);
        if t_off>=t_zoom && t_off<=ref.time(end)
            xline(t_off,':','PB off','HandleVisibility','off');
        end
    end
end
yline(0,'--','floor','HandleVisibility','off');
xlim([t_zoom ref.time(end)]);
ylabel('height above floor [m]'); xlabel('t [s]');
title('Landing and final hold'); legend('Location','best');
end

function value = window_rms(t,e2,requested_lo,requested_hi,tag,phase)
% Use only the measured interval. Do not extrapolate missing flight data.
value = NaN;
if requested_hi<=requested_lo || numel(t)<2
    warning('harness:task4:emptyMetricWindow', ...
        '%s %s: no positive-duration interval is available.',tag,phase);
    return
end
lo = max(requested_lo,t(1));
hi = min(requested_hi,t(end));
if hi<=lo
    warning('harness:task4:uncoveredMetricWindow', ...
        '%s %s: requested [%.12f %.12f] s is not covered by the log.', ...
        tag,phase,requested_lo,requested_hi);
    return
end
time_tol = 100*eps(max([1,abs(requested_lo),abs(requested_hi)]));
if lo-requested_lo>time_tol || requested_hi-hi>time_tol
    warning('harness:task4:partialMetricWindow', ...
        ['%s %s: RMS uses the logged interval [%.12f %.12f] s, ' ...
         'within the requested [%.12f %.12f] s. Check the StopTime precision ' ...
         'and logging coverage. No data are extrapolated.'], ...
        tag,phase,lo,hi,requested_lo,requested_hi);
end
tt = unique([lo; t(t>lo & t<hi); hi]);
yy = interp1(t,e2,tt,'linear');
if any(~isfinite(yy))
    warning('harness:task4:nonfiniteMetric', ...
        '%s %s: logged/interpolated error contains nonfinite values.',tag,phase);
    return
end
value = sqrt(trapz(tt,yy)/(hi-lo));
end

function d = box_clearance(p,buildings)
d = inf;
for j = 1:size(buildings,1)
    B = buildings(j,:);
    dx = max([B(1)-p(:,1),p(:,1)-B(2),zeros(size(p,1),1)],[],2);
    dy = max([B(3)-p(:,2),p(:,2)-B(4),zeros(size(p,1),1)],[],2);
    dz = max([B(5)-p(:,3),p(:,3)-B(6),zeros(size(p,1),1)],[],2);
    d = min(d,min(sqrt(dx.^2+dy.^2+dz.^2)));
end
end

function task4_motor_report(runs,tags,UAV)
%TASK4_MOTOR_REPORT Distinguish requested speed, clipping, and motor motion.
% The plant clips the speed REQUEST before its first-order motor dynamics.
% A stopped motor can be below Omega_min while it spins up; that is not an
% actuator upper-limit violation. Command clipping duration is time weighted.
needed = {'Omega_cmd','Omega_limited','Omega_actual'};
fprintf('\n=== TASK 4 MOTOR STATES AND COMMAND CLIPPING ===\n');
fprintf('%-12s %12s %12s %12s %12s %12s\n', ...
    'case','cmd max','limited max','actual min','actual max','clipped %');
have_logs = true;
for k=1:numel(runs)
    R = runs{k};
    if ~all(isfield(R,needed))
        warning('harness:task4:missingMotorLogs', ...
            '%s lacks motor-state logs; rerun with the v3 Task 4 model.',tags{k});
        have_logs = false;
        continue
    end
    [tc,ic] = unique(R.Omega_cmd_t(:),'last');
    cmd = R.Omega_cmd(ic,:);
    clipped = any(cmd>UAV.Omega_max+1e-7 | cmd<UAV.Omega_min-1e-7,2);
    fraction = 100*trapz(tc,double(clipped))/(tc(end)-tc(1));
    actual = R.Omega_actual;
    fprintf('%-12s %12.2f %12.2f %12.2f %12.2f %12.3f\n', ...
        tags{k},max(cmd(:)),max(R.Omega_limited(:)), ...
        min(actual(:)),max(actual(:)),fraction);
    [~,peak_sample] = max(max(cmd,[],2));
    fprintf('  largest speed request at %.6f s; initial actual speed range [%.6f %.6f] rad/s\n', ...
        tc(peak_sample),min(actual(1,:)),max(actual(1,:)));
    if max(actual(:))>UAV.Omega_max+0.01
        warning('harness:task4:actualMotorLimit', ...
            '%s actual motor speed exceeds the upper bound by more than 0.01 rad/s.',tags{k});
    end
end
fprintf('All speeds are rad/s. Clipped %% = time with at least one request outside the limits.\n');
fprintf('The existing actuator-envelope figure uses commanded speeds; actual states are shown here.\n');
if ~have_logs, return; end
figure('Name','Task 4 actual motor speeds and commands','Color','w');
tiledlayout(numel(runs),1);
for k=1:numel(runs)
    R = runs{k}; nexttile; hold on; grid on
    plot(R.Omega_cmd_t,max(R.Omega_cmd,[],2),'--','DisplayName','largest request');
    plot(R.Omega_limited_t,max(R.Omega_limited,[],2),':', ...
        'LineWidth',1.3,'DisplayName','largest clipped request');
    plot(R.Omega_actual_t,max(R.Omega_actual,[],2),'LineWidth',1.5, ...
        'DisplayName','largest actual speed');
    plot(R.Omega_actual_t,min(R.Omega_actual,[],2),'LineWidth',1.2, ...
        'DisplayName','smallest actual speed');
    yline(UAV.Omega_max,'--','maximum','HandleVisibility','off');
    yline(UAV.Omega_min,':','minimum running request','HandleVisibility','off');
    title(tags{k}); ylabel('rotor speed [rad/s]'); xlabel('t [s]');
    legend('Location','best');
end
end

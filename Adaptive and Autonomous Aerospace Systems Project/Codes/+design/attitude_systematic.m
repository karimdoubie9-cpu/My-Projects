function [att, report] = attitude_systematic(UAV, pos, opts)
%ATTITUDE_SYSTEMATIC Search existing geometric PD with ORIGINAL position LQR.
% Returns [] with explicit infeasibility if no candidate meets ALL requirements.
% Does not silently relax margins, retune position, or change the saved model.
if nargin<3,opts=struct;end
folder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'control_analysis');
oldpath=path;addpath(folder);guard=onCleanup(@()path(oldpath)); %#ok<NASGU>
o=diagnosis_options(opts);
check_attitude_signs(UAV,design.attitude(UAV,'as_flown'),o);
[grid,feasible,report]=sweep_attitude_gains(UAV,pos,o);report.candidates=grid;report.linear_candidates=feasible;
att=[];
if isempty(feasible),return;end
if ~o.run_nonlinear
    report.status='unvalidated_linear_candidates';
    report.linear_candidates=feasible;return
end
for k=1:numel(feasible)
    candidate=feasible{k};
    [test,runs]=attitude_isolated_validation(UAV,{candidate},o);
    report.isolated_attempts{k}=test;
    if ~all(test.stable),continue;end
    cascade=build_position_cascade_model(UAV,pos,candidate,o);
    if ~cascade.stable || any(cascade.phase_margin<o.position_PM_min),continue;end
    [full_test,full_runs]=nonlinear_cascade_validation(UAV,pos,{candidate},o);
    report.cascade_attempts{k}=full_test;
    if ~all(full_test.stable),continue;end
    att=candidate; m=cascade.attitude;
    att.metrics.att_wc=cellfun(@(x)x.wc,m);
    att.metrics.att_bw=cellfun(@(x)x.bandwidth,m);
    att.metrics.att_PM=cellfun(@(x)x.phase_margin,m);
    att.metrics.att_GM=cellfun(@(x)x.gain_margin_db,m);
    att.metrics.pos_wc=cascade.wc;att.metrics.pos_PM=cascade.phase_margin;
    att.metrics.isolated_peak_torque_utilization=max(test.max_torque_utilization);
    att.metrics.isolated_stable=true;
    report.selected=att;report.isolated_runs=runs;report.cascade_summary=full_test;report.cascade_runs=full_runs;
    report.status='nonlinear_validated';
    return
end
report.status='nonlinear_validation_rejected_all_tested_candidates';
end


function add_task4_ground_block(model)
%ADD_TASK4_GROUND_BLOCK Insert the ground-support MATLAB Function in both plants.
% Called once by harness.build_single_model. Adds Antonio's ground-support
% block (from Final_Antonio.slx) to both plants, with the same wiring, so that
% harness.build_single_model can then turn it into the contact model:
%
%   p_i (RB out 1), v_i (RB out 5), R (RB out 2), free wrench -> MATLAB Function
%   MATLAB Function out 1 (ground wrench) + free wrench       -> rigid body in 1
%
% The free wrench is whatever fed the rigid body before (gravity + aero +
% actuators + package), so the support cancels exactly that while resting.
% Branches that already have the block are left untouched.
arguments
    model (1,:) char
end
branches = {[model '/UAV/no package'],[model '/UAV/package']};
for k = 1:2
    branch = branches{k};
    if getSimulinkBlockHandle([branch '/MATLAB Function']) ~= -1, continue; end
    rb = find_system(branch,'SearchDepth',1,'BlockType','SubSystem');
    rb = rb(contains(rb,'Rigid Body ') & contains(rb,'dynamics&kinematics'));
    assert(isscalar(rb),'harness:task4:groundLayout', ...
        'Expected one rigid-body subsystem in %s.',branch);
    rp = get_param(rb{1},'PortHandles');

    % Free wrench currently driving the rigid body.
    line = get_param(rp.Inport(1),'Line');
    free = get_param(line,'SrcPortHandle');
    delete_line(line);

    mf = [branch '/MATLAB Function'];
    add_block('simulink/User-Defined Functions/MATLAB Function',mf, ...
        'Position',[480 560 640 680]);
    chart = find(sfroot,'-isa','Stateflow.EMChart','Path',mf);
    chart.Script = ground_support_source();
    sum_blk = add_block('simulink/Math Operations/Sum',[branch '/Task 4 ground sum'], ...
        'Inputs','++','Position',[700 560 725 600]);
    term = add_block('simulink/Sinks/Terminator',[branch '/Task 4 support flag'], ...
        'Position',[700 640 720 660]);
    mp = get_param(mf,'PortHandles');
    sp = get_param(sum_blk,'PortHandles');
    tp = get_param(term,'PortHandles');

    add_line(branch,rp.Outport(1),mp.Inport(1),'autorouting','on');   % p_i
    add_line(branch,rp.Outport(5),mp.Inport(2),'autorouting','on');   % v_i
    add_line(branch,rp.Outport(2),mp.Inport(3),'autorouting','on');   % R
    add_line(branch,free,mp.Inport(4),'autorouting','on');            % free wrench
    add_line(branch,mp.Outport(1),sp.Inport(1),'autorouting','on');
    add_line(branch,free,sp.Inport(2),'autorouting','on');
    add_line(branch,sp.Outport(1),rp.Inport(1),'autorouting','on');
    add_line(branch,mp.Outport(2),tp.Inport(1),'autorouting','on');
end
end

function s = ground_support_source()
% Antonio's original 4-input law; build_single_model replaces it with
% task4_ground_support (5 inputs, floor height) after checking it is there.
s = sprintf(['function [w_ground, contact] = ground_support(p_i, v_i, R, w_free)\n' ...
    '%% Ideal ground support (replaced by harness.task4_ground_support).\n' ...
    'w_ground = zeros(6,1);\n' ...
    'contact  = 0;\n' ...
    'F_free_i = R * w_free(1:3);\n' ...
    'if (p_i(3) <= 1e-6) && (v_i(3) <= 1e-6) && (F_free_i(3) <= 0)\n' ...
    '    w_ground = -w_free;\n' ...
    '    contact = 1;\n' ...
    'end\n' ...
    'end\n']);
end

function build_single_model(model)
%BUILD_SINGLE_MODEL  Put the Task 4 ground, PBMRAC flight gate and logs into Final.slx.
%
%   harness.build_single_model()          % patches Final.slx in place (backup first)
%   harness.build_single_model('Name')
%
%   One model for every task. The Task 4 additions are switched by parameters
%   (harness.config sets the defaults that leave Tasks 1-3 unchanged):
%
%   task4_ground_z  floor height [m]. -Inf (default): no floor, the added
%                   blocks are inert. Finite: ideal sticking floor -- the
%                   position state is limited at the floor, [v_b;omega_b] is
%                   reset on impact and the support law carries the vehicle at
%                   rest until the rotors lift it.
%   pb_flight_gate  variant control of the 'PBMRAC enable' Variant Subsystem.
%                   0 (default): Enable = MRAC_on, exactly as before; the gate
%                   is removed at compile time, so its discrete rate does not
%                   touch the variable-step solver. 1: measured-lift-off gate
%                   (harness.task4_flight_gate at pb_gate_Ts) AND no contact.
%   pb_gate_Ts, pb_gate_cfg   gate sample time and [ground_z; h; vz; confirm].
%
%   Also moves the f_bl / f_ad / f_total logs to the root (full-run logs even
%   while PBMRAC is disabled) and adds pb_active, pb_flight_phase,
%   ground_contact, Omega_cmd, Omega_limited and Omega_actual.
%
%   Blocks reused from the former Task 4 copy: harness.add_task4_ground_block,
%   harness.add_task4_impact, harness.task4_ground_support,
%   harness.task4_flight_gate. The velocity integrator keeps its
%   [vb_0;omegab_0] initial condition (the impact resets to it; Task 4 starts
%   at rest, so that is zero there).
arguments
    model (1,:) char = 'Final'
end
load_system(model);
assert(getSimulinkBlockHandle([model '/PBMRAC enable']) == -1, ...
    'harness:single:built', '%s already has the single-model blocks.', model);
assert(strcmp(get_param(model,'Dirty'),'off'), 'harness:single:unsaved', ...
    'Save %s first.', model);
file = get_param(model,'FileName');
[folder,name] = fileparts(file);
bdir = fullfile(folder,'backup'); if ~isfolder(bdir), mkdir(bdir); end
backup = fullfile(bdir,[name '_before_single_model.slx']);
[ok,msg] = copyfile(file,backup);
assert(ok,'harness:single:backup','Cannot back up %s: %s',file,msg);
fprintf('Backup: %s\n',backup);

try
    ad = [model '/adaptive controller'];
    ap = get_param(ad,'PortHandles');
    p_src = src(ap.Inport(1)); v_src = src(ap.Inport(2)); u_src = src(ap.Inport(3));
    req = src(ap.Enable);
    rb = get_param(req,'Parent');
    assert(strcmp(get_param(rb,'BlockType'),'Constant') && ...
        strcmp(strtrim(get_param(rb,'Value')),'MRAC_on'), 'harness:single:layout', ...
        'Expected the MRAC_on Constant on the adaptive Enable port.');
    enable_line = get_param(ap.Enable,'Line');
    assert(numel(get_param(enable_line,'DstPortHandle')) == 1, 'harness:single:layout', ...
        'MRAC_on has other branches.');

    %% ---- ground: support law, floor limit, impact reset, contact out ----
    harness.add_task4_ground_block(model);
    branches = {[model '/UAV/no package'],[model '/UAV/package']};
    for k = 1:2
        br = branches{k};
        rbs = find_system(br,'SearchDepth',1,'BlockType','SubSystem');
        rbs = rbs(contains(rbs,'Rigid Body ') & contains(rbs,'dynamics&kinematics'));
        assert(isscalar(rbs),'harness:single:layout','Expected one rigid body in %s.',br);
        body = rbs{1};
        ic = regexprep(get_param([body '/Integrator1'],'InitialCondition'),'\s','');
        assert(strcmp(ic,'[vb_0;omegab_0]'),'harness:single:layout', ...
            'Unexpected velocity IC in %s.',body);
        contact = harness.add_task4_impact(body);
        set_param([body '/Integrator1'],'InitialCondition','[vb_0;omegab_0]');
        outport(body,contact,7,[1180 450 1210 470]);
        bp = get_param(body,'PortHandles');
        outport(br,bp.Outport(7),7,[1050 610 1080 630]);

        chart = emchart([br '/MATLAB Function']);
        chart.Script = fileread(fullfile(fileparts(mfilename('fullpath')),'task4_ground_support.m'));
        sp = get_param([br '/MATLAB Function'],'PortHandles');
        assert(numel(sp.Inport) == 5,'harness:single:layout','Support law floor input missing.');
        fl = add_block('simulink/Sources/Constant',[br '/Task 4 floor'], ...
            'Value','task4_ground_z','Position',[480 680 590 710]);
        fp = get_param(fl,'PortHandles');
        add_line(br,fp.Outport(1),sp.Inport(5),'autorouting','on');
    end
    uav = [model '/UAV'];
    assert(strcmp(strtrim(get_param([uav '/Constant1'],'Value')),'package_exist'), ...
        'harness:single:layout','Expected the package selector UAV/Constant1.');
    sel = add_block('simulink/Signal Routing/Switch',[uav '/Task 4 contact select'], ...
        'Criteria','u2 > Threshold','Threshold','0.5','Position',[720 650 760 710]);
    cp = get_param(sel,'PortHandles');
    bare = get_param(branches{1},'PortHandles'); loaded = get_param(branches{2},'PortHandles');
    pk = get_param([uav '/Constant1'],'PortHandles');
    add_line(uav,loaded.Outport(7),cp.Inport(1),'autorouting','on');
    add_line(uav,pk.Outport(1),cp.Inport(2),'autorouting','on');
    add_line(uav,bare.Outport(7),cp.Inport(3),'autorouting','on');
    up0 = get_param(uav,'PortHandles');
    outport(uav,cp.Outport(1),numel(up0.Outport)+1,[820 665 850 685]);
    up = get_param(uav,'PortHandles');
    contact_src = up.Outport(end);

    %% ---- PBMRAC enable: Variant Subsystem (0 direct, 1 flight gate) -----
    vss = [model '/PBMRAC enable'];
    add_block('simulink/Ports & Subsystems/Variant Subsystem',vss,'Position',[210 700 410 825]);
    delete_block([vss '/Subsystem']); delete_block([vss '/In1']); delete_block([vss '/Out1']);
    ins = {'p_i','v_i','request','contact'}; outs = {'active','phase'};
    for k = 1:4, add_block('built-in/Inport',[vss '/' ins{k}],'Port',num2str(k),'Position',[20 40*k 50 40*k+14]); end
    for k = 1:2, add_block('built-in/Outport',[vss '/' outs{k}],'Port',num2str(k),'Position',[400 40*k 430 40*k+14]); end
    set_param(vss,'VariantActivationTime','update diagram');

    direct = [vss '/Direct'];
    add_block('built-in/SubSystem',direct,'Position',[150 30 250 110]);
    ports(direct,ins,outs);
    for k = [1 2 4]
        t = add_block('built-in/Terminator',[direct '/T' num2str(k)],'Position',[150 40*k 170 40*k+14]);
        add_line(direct,[ins{k} '/1'],[get_param(t,'Name') '/1']);
    end
    add_line(direct,'request/1','active/1');
    add_block('built-in/Constant',[direct '/No phase'],'Value','0','Position',[150 200 200 220]);
    add_line(direct,'No phase/1','phase/1');
    set_param(direct,'VariantControl','pb_flight_gate == 0');

    gsys = [vss '/Flight gate'];
    add_block('built-in/SubSystem',gsys,'Position',[150 150 250 230]);
    ports(gsys,ins,outs);
    gate = [gsys '/PBMRAC flight gate'];
    add_block('simulink/User-Defined Functions/MATLAB Function',gate,'Position',[200 40 400 200]);
    chart = emchart(gate);
    chart.Script = fileread(fullfile(fileparts(mfilename('fullpath')),'task4_flight_gate.m'));
    chart.ChartUpdate = 'DISCRETE';
    chart.SampleTime = 'pb_gate_Ts';
    add_block('simulink/Sources/Clock',[gsys '/PB gate clock'],'Position',[100 250 130 270]);
    add_block('simulink/Sources/Constant',[gsys '/PB gate settings'],'Value','pb_gate_cfg', ...
        'Position',[80 300 180 330]);
    add_block('simulink/Logic and Bit Operations/Logical Operator',[gsys '/PB no contact'], ...
        'Operator','NOT','Position',[300 260 350 290]);
    add_block('simulink/Logic and Bit Operations/Logical Operator',[gsys '/PB contact enable'], ...
        'Operator','AND','Inputs','2','Position',[460 60 510 120]);
    add_line(gsys,'p_i/1','PBMRAC flight gate/1','autorouting','on');
    add_line(gsys,'v_i/1','PBMRAC flight gate/2','autorouting','on');
    add_line(gsys,'request/1','PBMRAC flight gate/3','autorouting','on');
    add_line(gsys,'PB gate clock/1','PBMRAC flight gate/4','autorouting','on');
    add_line(gsys,'PB gate settings/1','PBMRAC flight gate/5','autorouting','on');
    add_line(gsys,'contact/1','PBMRAC flight gate/6','autorouting','on');
    add_line(gsys,'contact/1','PB no contact/1','autorouting','on');
    add_line(gsys,'PBMRAC flight gate/1','PB contact enable/1','autorouting','on');
    add_line(gsys,'PB no contact/1','PB contact enable/2','autorouting','on');
    add_line(gsys,'PB contact enable/1','active/1','autorouting','on');
    add_line(gsys,'PBMRAC flight gate/2','phase/1','autorouting','on');
    set_param(gsys,'VariantControl','pb_flight_gate == 1');

    vp = get_param(vss,'PortHandles');
    delete_line(enable_line);
    add_line(model,p_src,vp.Inport(1),'autorouting','on');
    add_line(model,v_src,vp.Inport(2),'autorouting','on');
    add_line(model,req,vp.Inport(3),'autorouting','on');
    add_line(model,contact_src,vp.Inport(4),'autorouting','on');
    add_line(model,vp.Outport(1),ap.Enable,'autorouting','on');

    %% ---- predictor re-seeded from the measured [p;v] on every enable -----
    % With StatesWhenEnabling = reset the predictor restarts from this IC. For
    % MRAC_on from t = 0 it is the plant's own initial [p;v] (= [p_i_0;v_i_0]);
    % with the flight gate it is the state at lift-off, not the pad.
    pred = [ad '/Integrator'];
    ic = regexprep(get_param(pred,'InitialCondition'),'\s','');
    assert(strcmp(ic,'[p_i_0;v_i_0]') && strcmp(get_param(pred,'InitialConditionSource'),'internal'), ...
        'harness:single:layout','Unexpected predictor initial condition.');
    set_param(pred,'InitialConditionSource','external');
    rs = add_block('simulink/Math Operations/Reshape',[ad '/Predictor IC column'], ...
        'OutputDimensionality','Customize','OutputDimensions','[6 1]','Position',[1370 135 1480 175]);
    rsp = get_param(rs,'PortHandles'); ip = get_param(pred,'PortHandles');
    mx = get_param([ad '/Mux'],'PortHandles');
    add_line(ad,mx.Outport(1),rsp.Inport(1),'autorouting','on');
    add_line(ad,rsp.Outport(1),ip.Inport(2),'autorouting','on');

    %% ---- logs -----------------------------------------------------------
    for nm = {'f_bl','f_ad','f_total'}
        s = find_system(ad,'SearchDepth',1,'BlockType','ToWorkspace','VariableName',nm{1});
        assert(isscalar(s),'harness:single:layout','Expected one %s logger in the adaptive controller.',nm{1});
        delete_block(s{1});
    end
    sp = get_param([model '/Sum'],'PortHandles');
    assert(src(sp.Inport(1)) == ap.Outport(1) && src(sp.Inport(2)) == u_src, ...
        'harness:single:layout','Expected root Sum = adaptive output + baseline force.');
    logger(model,'PB log baseline','f_bl',u_src,[460 840 565 865]);
    logger(model,'PB log adaptive','f_ad',ap.Outport(1),[460 885 565 910]);
    logger(model,'PB log total','f_total',sp.Outport(1),[460 930 565 955]);
    logger(model,'PB log active','pb_active',vp.Outport(1),[460 700 565 725]);
    logger(model,'PB log phase','pb_flight_phase',vp.Outport(2),[460 745 565 770]);
    logger(model,'Task 4 contact log','ground_contact',contact_src,[900 1020 1015 1045]);
    motors = [model '/actuators/actuators'];
    logger(motors,'Task 4 command log','Omega_cmd',outp([motors '/omega_prop_cmd']),[240 325 350 350]);
    logger(motors,'Task 4 clipped command log','Omega_limited',outp([motors '/Saturation1']),[425 325 535 350]);
    logger(motors,'Task 4 motor state log','Omega_actual',outp([motors '/Integrator']),[610 325 720 350]);

    save_system(model);
catch err
    close_system(model,0);
    copyfile(backup,file,'f');
    rethrow(err)
end
fprintf('%s: single model built (ground, PBMRAC enable variant, logs).\n',model);
end

function ports(sys,ins,outs)
for k = 1:numel(ins), add_block('built-in/Inport',[sys '/' ins{k}],'Port',num2str(k),'Position',[20 40*k 50 40*k+14]); end
for k = 1:numel(outs), add_block('built-in/Outport',[sys '/' outs{k}],'Port',num2str(k),'Position',[600 40*k 630 40*k+14]); end
end

function outport(sys,source,n,pos)
o = add_block('built-in/Outport',[sys '/ground_contact'],'Port',num2str(n),'Position',pos);
p = get_param(o,'PortHandles');
add_line(sys,source,p.Inport(1),'autorouting','on');
end

function logger(sys,name,var,source,pos)
b = add_block('simulink/Sinks/To Workspace',[sys '/' name],'VariableName',var, ...
    'SaveFormat','Timeseries','MaxDataPoints','inf','SampleTime','-1','Position',pos);
p = get_param(b,'PortHandles');
add_line(sys,source,p.Inport(1),'autorouting','on');
end

function h = outp(block)
p = get_param(block,'PortHandles'); h = p.Outport(1);
end

function h = src(port)
l = get_param(port,'Line');
assert(isscalar(l) && l ~= -1,'harness:single:layout','Unconnected input port.');
h = get_param(l,'SrcPortHandle');
end

function chart = emchart(path)
chart = find(sfroot,'-isa','Stateflow.EMChart','Path',path);
assert(isscalar(chart),'harness:single:layout','Expected one MATLAB Function at %s.',path);
end

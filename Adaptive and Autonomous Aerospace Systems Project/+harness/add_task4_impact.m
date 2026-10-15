function contact_source = add_task4_impact(system)
%ADD_TASK4_IMPACT Add the ideal floor constraint and inelastic impact reset.
% Shared by the real plant patch and the small Simulink contact smoke test.
% Integrator holds p_i [3x1]; Integrator1 holds [v_b;omega_b] [6x1].
% This Task 4 plant starts at rest. Its velocity reset value is also zero.

% Constrain the POSITION STATE, not a plotted signal downstream of it.
set_param([system '/Integrator'],'LimitOutput','on', ...
    'UpperSaturationLimit','[inf;inf;inf]', ...
    'LowerSaturationLimit','[-inf;-inf;task4_ground_z]','ZeroCross','on');
set_param([system '/Integrator1'],'ExternalReset','rising', ...
    'InitialCondition','[0;0;0;0;0;0]','ZeroCross','on');
pp = get_param([system '/Integrator'],'PortHandles');
vp = get_param([system '/Integrator1'],'PortHandles');
assert(numel(vp.Inport)==2,'harness:task4:contactLayout', ...
    'Expected derivative and external-reset inputs on the velocity integrator.');
height = add_block('simulink/Math Operations/Gain',[system '/Task 4 height'], ...
    'Gain','[0 0 1]','Multiplication','Matrix(K*u)', ...
    'Position',[900 430 990 460]);
floor_block = add_block('simulink/Sources/Constant',[system '/Task 4 floor'], ...
    'Value','task4_ground_z','Position',[900 495 990 525]);
contact = add_block('simulink/Logic and Bit Operations/Relational Operator', ...
    [system '/Task 4 floor contact'],'Operator','<=','ZeroCross','on', ...
    'OutDataTypeStr','boolean','Position',[1060 435 1120 495]);
hp = get_param(height,'PortHandles'); fp = get_param(floor_block,'PortHandles');
cp = get_param(contact,'PortHandles');
add_line(system,pp.Outport(1),hp.Inport(1),'autorouting','on');
add_line(system,hp.Outport(1),cp.Inport(1),'autorouting','on');
add_line(system,fp.Outport(1),cp.Inport(2),'autorouting','on');
% The position integrator is non-feedthrough, breaking the path from
% velocity reset -> p_dot -> p -> reset. No artificial delay is used.
add_line(system,cp.Outport(1),vp.Inport(2),'autorouting','on');
contact_source = cp.Outport(1);
end

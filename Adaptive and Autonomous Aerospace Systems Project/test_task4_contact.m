function test_task4_contact()
%TEST_TASK4_CONTACT Simulate a drop, supported rest, lift-off, and second landing.
% Small unit-mass test of the SAME impact blocks/support law used by the patch.
% This is not a full UAV or controller validation. No SLX file is written.
load_system('simulink');
[~,suffix] = fileparts(tempname);
model = ['t4c_' suffix(end-7:end)];   % short: Windows 260-char build paths
new_system(model);
cleanup_guard = onCleanup(@() close_test_model(model)); %#ok<NASGU>
add_block('simulink/Continuous/Integrator',[model '/Integrator'], ...
    'InitialCondition','p_i_0','Position',[620 120 660 160]);
add_block('simulink/Continuous/Integrator',[model '/Integrator1'], ...
    'InitialCondition','zeros(6,1)','Position',[420 120 460 160]);
add_block('simulink/Math Operations/Gain',[model '/Translation'], ...
    'Gain','[eye(3) zeros(3)]','Multiplication','Matrix(K*u)', ...
    'Position',[510 120 570 160]);
add_block('simulink/Math Operations/Sum',[model '/Wrench'], ...
    'Inputs','++','Position',[330 120 355 165]);
add_block('simulink/Sources/From Workspace',[model '/Free wrench'], ...
    'VariableName','contact_test_wrench','Position',[35 115 145 155]);
add_block('simulink/Math Operations/Reshape',[model '/Wrench column'], ...
    'OutputDimensionality','Customize','OutputDimensions','[6 1]', ...
    'Position',[195 115 255 155]);
add_block('simulink/Sources/Constant',[model '/Rotation'], ...
    'Value','eye(3)','Position',[40 260 110 295]);
add_block('simulink/Sources/Constant',[model '/Floor'], ...
    'Value','task4_ground_z','Position',[40 325 130 360]);
support = [model '/Support'];
add_block('simulink/User-Defined Functions/MATLAB Function',support, ...
    'Position',[205 225 345 350]);
chart = find(sfroot,'-isa','Stateflow.EMChart','Path',support);
chart.Script = fileread(which('harness.task4_ground_support'));
add_line(model,'Free wrench/1','Wrench column/1');
add_line(model,'Wrench column/1','Wrench/1');
add_line(model,'Wrench column/1','Support/4');
add_line(model,'Support/1','Wrench/2');
add_line(model,'Wrench/1','Integrator1/1');
add_line(model,'Integrator1/1','Translation/1');
add_line(model,'Translation/1','Integrator/1');
add_line(model,'Translation/1','Support/2');
add_line(model,'Integrator/1','Support/1');
add_line(model,'Rotation/1','Support/3');
add_line(model,'Floor/1','Support/5');
harness.add_task4_impact(model);
for item={'p_test','nu_test'}
    add_block('simulink/Sinks/To Workspace',[model '/' item{1}], ...
        'VariableName',item{1},'SaveFormat','Timeseries','MaxDataPoints','inf');
end
add_line(model,'Integrator/1','p_test/1');
add_line(model,'Integrator1/1','nu_test/1');

wrench.time = [0;1.99;2;3;3.01;5];
wrench.signals.values = zeros(6,6);
wrench.signals.values(:,1) = 0.2; % horizontal load is supported at rest
wrench.signals.values(:,3) = [-9.81;-9.81;2;2;-9.81;-9.81];
wrench.signals.dimensions = 6;
floor_z = 0.4; % Exercise a nonzero ground plane as well.
in = Simulink.SimulationInput(model);
in = in.setVariable('contact_test_wrench',wrench);
in = in.setVariable('task4_ground_z',floor_z);
in = in.setVariable('p_i_0',[0;0;floor_z+1]);
in = in.setModelParameter('StopTime','5','Solver','ode45', ...
    'MaxStep','0.01','RelTol','1e-7','AbsTol','1e-9','ReturnWorkspaceOutputs','on');
out = sim(in);
t = out.p_test.Time(:);
p = samples(out.p_test); nu = samples(out.nu_test);
assert(numel(t)==size(nu,1),'task4:contactTestGrid','Unexpected contact-test log grids.');
h = p(:,3)-floor_z;
assert(min(h)>=-1e-7,'task4:contactPenetration','Contact test penetrated the floor.');
first = find(h<=1e-7,1,'first');
assert(~isempty(first) && abs(t(first)-sqrt(2/9.81))<0.03, ...
    'task4:contactImpact','Drop did not reach the floor at the expected time.');
rest = t>=0.8 & t<=1.8;
assert(any(rest) && max(abs(h(rest)))<1e-6 && ...
    max(vecnorm(nu(rest,:),2,2))<1e-6, ...
    'task4:contactRest','Position and velocity must both remain at supported rest.');
assert(p(find(rest,1),1)>0.005 && ...
    max(vecnorm(diff(p(rest,1:2),1,1),2,2))<1e-6, ...
    'task4:contactPosition','Impact must preserve horizontal position and stop sliding.');
assert(max(h(t>=2 & t<=3.5))>0.5,'task4:contactRelease', ...
    'Contact must release when the free vertical force becomes upward.');
assert(abs(h(end))<1e-6 && norm(nu(end,:))<1e-6, ...
    'task4:contactRelanding','Second touchdown must again stop at the floor.');
fprintf('Task 4 contact simulation passed: drop, rest, release, second landing.\n');
end

function values = samples(ts)
values = squeeze(ts.Data);
if size(values,1)~=numel(ts.Time), values=values.'; end
end

function close_test_model(model)
if bdIsLoaded(model), close_system(model,0); end
end

function fig = kinematics(ref, opts)
%KINEMATICS  Reference position, velocity and acceleration profiles.
%
%   viz.kinematics(ref)
%   viz.kinematics(ref, phases=[t1 t2], name='Mission reference')
%
%   Also prints the peak speed, which is the quantity the exam constrains
%   (no larger than 4-5 m/s) -- worth having on screen every time a new
%   trajectory is generated rather than discovering a violation in the results.

arguments
    ref          struct
    opts.phases  double = []
    opts.name    (1,:) char = 'Reference kinematics'
end

t = ref.time;
p = traj.pos(ref); v = traj.vel(ref); a = traj.acc(ref);

fig = figure('Name', opts.name, 'Color', 'w', 'Position', [120 120 840 640]);

data = {p, v, a};
labs = {'position  [m]','velocity  [m/s]','acceleration  [m/s^2]'};
legs = {{'x','y','z'}, {'v_x','v_y','v_z'}, {'a_x','a_y','a_z'}};

for i = 1:3
    subplot(3,1,i)
    plot(t, data{i}.', 'LineWidth', 2); grid on
    ylabel(labs{i}); legend(legs{i}, 'Location','best');
    for k = 1:numel(opts.phases)
        xline(opts.phases(k), '--k', 'HandleVisibility','off');
    end
    if i == 1, title(opts.name); end
    if i == 3, xlabel('t  [s]'); end
end

speed = vecnorm(v, 2, 1);
fprintf('%s: duration %.1f s, peak speed %.2f m/s, peak accel %.2f m/s^2\n', ...
        opts.name, t(end), max(speed), max(vecnorm(a,2,1)));
if max(speed) > 5
    warning('viz:kinematics:speedLimit', ...
            'Peak speed %.2f m/s exceeds the 4-5 m/s limit set in the exam text.', ...
            max(speed));
end
end

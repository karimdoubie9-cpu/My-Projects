function fig = path3d(runs, tags, opts)
%PATH3D  3-D reference and flown paths, optionally over a building set.
%
%   viz.path3d({R1,R2}, {'Baseline','PBMRAC'})
%   viz.path3d(runs, tags, buildings=B, markers=struct('Takeoff',p1,'Landing',p2))
%
%   buildings is the same M-by-6 [xmin xmax ymin ymax zmin zmax] matrix used by
%   traj.mission, so the Task 4 figures and the planner share one description of
%   the environment instead of two that can drift apart.

arguments
    runs           cell
    tags           cell
    opts.buildings double = zeros(0,6)
    opts.markers   struct = struct()
    opts.view      (:,2) double = [40 25]   % one row per view; N rows -> N panels
    opts.name      (1,:) char = '3D path'
    opts.limits    double = zeros(0,2)
end

% Several azimuths of the same scene read far better side by side than as
% separate figures: a single view of an urban corridor hides whether the path
% passes in front of or behind a building.
s  = viz.style();
nv = size(opts.view, 1);

fig = figure('Name', opts.name, 'Color', 'w', ...
             'Position', [80 80 min(1400, 780*nv) 660]);

for vv = 1:nv
    if nv > 1, subplot(1, nv, vv); end
    hold on; grid on; view(opts.view(vv,1), opts.view(vv,2));

    for b = 1:size(opts.buildings,1)
        viz.buildings(opts.buildings(b,:));
    end

    plot3(runs{1}.p_d(:,1), runs{1}.p_d(:,2), runs{1}.p_d(:,3), s.ref_style, ...
          'Color', s.ref_color, 'LineWidth', s.lw_ref, 'DisplayName', 'Reference');

    for k = 1:numel(runs)
        plot3(runs{k}.p(:,1), runs{k}.p(:,2), runs{k}.p(:,3), 'LineWidth', s.lw, ...
              'Color', s.colors{1+mod(k-1,numel(s.colors))}, 'DisplayName', tags{k});
    end

    names = fieldnames(opts.markers);
    shapes = {'^','v','s','o'};
    for i = 1:numel(names)
        q = opts.markers.(names{i});
        plot3(q(1), q(2), q(3), shapes{1+mod(i-1,4)}, 'MarkerSize', 11, ...
              'MarkerFaceColor', s.colors{2+mod(i,3)}, 'MarkerEdgeColor','k', ...
              'DisplayName', names{i});
    end

    xlabel('x  [m]'); ylabel('y  [m]'); zlabel('z  [m]');
    axis equal
    if ~isempty(opts.limits)
        xlim(opts.limits(1,:)); ylim(opts.limits(2,:)); zlim(opts.limits(3,:));
    end
    if nv > 1
        title(sprintf('view %d  (az %d, el %d)', vv, opts.view(vv,1), opts.view(vv,2)));
    end
    if vv == 1, legend('Location','best'); end
end
end

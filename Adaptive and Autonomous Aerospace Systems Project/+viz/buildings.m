function buildings(B, color, alpha)
%BUILDINGS  Draw one axis-aligned box from [xmin xmax ymin ymax zmin zmax].
%
%   viz.buildings(B)
%   viz.buildings(B, color, alpha)
%
%   Boxes are excluded from the legend (HandleVisibility off) so a scene with
%   three buildings does not push three meaningless entries into it.

arguments
    B     (1,6) double
    color (1,3) double = [0.55 0.55 0.62]
    alpha (1,1) double = 0.85
end

x = B(1:2); y = B(3:4); z = B(5:6);

v = [x(1) y(1) z(1); x(2) y(1) z(1); x(2) y(2) z(1); x(1) y(2) z(1);
     x(1) y(1) z(2); x(2) y(1) z(2); x(2) y(2) z(2); x(1) y(2) z(2)];
f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

patch('Vertices', v, 'Faces', f, 'FaceColor', color, 'FaceAlpha', alpha, ...
      'EdgeColor', [0.25 0.25 0.25], 'HandleVisibility', 'off');
end

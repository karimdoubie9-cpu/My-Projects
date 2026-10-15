function s = style()
%STYLE  One place to define the look of every figure in the project.
%
%   Keeping colours and line widths here means the slides are visually
%   consistent and a change of palette is a one-line edit rather than a
%   search-and-replace across every plotting script.
%
%   Also applies viz.defaults (Times New Roman, heavier axes), so any viz.*
%   function that asks for the style gets the project look even when the
%   calling script forgot to set it.

viz.defaults();

s.colors = { [0.35 0.35 0.35],    ... % 1 baseline   (grey)
             [0.80 0.00 0.80],    ... % 2 adaptive   (magenta)
             [0.00 0.45 0.74],    ... % 3 third case (blue)
             [0.85 0.33 0.10] };      % 4 fourth     (orange)
s.ref_color   = 'k';
s.ref_style   = '--';
s.lw          = 2.2;
s.lw_ref      = 1.7;
s.axis_labels = {'x','y','z'};
s.event_alpha = 0.12;
end

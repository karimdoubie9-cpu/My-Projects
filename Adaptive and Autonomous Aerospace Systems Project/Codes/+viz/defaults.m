function defaults(opts)
%DEFAULTS  Project-wide plot look: Times New Roman and heavier lines.
%
%   viz.defaults()              % apply (idempotent, cheap to call again)
%   viz.defaults(reset=true)    % back to MATLAB's factory defaults
%
%   Sets GRAPHICS-ROOT (groot) defaults, so every figure created afterwards in
%   the MATLAB session picks them up -- plots drawn directly in the task scripts
%   as well as the viz.* functions. It is called by viz.style (which almost every
%   viz.* function uses) and at the top of each run_task*.m script.
%
%   Two things the defaults cannot reach:
%     * values passed explicitly in a plot call win over any default -- the
%       project's explicit line widths are sized to match (viz.style s.lw);
%     * tiledlayout titles ignore groot fonts; viz.export sets those before
%       saving, and control_Design.m passes the font explicitly.
%
%   The defaults last for the whole MATLAB session, so they also apply to other
%   projects opened afterwards. viz.defaults(reset=true) removes them.

arguments
    opts.reset (1,1) logical = false
    opts.font  (1,:) char    = 'Times New Roman'
end

f = opts.font;
settings = {
    % ---- fonts
    'defaultAxesFontName',           f
    'defaultTextFontName',           f
    'defaultLegendFontName',         f
    'defaultColorbarFontName',       f
    'defaultConstantlineFontName',   f     % xline / yline labels
    'defaultTextarrowshapeFontName', f     % annotation(...)
    'defaultTextboxshapeFontName',   f
    'defaultPolaraxesFontName',      f
    'defaultGeoaxesFontName',        f
    'defaultAxesFontSize',           11    % factory 10
    % ---- thickness
    'defaultAxesLineWidth',          1.0   % factory 0.5: axes box and ticks
    'defaultAxesGridLineWidth',      0.6   % factory 0.5
    'defaultLineLineWidth',          1.5   % factory 0.5: lines with no explicit width
    'defaultStairLineWidth',         1.5
    'defaultErrorbarLineWidth',      1.2
    'defaultConstantlineLineWidth',  1.2   % xline / yline
    'defaultHistogramLineWidth',     0.8   % bar edges
};

for k = 1:size(settings, 1)
    if opts.reset
        set(groot, settings{k,1}, 'remove');
    else
        set(groot, settings{k,1}, settings{k,2});
    end
end
end

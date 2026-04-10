%% plot_orth_gap.m — CQRRT linop vs explicit orthogonality gap
%
% Reads diagnostic_*.csv files produced by CQRRT_diagnostic (both file mode
% and generate mode) and produces a bar chart comparing:
%   CQRRT_expl  — backward-stable TRSM path (orth_Q1)
%   CQRRT_linop — linop default: TRSM-with-identity (orth_Q2)
%
% Usage:
%   plot_orth_gap(data_dir)              — standalone, reads from data_dir
%   plot_orth_gap(data_dir, parent_fig)  — embed in existing figure/tab
%
% data_dir should contain diagnostic_*.csv files, one per matrix.
% Files with kappa_target in header are treated as synthetic (sorted by kappa);
% files without are treated as real matrices (sorted last, by filename).

function plot_orth_gap(data_dir, parent_fig)

if nargin < 1 || isempty(data_dir)
    data_dir = fullfile(fileparts(mfilename('fullpath')), ...
                        'benchmark-output', 'orth-gap');
end
if nargin < 2, parent_fig = []; end

% Colorblind-friendly palette (Wong 2011)
color_expl  = [  0 114 178] / 255;   % blue  — [1] expl_trsm
color_linop = [213  94   0] / 255;   % vermilion — [2] expl_inv_trsm (linop default)
default_fontsize = 13;

%% -----------------------------------------------------------------------
%  Collect all diagnostic_*.csv files
%  -----------------------------------------------------------------------
files = dir(fullfile(data_dir, 'diagnostic_*.csv'));
if isempty(files)
    warning('plot_orth_gap: no diagnostic_*.csv found in %s', data_dir);
    return;
end

n_files = numel(files);
labels       = cell(n_files, 1);
kappa_vals   = zeros(n_files, 1);   % NaN for real matrices
expl_orth    = zeros(n_files, 1);
linop_orth   = zeros(n_files, 1);

for i = 1:n_files
    fpath = fullfile(data_dir, files(i).name);

    % ---- Parse header comments ----
    fid = fopen(fpath, 'r');
    kappa_target = NaN;
    cond_A = NaN;
    meta_m = 0; meta_n = 0;
    while true
        line = fgetl(fid);
        if ~ischar(line) || isempty(line) || line(1) ~= '#', break; end
        tok = regexp(line, 'kappa_target=([\d.e+\-]+)', 'tokens');
        if ~isempty(tok), kappa_target = str2double(tok{1}{1}); end
        tok = regexp(line, 'cond_A=([\d.e+\-]+)', 'tokens');
        if ~isempty(tok), cond_A = str2double(tok{1}{1}); end
        tok = regexp(line, 'm=(\d+)', 'tokens');
        if ~isempty(tok), meta_m = str2double(tok{1}{1}); end
        tok = regexp(line, 'n=(\d+)', 'tokens');
        if ~isempty(tok), meta_n = str2double(tok{1}{1}); end
    end
    fclose(fid);

    % ---- Read CSV (skip comment lines) ----
    n_comments = count_comment_lines(fpath);
    opts = detectImportOptions(fpath, 'NumHeaderLines', n_comments);
    T_data = readtable(fpath, opts);

    % Use worst (maximum) across runs — conservative, consistent with diagnostic figures
    expl_orth(i)  = max(T_data.orth_Q1);   % path [1] expl_trsm
    linop_orth(i) = max(T_data.orth_Q2);   % path [2] expl_inv_trsm (linop default)

    kappa_vals(i) = kappa_target;

    if ~isnan(kappa_target)
        exp_k = round(log10(kappa_target));
        labels{i} = sprintf('synthetic%s\\kappa = 10^{%d}', char(10), exp_k);
    else
        % Use matrix name from the Matrix: header comment
        fid2 = fopen(fpath, 'r');
        mat_name = '';
        while true
            ln = fgetl(fid2);
            if ~ischar(ln) || isempty(ln) || ln(1) ~= '#', break; end
            tok = regexp(ln, '# Matrix:\s*(.+)', 'tokens');
            if ~isempty(tok)
                mat_name = strtrim(tok{1}{1});
                [~, mat_name] = fileparts(mat_name);   % basename only
                break;
            end
        end
        fclose(fid2);
        if isempty(mat_name)
            mat_name = strrep(files(i).name, '.csv', '');
        end
        mat_label = strrep(mat_name, '_', '\_');
        if ~isnan(cond_A)
            exp_c    = floor(log10(cond_A));
            mantissa = cond_A / 10^exp_c;
            labels{i} = sprintf('%s%s\\kappa \\approx %.2f\\times10^{%d}', ...
                mat_label, char(10), mantissa, exp_c);
        else
            labels{i} = mat_label;
        end
    end
end

%% -----------------------------------------------------------------------
%  Sort: synthetic by ascending kappa, then real matrices alphabetically
%  -----------------------------------------------------------------------
is_synthetic = ~isnan(kappa_vals);
synth_idx = find(is_synthetic);
real_idx  = find(~is_synthetic);

% Sort synthetic by kappa
[~, ks] = sort(kappa_vals(synth_idx));
synth_idx = synth_idx(ks);

% Sort real by label
real_labels = labels(real_idx);
[~, kr] = sort(real_labels);
real_idx = real_idx(kr);

idx = [synth_idx; real_idx];
labels     = labels(idx);
expl_orth  = expl_orth(idx);
linop_orth = linop_orth(idx);
n_mats     = numel(idx);

%% -----------------------------------------------------------------------
%  Plot
%  -----------------------------------------------------------------------
if isempty(parent_fig)
    figure('Position', [100 100 max(700, n_mats*130) 580]);
    parent_fig = gcf;
end
% 'loose' padding leaves space below the axes for multi-line tick labels
tl = tiledlayout(parent_fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
ax = nexttile(tl);
hold(ax, 'on');

x = 1:n_mats;
bar(ax, x - 0.18, expl_orth,  0.32, 'FaceColor', color_expl,  'EdgeColor', 'none', ...
    'DisplayName', 'CQRRT\_expl');
bar(ax, x + 0.18, linop_orth, 0.32, 'FaceColor', color_linop, 'EdgeColor', 'none', ...
    'DisplayName', 'CQRRT\_linop');

eps_val = eps('double');
yline(ax, eps_val, '--', '\epsilon_{mach}', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'right', 'FontSize', default_fontsize - 2, ...
      'HandleVisibility', 'off');

set(ax, 'YScale', 'log', 'YLim', [eps_val/10, 2], ...
    'XTick', x, 'XTickLabel', {}, 'FontSize', default_fontsize);
ax.Clipping = 'off';
% Multi-line tick labels via text(); Clipping='off' lets them render below the axis
yl = get(ax, 'YLim');
for k = 1:n_mats
    text(ax, x(k), yl(1), labels{k}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', default_fontsize - 1, 'Interpreter', 'tex', 'Clipping', 'off');
end

ylabel(ax, '$\|Q^\top Q - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', default_fontsize + 1);
title(ax, 'CQRRT: Explicit vs LinOp Orthogonality Error', ...
      'FontSize', default_fontsize + 1, 'Interpreter', 'none');
legend(ax, 'Location', 'northwest', 'FontSize', default_fontsize - 1);
grid(ax, 'on');

end  % function


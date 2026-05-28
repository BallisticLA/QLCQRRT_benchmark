%% plot_diagnostic_results.m — Plot CQRRT_diagnostic benchmark results
%
% Usage:
%   plot_diagnostic_results(csv_path)
%   plot_diagnostic_results(csv_path, title_suffix)
%   plot_diagnostic_results(csv_path, title_suffix, orth_tab, step_tab)
%
% Arguments:
%   csv_path     — full path to diagnostic_YYYYMMDD_HHMMSS.csv
%   title_suffix — optional string appended to figure titles
%   orth_tab     — optional uitab for the orthogonality + cond figure
%   step_tab     — optional uitab for the step-by-step divergence figure
%
% Produces two figures:
%   Figure 1: 2x1 — (a) orth_error for all 4 paths  (b) cond(A_pre) and cond(R_sk)
%   Figure 2: 1x1 — step-by-step divergence: MR^pre, G, R^chol, R (paths [1] vs [2])

function [fig_orth, fig_step] = plot_diagnostic_results(csv_path, title_suffix, orth_tab, step_tab)

fig_orth = [];
fig_step = [];
if nargin < 2, title_suffix = ''; end
if nargin < 3, orth_tab = []; end
if nargin < 4, step_tab = []; end

%% ------------------------------------------------------------------
%  Load CSV
%  ------------------------------------------------------------------
n_comments = count_comment_lines(csv_path);
opts = detectImportOptions(csv_path, 'NumHeaderLines', n_comments);
T = readtable(csv_path, opts);

% Read header comments for metadata
fid = fopen(csv_path, 'r');
meta_m = 0; meta_n = 0; meta_cond_A = 0;
for k = 1:n_comments
    line = fgetl(fid);
    if contains(line, 'm=')
        tok = regexp(line, 'm=(\d+)', 'tokens'); if ~isempty(tok), meta_m = str2double(tok{1}{1}); end
        tok = regexp(line, 'n=(\d+)', 'tokens'); if ~isempty(tok), meta_n = str2double(tok{1}{1}); end
    end
    if contains(line, 'cond_A=')
        tok = regexp(line, 'cond_A=([\d.e+\-]+)', 'tokens'); if ~isempty(tok), meta_cond_A = str2double(tok{1}{1}); end
    end
end
fclose(fid);

num_runs = height(T);

% Collect per-run values (4 paths per metric)
orth_Q    = [T.orth_Q1, T.orth_Q2, T.orth_Q3, T.orth_Q4];
cond_Apre = [T.cond_Apre1, T.cond_Apre2, T.cond_Apre3, T.cond_Apre4];
cond_Rsk  = T.cond_Rsk;

rd_Msk_12       = T.rd_Msk_12;
rd_Rsk_12       = T.rd_Rsk_12;
rd_Apre_12      = T.rd_Apre_12_step;
rd_G_12         = T.rd_G_12;
rd_Rchol_12     = T.rd_Rchol_12;
rd_Rfinal_12    = T.rd_Rfinal_12;

% Use worst (maximum) across runs — shows worst-case divergence between methods
orth_best    = max(orth_Q,    [], 1);
cond_A_best  = max(cond_Apre, [], 1);
cond_Rsk_best = max(cond_Rsk);

path_labels = {'[1] expl\_trsm', '[2] expl\_inv\_trsm', '[3] expl\_inv\_trtri', '[4] expl\_inv\_geqp3'};
% Wong colorblind-friendly palette
path_colors = [
      0/255  114/255  178/255;   % #0072B2 blue       - [1] expl_trsm   (CQRRT_expl, stable)
    213/255   94/255    0/255;   % #D55E00 vermilion  - [2] expl_inv_trsm (CQRRT_linop, unstable)
    230/255  159/255    0/255;   % #E69F00 orange     - [3] expl_inv_trtri
      0/255  158/255  115/255;   % #009E73 green      - [4] expl_inv_geqp3 (stable)
];

base_title = sprintf('%d x %d', meta_m, meta_n);
if ~isempty(title_suffix), base_title = [base_title, ' - ', title_suffix]; end

%% ------------------------------------------------------------------
%  Figure 1: Orthogonality + cond
%  ------------------------------------------------------------------
if isempty(orth_tab)
    fig_orth = figure('Position', [100 100 1000 650]);
    parent_orth = fig_orth;
else
    parent_orth = orth_tab;
end
tl1 = tiledlayout(parent_orth, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(tl1, ['CQRRT Diagnostic: Orthogonality - ', base_title], ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');

% --- Subplot 1: orth_error ---
ax1 = nexttile(tl1);
hold(ax1, 'on');
for p = 1:4
    bh = bar(ax1, p, orth_best(p));
    bh.FaceColor = path_colors(p, :);
    bh.FaceAlpha = 0.75;
    % individual run dots
    scatter(ax1, p * ones(num_runs, 1), orth_Q(:, p), 30, path_colors(p, :), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
y1_floor = 1e-16; if min(orth_Q(:)) < 1e-16, y1_floor = min(orth_Q(:)) / 10; end
set(ax1, 'YScale', 'log', 'YLim', [y1_floor, 1], 'YTick', [1e-15, 1e-10, 1e-5, 1], ...
    'XTick', 1:4, 'XTickLabel', path_labels, 'TickLabelInterpreter', 'tex', 'FontSize', 10);
ylabel(ax1, '$\|Q^TQ - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 11);
title(ax1, 'Full-pipeline Orthogonality Error', 'FontSize', 12);
grid(ax1, 'on');

% --- Subplot 2: cond(A_pre) across paths, plus cond(R_sk) as reference ---
ax2 = nexttile(tl1);
hold(ax2, 'on');
for p = 1:4
    bh = bar(ax2, p, cond_A_best(p));
    bh.FaceColor = path_colors(p, :);
    bh.FaceAlpha = 0.75;
end
yline(ax2, meta_cond_A, '--k', sprintf('cond(M) \\approx %.2e', meta_cond_A), ...
      'LineWidth', 1.2, 'LabelHorizontalAlignment', 'right', 'FontSize', 9);
set(ax2, 'YScale', 'log', 'XTick', 1:4, 'XTickLabel', path_labels, ...
         'TickLabelInterpreter', 'tex', 'FontSize', 10);
ylabel(ax2, 'cond(MR^{pre})', 'Interpreter', 'tex', 'FontSize', 11);
title(ax2, 'Condition Number of Preconditioned Matrix MR^{pre}', 'FontSize', 12);
grid(ax2, 'on');

if ~isempty(orth_tab)
    orth_tab.Title = ['Diagnostic: Ortho  ', base_title];
end

%% ------------------------------------------------------------------
%  Figure 2: Step-by-step divergence (paths [1] vs [2])
%  ------------------------------------------------------------------
% Pipeline stages (x-axis): M^sk  R^sk  MR^pre  G  R^chol  R
% M^sk and R^sk: different sketch code paths (sketch_general vs SpGEMM) from same seed.

stage_labels = {'M^{sk}', 'R^{sk}', 'MR^{pre}', 'G', 'R^{chol}', 'R'};
stage_colors = {
    [0.65 0.65 0.65],        % M^sk   - gray (sketch code path diff)
    [0.65 0.65 0.65],        % R^sk   - gray (QR of above)
    [213   94    0] / 255,   % MR^pre - #D55E00 vermilion
    [230  159    0] / 255,   % G      - #E69F00 orange
    [  0  158  115] / 255,   % R^chol - #009E73 green
    [  0  114  178] / 255,   % R      - #0072B2 blue
};

% Collect per-run values for each stage
stage_data = [rd_Msk_12, rd_Rsk_12, rd_Apre_12, rd_G_12, rd_Rchol_12, rd_Rfinal_12];
stage_best = max(stage_data, [], 1);

if isempty(step_tab)
    fig_step = figure('Position', [150 150 1000 500]);
    parent_step = fig_step;
else
    parent_step = step_tab;
end
tl2 = tiledlayout(parent_step, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(tl2, ['Step-by-step Divergence: CQRRT\_expl vs CQRRT\_linop - ', base_title], ...
        'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');

ax3 = nexttile(tl2);
hold(ax3, 'on');
for s = 1:6
    bh = bar(ax3, s, max(stage_best(s), eps));
    bh.FaceColor = stage_colors{s};
    bh.FaceAlpha = 0.8;
    % scatter individual runs
    vals = max(stage_data(:, s), eps);
    scatter(ax3, s * ones(num_runs, 1), vals, 30, stage_colors{s}, 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
set(ax3, 'YScale', 'log', 'XTick', 1:6, 'XTickLabel', stage_labels, ...
         'TickLabelInterpreter', 'tex', 'FontSize', 11);
ylabel(ax3, '$\|x_1 - x_2\| / \|x_1\|$', 'Interpreter', 'latex', 'FontSize', 11);
title(ax3, 'CQRRT_expl vs CQRRT_linop: divergence at each pipeline stage', ...
      'FontSize', 12, 'Interpreter', 'none');
yline(ax3, eps, '--k', '\epsilon_{mach}', 'LineWidth', 1, ...
      'LabelHorizontalAlignment', 'right', 'FontSize', 9);
grid(ax3, 'on');

if ~isempty(step_tab)
    step_tab.Title = ['Diagnostic: Step-by-step  ', base_title];
end

fprintf('Diagnostic plots generated for %s\n', csv_path);
fprintf('  Runs: %d, Matrix: %d x %d, cond(M)=%.2e, cond(R_sk)=%.2e\n', ...
        num_runs, meta_m, meta_n, meta_cond_A, cond_Rsk_best);
fprintf('  orth_error (best): [1]=%.2e  [2]=%.2e  [3]=%.2e  [4]=%.2e\n', ...
        orth_best(1), orth_best(2), orth_best(3), orth_best(4));

end


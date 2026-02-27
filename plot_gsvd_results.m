%% plot_gsvd_results.m — Plot GSVD / Generalized LS benchmark results
%
% Usage:
%   plot_gsvd_results(data_dir, results_csv, svals_csv, breakdown_csv)
%   plot_gsvd_results(data_dir, results_csv, svals_csv, breakdown_csv, plot_mode)
%
% Arguments:
%   data_dir      — directory containing the CSV files
%   results_csv   — filename of *_gsvd_results.csv
%   svals_csv     — filename of *_gsvd_svals.csv
%   breakdown_csv — filename of *_gsvd_breakdown.csv
%   plot_mode     — 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%
% Produces two figures:
%   Figure 1 (2x1): Timing bars, orthogonality
%   Figure 2 (1x3): Runtime breakdown stacked bars per algorithm

function plot_gsvd_results(data_dir, results_csv, svals_csv, breakdown_csv, plot_mode)

if nargin < 5
    plot_mode = 'best_speed';
end

% Mode display string for titles
if strcmp(plot_mode, 'best_speed')
    mode_str = 'Best Speed';
elseif strcmp(plot_mode, 'worst_ortho')
    mode_str = 'Worst Ortho';
else
    mode_str = 'Best Ortho';
end

results_path   = fullfile(data_dir, results_csv);
svals_path     = fullfile(data_dir, svals_csv);
breakdown_path = fullfile(data_dir, breakdown_csv);

%% ------------------------------------------------------------------
%  Load results CSV
%  ------------------------------------------------------------------
n_comments = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_comments);
T = readtable(results_path, opts);

algorithms = T.algorithm;
qr_time    = T.qr_time_us;
orth_error = T.orth_error;
max_orth   = T.max_orth_cols;
total_a    = T.total_a_time_us;
total_b    = T.total_b_time_us;
total_c    = T.total_c_time_us;
ls_err     = T.ls_rel_error;
m_val      = T.m(1);
n_val      = T.n(1);

unique_algs = unique(algorithms, 'stable');
n_algs = numel(unique_algs);

%% ------------------------------------------------------------------
%  Select representative run per algorithm based on plot_mode
%  ------------------------------------------------------------------
sel_idx = zeros(n_algs, 1);

for a = 1:n_algs
    mask = strcmp(algorithms, unique_algs{a});
    indices = find(mask);

    switch plot_mode
        case 'best_speed'
            [~, best] = min(total_a(mask));
            sel_idx(a) = indices(best);
        case 'worst_ortho'
            [~, worst] = max(orth_error(mask));
            sel_idx(a) = indices(worst);
        case 'best_ortho'
            [~, best] = min(orth_error(mask));
            sel_idx(a) = indices(best);
        otherwise
            error('Unknown plot_mode: %s', plot_mode);
    end
end

%% ------------------------------------------------------------------
%  Figure 1: Main comparison (2x1 grid)
%  ------------------------------------------------------------------
figure('Position', [100, 100, 800, 600]);

% --- Subplot 1: Timing bar chart ---
subplot(2, 1, 1);

bar_data = zeros(n_algs, 4);
for a = 1:n_algs
    i = sel_idx(a);
    bar_data(a, 1) = qr_time(i) / 1000;    % ms
    bar_data(a, 2) = total_a(i) / 1000;     % ms
    bar_data(a, 3) = total_b(i) / 1000;     % ms
    bar_data(a, 4) = total_c(i) / 1000;     % ms
end

b = bar(bar_data);
set(gca, 'XTickLabel', unique_algs);
ylabel('Time (ms)', 'FontSize', 12);
title(sprintf('Timing (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
legend({'Q-less QR', 'App (a): Gen. LS', 'App (b): Gen. svals', 'App (c): Gen. svecs'}, ...
       'Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% Color scheme
colors = [0.2 0.4 0.8;   % blue - QR
          0.8 0.3 0.2;   % red - App a
          0.3 0.7 0.3;   % green - App b
          0.7 0.4 0.8];  % purple - App c
for k = 1:4
    b(k).FaceColor = colors(k, :);
end

% --- Subplot 2: Orthogonality comparison ---
subplot(2, 1, 2);

orth_vals = zeros(n_algs, 1);
for a = 1:n_algs
    orth_vals(a) = orth_error(sel_idx(a));
end

bh = bar(orth_vals);
bh.FaceColor = [0.2 0.4 0.8];
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', unique_algs);
ylabel('$\|Q^TQ - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Orthogonality (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);

% Add text labels on bars
for a = 1:n_algs
    text(a, orth_vals(a) * 2, sprintf('%.1e', orth_vals(a)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

sgtitle(sprintf('GSVD Benchmark (%d \\times %d)', m_val, n_val), ...
    'FontSize', 15, 'FontWeight', 'bold');

%% ------------------------------------------------------------------
%  Figure 2: Runtime breakdown (1x3 tiled layout)
%  ------------------------------------------------------------------
if ~isempty(breakdown_csv)
    n_comments_bd = count_comment_lines(breakdown_path);
    bd_opts = detectImportOptions(breakdown_path, 'NumHeaderLines', n_comments_bd);
    T_bd = readtable(breakdown_path, bd_opts);

    bd_algs = T_bd.algorithm;
    % Breakdown times are in columns t0..t12
    bd_times = T_bd{:, 5:end};  % skip m, n, run, algorithm

    % First pass: compute max total time across algorithms for uniform y-axis
    max_total_ms = 0;
    for a = 1:n_algs
        alg_name = unique_algs{a};
        mask = strcmp(bd_algs, alg_name);
        bd_indices = find(mask);
        sel_run = T.run(sel_idx(a));
        bd_runs = T_bd.run(bd_indices);
        match = find(bd_runs == sel_run, 1);
        if isempty(match), match = 1; end
        row = bd_indices(match);
        % Last non-zero entry is the total
        if strcmp(alg_name, 'CQRRT_linop')
            total_ms = bd_times(row, 11) / 1000;
        elseif strcmp(alg_name, 'CholQR')
            total_ms = bd_times(row, 6) / 1000;
        elseif strcmp(alg_name, 'sCholQR3')
            total_ms = bd_times(row, 13) / 1000;
        else
            total_ms = 0;
        end
        max_total_ms = max(max_total_ms, total_ms);
    end
    bd_ylim = [0, max_total_ms * 1.08];

    figure('Position', [150, 150, 1600, 450]);
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for a = 1:n_algs
        alg_name = unique_algs{a};
        mask = strcmp(bd_algs, alg_name);
        bd_indices = find(mask);

        % Use the selected run for this algorithm
        sel_run = T.run(sel_idx(a));
        bd_runs = T_bd.run(bd_indices);
        match = find(bd_runs == sel_run, 1);
        if isempty(match)
            match = 1;  % fallback to first
        end
        row = bd_indices(match);
        times_us = bd_times(row, :);

        nexttile
        if strcmp(alg_name, 'CQRRT_linop')
            % alloc, saso, qr, trtri, linop_precond, linop_gram, trmm_gram, potrf, finalize, rest, total
            vals = times_us(1:10) / 1000;  % ms (exclude total)
            labels_bd = {'Alloc', 'SASO', 'QR', 'TRTRI', 'LinOp(Prec)', 'LinOp(Gram)', 'TRMM', 'Chol', 'Finalize', 'Rest'};
            bd_colors = {[0.6 0.6 0.6], '#2CA02C', '#1F77B4', '#E377C2', '#B15928', '#EDB120', '#FF7F0E', '#7E2F8E', '#17BECF', [0.8 0.8 0.8]};
        elseif strcmp(alg_name, 'CholQR')
            % alloc, materialize, gram, potrf, rest, total
            vals = times_us(1:5) / 1000;
            labels_bd = {'Alloc', 'LinOp', 'Gram', 'Chol', 'Rest'};
            bd_colors = {[0.6 0.6 0.6], '#B15928', '#EDB120', '#7E2F8E', [0.8 0.8 0.8]};
        elseif strcmp(alg_name, 'sCholQR3')
            % alloc, materialize, gram1, potrf1, trsm1, syrk2, potrf2, update2, syrk3, potrf3, update3, rest, total
            vals = times_us(1:12) / 1000;
            labels_bd = {'Alloc', 'LinOp', 'Gram_1', 'Chol_1', 'Trsm_1', 'Syrk_2', 'Chol_2', 'Upd_2', 'Syrk_3', 'Chol_3', 'Upd_3', 'Rest'};
            bd_colors = {[0.6 0.6 0.6], '#B15928', '#EDB120', '#CAB2D6', '#17BECF', '#FDBF6F', '#7E2F8E', '#FB9A99', '#FF7F0E', '#984EA3', '#E31A1C', [0.8 0.8 0.8]};
        else
            continue;
        end

        % Filter out zero-time entries for cleaner legend
        nonzero = vals > 0;
        vals_nz = vals(nonzero);
        labels_nz = labels_bd(nonzero);
        colors_nz = bd_colors(nonzero);

        bplot = bar(1, vals_nz, 'stacked');
        for i = 1:length(colors_nz)
            bplot(i).FaceColor = colors_nz{i};
            bplot(i).FaceAlpha = 0.9;
        end
        ylim(bd_ylim);
        if a == 1
            ylabel('Time (ms)', 'FontSize', 10);
        end
        title(alg_name, 'FontSize', 12, 'FontWeight', 'bold');
        lgd = legend(labels_nz{:});
        lgd.FontSize = 7;
        lgd.Location = 'eastoutside';
        set(gca, 'XTick', []);
        grid on;
        set(gca, 'FontSize', 10);
    end

    sgtitle(sprintf('QR Runtime Breakdown (%d \\times %d, %s)', m_val, n_val, mode_str), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

fprintf('Plots generated for %s\n', results_path);
fprintf('  Mode: %s\n', plot_mode);
fprintf('  Algorithms: %s\n', strjoin(unique_algs, ', '));
fprintf('  Matrix size: %d x %d\n', m_val, n_val);

end

%% ------------------------------------------------------------------
%  Helper: count comment lines at start of file
%  ------------------------------------------------------------------
function n = count_comment_lines(filepath)
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end
    n = 0;
    while true
        line = fgetl(fid);
        if line == -1, break; end
        if line(1) == '#'
            n = n + 1;
        else
            break;
        end
    end
    fclose(fid);
end

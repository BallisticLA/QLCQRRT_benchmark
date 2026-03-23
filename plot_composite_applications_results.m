%% plot_composite_applications_results.m — Plot CQRRT_linop_composite_applications benchmark results
%
% Usage:
%   plot_composite_applications_results(data_dir, results_csv, breakdown_csv)
%   plot_composite_applications_results(data_dir, results_csv, breakdown_csv, plot_mode)
%
% Arguments:
%   data_dir      — directory containing the CSV files
%   results_csv   — filename of *_gsvd_results.csv
%   breakdown_csv — filename of *_gsvd_breakdown.csv
%   plot_mode     — 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%
% Produces two figures:
%   Figure 1 (3x1): Timing bars, orthogonality, memory
%   Figure 2 (1x3): Runtime breakdown stacked bars per algorithm

function plot_composite_applications_results(data_dir, results_csv, breakdown_csv, plot_mode)

if nargin < 4
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
peak_rss   = T.peak_rss_kb;
analytical = T.analytical_kb;
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
%  Figure 1: Main comparison (3x1 grid)
%  ------------------------------------------------------------------
figure('Position', [100, 100, 800, 900]);

% --- Subplot 1: Timing bar chart ---
subplot(3, 1, 1);

bar_data = zeros(n_algs, 4);
for a = 1:n_algs
    i = sel_idx(a);
    bar_data(a, 1) = qr_time(i) / 1000;    % ms
    bar_data(a, 2) = total_a(i) / 1000;     % ms
    bar_data(a, 3) = total_b(i) / 1000;     % ms
    bar_data(a, 4) = total_c(i) / 1000;     % ms
end

b = bar(bar_data);
set(gca, 'XTickLabel', unique_algs, 'TickLabelInterpreter', 'none');
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
subplot(3, 1, 2);

orth_vals = zeros(n_algs, 1);
for a = 1:n_algs
    orth_vals(a) = orth_error(sel_idx(a));
end

bh = bar(orth_vals);
bh.FaceColor = [0.2 0.4 0.8];
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', unique_algs, 'TickLabelInterpreter', 'none');
ylabel('$\|Q^TQ - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('Orthogonality (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11);
% Expand y-axis to fit all bars with headroom (log scale)
yl = ylim;
ylim([yl(1) / 10, max(orth_vals) * 10]);

% --- Subplot 3: Memory comparison (peak RSS vs analytical) ---
subplot(3, 1, 3);

mem_data = zeros(n_algs, 2);
for a = 1:n_algs
    i = sel_idx(a);
    mem_data(a, 1) = peak_rss(i) / 1024;      % MB
    mem_data(a, 2) = analytical(i) / 1024;     % MB
end

bm = bar(mem_data);
set(gca, 'XTickLabel', unique_algs, 'TickLabelInterpreter', 'none');
ylabel('Memory (MB)', 'FontSize', 12);
title('Peak Working Memory', 'FontSize', 14, 'FontWeight', 'bold');
legend({'Peak RSS', 'Analytical'}, 'Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% Color scheme for memory bars
bm(1).FaceColor = [0.2 0.6 0.4];   % teal - RSS
bm(2).FaceColor = [0.8 0.5 0.2];   % orange - analytical

% Ensure y-axis has headroom for all bars
ylim([0, max(mem_data(:)) * 1.15]);

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

    % Per-algorithm breakdown definitions using display groups.
    % Each group is {label, [1-based column indices in bd_times], color}.
    % Columns within a group are summed for the stacked bar.
    % total_col gives the 1-based column index of the total time.
    %
    % Unified color palette (same operation type = same color across all
    % subplots for easy visual comparison):
    %   Alloc       = [0.6 0.6 0.6]  (medium gray)
    %   Sketch      = #1F77B4         (blue)
    %   QR+Tri.Inv  = #17BECF         (teal)
    %   Fwd (A*x)   = #2CA02C         (green)
    %   Adj (A^Tx)  = #EDB120         (gold)
    %   Gemm/TRMM   = #D95F02         (dark orange)
    %   Syrk        = #FDBF6F         (light gold)
    %   Chol        = #7E2F8E         (purple)
    %   Update      = #E31A1C         (red)
    %   Q_Mat       = #4DBEEE         (cyan)
    %   Rest        = [0.8 0.8 0.8]   (light gray)
    %
    % CSV column layouts (1-based indices into bd_times):
    %   CQRRT_linop (11): 1=alloc 2=sketch 3=qr 4=tri_inv 5=fwd 6=adj 7=trmm 8=chol 9=finalize 10=rest 11=total
    %   CholQR      (6):  1=alloc 2=fwd 3=adj 4=chol 5=rest 6=total
    %   sCholQR3    (18): 1=alloc 2=fwd1 3=adj1 4=chol1 5=upd1 6=fwd2 7=adj2 8=gemm2 9=chol2 10=upd2 11=fwd3 12=adj3 13=gemm3 14=chol3 15=upd3 16=q_mat 17=rest 18=total
    %   sCholQR3_basic (15): 1=alloc 2=fwd1 3=adj1 4=chol1 5=trsm1 6=fwd_q 7=syrk2 8=chol2 9=upd2 10=syrk3 11=chol3 12=upd3 13=q_mat 14=rest 15=total

    bd_defs = struct( ...
        'CQRRT_linop', struct('total_col', 11, ...
            'groups', {{
                {'Alloc',      [1],    [0.6 0.6 0.6]}
                {'Sketch',     [2],    '#1F77B4'}
                {'QR+Tri.Inv', [3 4],  '#17BECF'}
                {'Fwd',        [5],    '#2CA02C'}
                {'Adj',        [6],    '#EDB120'}
                {'TRMM',       [7],    '#D95F02'}
                {'Chol',       [8],    '#7E2F8E'}
                {'Finalize+Rest', [9 10], '#E31A1C'}
            }}), ...
        'CholQR', struct('total_col', 6, ...
            'groups', {{
                {'Alloc',      [1],    [0.6 0.6 0.6]}
                {'Fwd',        [2],    '#2CA02C'}
                {'Adj',        [3],    '#EDB120'}
                {'Chol+Rest',  [4 5],  '#7E2F8E'}
            }}), ...
        'sCholQR3', struct('total_col', 18, ...
            'groups', {{
                {'Alloc',   [1],          [0.6 0.6 0.6]}
                {'Fwd',     [2 6 11],     '#2CA02C'}
                {'Adj',     [3 7 12],     '#EDB120'}
                {'Gemm',    [8 13],       '#D95F02'}
                {'Chol',    [4 9 14],     '#7E2F8E'}
                {'Update',  [5 10 15],    '#E31A1C'}
                {'Q Mat',  [16],         '#4DBEEE'}
                {'Rest',    [17],         [0.8 0.8 0.8]}
            }}), ...
        'sCholQR3_basic', struct('total_col', 15, ...
            'groups', {{
                {'Alloc',   [1],          [0.6 0.6 0.6]}
                {'Fwd',     [2 6],        '#2CA02C'}
                {'Adj',     [3],          '#EDB120'}
                {'Syrk',    [7 10],       '#FDBF6F'}
                {'Chol',    [4 8 11],     '#7E2F8E'}
                {'Update',  [5 9 12],     '#E31A1C'}
                {'Q Mat',  [13],         '#4DBEEE'}
                {'Rest',    [14],         [0.8 0.8 0.8]}
            }}) ...
    );

    % First pass: compute max total time across algorithms for uniform y-axis
    max_total_ms = 0;
    for a = 1:n_algs
        alg_name = unique_algs{a};
        if ~isfield(bd_defs, alg_name), continue; end
        def = bd_defs.(alg_name);
        mask = strcmp(bd_algs, alg_name);
        bd_indices = find(mask);
        sel_run = T.run(sel_idx(a));
        bd_runs = T_bd.run(bd_indices);
        match = find(bd_runs == sel_run, 1);
        if isempty(match), match = 1; end
        row = bd_indices(match);
        total_ms = bd_times(row, def.total_col) / 1000;
        max_total_ms = max(max_total_ms, total_ms);
    end
    bd_ylim = [0, max_total_ms * 1.08];

    % Count how many algorithms have breakdown definitions
    n_bd_algs = 0;
    for a = 1:n_algs
        if isfield(bd_defs, unique_algs{a}), n_bd_algs = n_bd_algs + 1; end
    end
    figure('Position', [150, 150, 450 * n_bd_algs, 450]);
    tiledlayout(1, n_bd_algs, 'TileSpacing', 'compact', 'Padding', 'compact');

    for a = 1:n_algs
        alg_name = unique_algs{a};
        if ~isfield(bd_defs, alg_name), continue; end
        def = bd_defs.(alg_name);
        groups = def.groups;

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

        % Sum columns per group
        n_groups = numel(groups);
        group_vals   = zeros(1, n_groups);
        group_labels = cell(1, n_groups);
        group_colors = cell(1, n_groups);
        for g = 1:n_groups
            grp = groups{g};
            group_labels{g} = grp{1};
            group_colors{g} = grp{3};
            cols = grp{2};
            group_vals(g) = sum(times_us(cols)) / 1000;  % ms
        end

        % Filter out zero-sum groups for cleaner legend
        nonzero = group_vals > 0;
        vals_nz    = group_vals(nonzero);
        labels_nz  = group_labels(nonzero);
        colors_nz  = group_colors(nonzero);

        nexttile
        % Flip data so first operation (Alloc) stacks on top of bar,
        % visually matching the legend which reads top-to-bottom.
        vals_bar   = fliplr(vals_nz);
        colors_bar = fliplr(colors_nz);
        bplot = bar(1, vals_bar, 'stacked');
        for i = 1:length(colors_bar)
            bplot(i).FaceColor = colors_bar{i};
            bplot(i).FaceAlpha = 0.9;
        end
        ylim(bd_ylim);
        if a == 1
            ylabel('Time (ms)', 'FontSize', 10);
        end
        title(alg_name, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
        % Legend: flip handles so top-of-bar (Alloc) = top-of-legend
        lgd = legend(flip(bplot), labels_nz);
        lgd.FontSize = 7;
        lgd.Location = 'eastoutside';
        lgd.Interpreter = 'none';
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

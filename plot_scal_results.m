function [] = plot_scal_results(data_dir, results_file, breakdown_file, plot_mode)
% PLOT_SCALING_RESULTS Plot CQRRT\_linop vs CholQR vs sCholQR3 scaling study results
%
% Usage:
%   plot_scal_results(data_dir, results_file, breakdown_file)
%   plot_scal_results(data_dir, results_file, breakdown_file, plot_mode)
%
% plot_mode: 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%   'best_speed'  - for each (m, algorithm), select the run with fastest time
%   'worst_ortho' - for each (m, algorithm), select the run with worst orthogonality error
%   'best_ortho'  - for each (m, algorithm), select the run with best orthogonality error
%
% Example:
%   plot_scal_results('C:\Users\mymel\data\', ...
%       '20260126_143052_scaling_results.csv', ...
%       '20260126_143052_scaling_breakdown.csv', ...
%       'worst_ortho')

    if nargin < 4, plot_mode = 'best_speed'; end

    % Mode display string for titles
    if strcmp(plot_mode, 'best_speed')
        mode_str = 'Best Speed';
    elseif strcmp(plot_mode, 'worst_ortho')
        mode_str = 'Worst Ortho';
    else
        mode_str = 'Best Ortho';
    end

    % Load main results
    % 12 comment lines + 1 column header = 13 header lines
    % (includes prepended "# Total benchmark runtime: ..." line)
    data = readmatrix(fullfile(data_dir, results_file), 'NumHeaderLines', 13);

    % Load breakdown data
    % 16 comment lines + 1 column header = 17 header lines
    % (includes prepended "# Total benchmark runtime: ..." line
    %  and 4 lines describing CQRRT\_linop/CholQR/sCholQR3/CQRRT\_expl columns)
    breakdown = readmatrix(fullfile(data_dir, breakdown_file), 'NumHeaderLines', 17);

    %% Select one row per (m, algorithm) based on plot_mode
    % Results columns (30 total):
    % 1:m, 2:n, 3:run, 4:aspect_ratio, 5:cond_num, 6:density
    % CQRRT\_linop (7-10): orth_error, max_orth_cols, is_orth, time_us
    % CholQR (11-14): orth_error, max_orth_cols, is_orth, time_us
    % sCholQR3 (15-18): orth_error, max_orth_cols, is_orth, time_us
    % CQRRT\_expl (19-22): orth_error, max_orth_cols, is_orth, time_us
    % Memory (23-30): cqrrt_peak_rss_kb, cqrrt_analytical_kb,
    %   cholqr_peak_rss_kb, cholqr_analytical_kb,
    %   scholqr3_peak_rss_kb, scholqr3_analytical_kb,
    %   dense_cqrrt_peak_rss_kb, dense_cqrrt_analytical_kb

    algo_time_cols = [10, 14, 18, 22];   % CQRRT, CholQR, sCholQR3, CQRRT_expl
    algo_orth_cols = [7, 11, 15, 19];

    unique_m = unique(data(:, 1));
    num_sizes = length(unique_m);
    num_algos = 4;

    % sel_idx(i, a) = row index into data for size i, algorithm a
    sel_idx = zeros(num_sizes, num_algos);
    for i = 1:num_sizes
        mask = data(:, 1) == unique_m(i);
        rows = find(mask);
        sub = data(mask, :);
        for a = 1:num_algos
            if strcmp(plot_mode, 'best_speed')
                [~, j] = min(sub(:, algo_time_cols(a)));
            elseif strcmp(plot_mode, 'worst_ortho')
                [~, j] = max(sub(:, algo_orth_cols(a)));
            else  % best_ortho
                [~, j] = min(sub(:, algo_orth_cols(a)));
            end
            sel_idx(i, a) = rows(j);
        end
    end

    m = unique_m;
    n = data(sel_idx(:, 1), 2);
    aspect_ratio = data(1, 4);

    % CQRRT
    cqrrt_orth_error = data(sel_idx(:,1), 7);
    cqrrt_max_orth_cols = data(sel_idx(:,1), 8);
    cqrrt_time = data(sel_idx(:,1), 10) / 1e6;  % Convert microseconds to seconds

    % CholQR
    cholqr_orth_error = data(sel_idx(:,2), 11);
    cholqr_max_orth_cols = data(sel_idx(:,2), 12);
    cholqr_time = data(sel_idx(:,2), 14) / 1e6;

    % sCholQR3
    scholqr3_orth_error = data(sel_idx(:,3), 15);
    scholqr3_max_orth_cols = data(sel_idx(:,3), 16);
    scholqr3_time = data(sel_idx(:,3), 18) / 1e6;

    % CQRRT\_expl
    dense_cqrrt_orth_error = data(sel_idx(:,4), 19);
    dense_cqrrt_max_orth_cols = data(sel_idx(:,4), 20);
    dense_cqrrt_time = data(sel_idx(:,4), 22) / 1e6;

    % Memory (KB -> MB). Use first run of largest m (memory is constant across runs).
    last_m_rows = find(data(:, 1) == unique_m(end));
    last = last_m_rows(1);
    mem_peak_rss = [data(last,23), data(last,25), data(last,27), data(last,29)] / 1024;
    mem_analytical = [data(last,24), data(last,26), data(last,28), data(last,30)] / 1024;
    largest_m = unique_m(end);

    %% Figure 1: Main Comparison (3x1 grid)
    figure('Position', [100, 100, 800, 900]);

    % Plot 1: Orthogonality error vs matrix size
    subplot(3, 1, 1);
    loglog(m, cqrrt_orth_error, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
    hold on;
    loglog(m, cholqr_orth_error, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
    loglog(m, scholqr3_orth_error, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
    loglog(m, dense_cqrrt_orth_error, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
    xlabel('Number of Rows (m)', 'FontSize', 12);
    ylabel('$\|Q^TQ - I\|_F / \sqrt{n}$', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Orthogonality Error vs Matrix Size (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);

    % Plot 2: Execution time vs matrix size
    subplot(3, 1, 2);
    loglog(m, cqrrt_time, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
    hold on;
    loglog(m, cholqr_time, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
    loglog(m, scholqr3_time, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
    loglog(m, dense_cqrrt_time, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
    xlabel('Number of Rows (m)', 'FontSize', 12);
    ylabel('Execution Time (seconds)', 'FontSize', 12);
    title(sprintf('Runtime vs Matrix Size (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 11);

    % Plot 3: Peak Working Memory (bar plot at largest m)
    subplot(3, 1, 3);
    mem_data = [mem_peak_rss; mem_analytical]';  % 4x2 matrix
    b = bar(mem_data);
    b(1).FaceColor = [0.2 0.4 0.8];  % Blue for Peak RSS
    b(2).FaceColor = [1.0 0.5 0.0];  % Orange for Analytical
    set(gca, 'XTickLabel', {'CQRRT\_linop', 'CholQR', 'sCholQR3', 'CQRRT\_expl'});
    ylabel('Memory (MB)', 'FontSize', 12);
    title(sprintf('Peak Working Memory (m = %d)', largest_m), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Peak RSS', 'Analytical', 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);

    %% Runtime Breakdown
    % Breakdown columns (52 total):
    % 1:m, 2:n, 3:run
    % CQRRT\_linop (4-14): alloc, saso, qr, trtri, linop_precond, linop_gram, trmm_gram, potrf, finalize, rest, total
    % CholQR (15-20): alloc, materialize, gram, potrf, rest, total
    % sCholQR3 (21-33): alloc, materialize, gram1, potrf1, trsm1, syrk2, potrf2, update2, syrk3, potrf3, update3, rest, total
    % CQRRT\_expl (34-44): materialize, saso, qr, trtri(=0), precond, gram, trmm_gram(=0), potrf, finalize, rest, total
    % Memory (45-52): cqrrt_peak_rss_kb, cqrrt_analytical_kb, cholqr_peak_rss_kb, cholqr_analytical_kb,
    %   scholqr3_peak_rss_kb, scholqr3_analytical_kb, dense_cqrrt_peak_rss_kb, dense_cqrrt_analytical_kb

    % Select per-algorithm breakdown rows (same indices as main results)
    cqrrt_breakdown = breakdown(sel_idx(:, 1), :);
    cholqr_breakdown = breakdown(sel_idx(:, 2), :);
    scholqr3_breakdown = breakdown(sel_idx(:, 3), :);
    dense_breakdown = breakdown(sel_idx(:, 4), :);

    % Create x-tick labels from m values
    num_labels = min(10, num_sizes);
    label_indices = round(linspace(1, num_sizes, num_labels));
    xtick_labels = cell(1, length(label_indices));
    for i = 1:length(label_indices)
        idx = label_indices(i);
        xtick_labels{i} = sprintf('%d', m(idx));
    end

    % Get breakdown data for all algorithms (vectorized)
    [cqrrt_pct, cqrrt_abs, cqrrt_colors, cqrrt_labels, cqrrt_lgd_fs] = get_cqrrt_scaling(cqrrt_breakdown);
    [cholqr_pct, cholqr_abs, cholqr_colors, cholqr_labels, cholqr_lgd_fs] = get_cholqr_scaling(cholqr_breakdown);
    [scholqr3_pct, scholqr3_abs, scholqr3_colors, scholqr3_labels, scholqr3_lgd_fs] = get_scholqr3_scaling(scholqr3_breakdown);
    [dense_pct, dense_abs, dense_colors, dense_labels, dense_lgd_fs] = get_dense_cqrrt_scaling(dense_breakdown);

    %% Figure 2: Runtime Breakdown - Percentage (1x4 horizontal layout)
    figure('Position', [150, 150, 2000, 450]);
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    render_scaling_breakdown(cqrrt_pct, cqrrt_colors, cqrrt_labels, 'CQRRT\_linop', true, 'Runtime %', [0 100], label_indices, xtick_labels, cqrrt_lgd_fs);
    nexttile
    render_scaling_breakdown(cholqr_pct, cholqr_colors, cholqr_labels, 'CholQR', false, '', [0 100], label_indices, xtick_labels, cholqr_lgd_fs);
    nexttile
    render_scaling_breakdown(scholqr3_pct, scholqr3_colors, scholqr3_labels, 'sCholQR3', false, '', [0 100], label_indices, xtick_labels, scholqr3_lgd_fs);
    nexttile
    render_scaling_breakdown(dense_pct, dense_colors, dense_labels, 'CQRRT\_expl', false, '', [0 100], label_indices, xtick_labels, dense_lgd_fs);

    sgtitle(sprintf('Runtime Breakdown %% (%s, Aspect Ratio %.0f:1)', mode_str, aspect_ratio), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% Figure 3: Runtime Breakdown - Absolute Time (1x4 horizontal layout)
    max_abs_time = max([max(sum(cqrrt_abs, 2)), max(sum(cholqr_abs, 2)), ...
                        max(sum(scholqr3_abs, 2)), max(sum(dense_abs, 2))]);
    abs_ylim = [0, max_abs_time * 1.05];

    figure('Position', [200, 200, 2000, 450]);
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    render_scaling_breakdown(cqrrt_abs, cqrrt_colors, cqrrt_labels, 'CQRRT\_linop', true, 'Runtime (seconds)', abs_ylim, label_indices, xtick_labels, cqrrt_lgd_fs);
    nexttile
    render_scaling_breakdown(cholqr_abs, cholqr_colors, cholqr_labels, 'CholQR', false, '', abs_ylim, label_indices, xtick_labels, cholqr_lgd_fs);
    nexttile
    render_scaling_breakdown(scholqr3_abs, scholqr3_colors, scholqr3_labels, 'sCholQR3', false, '', abs_ylim, label_indices, xtick_labels, scholqr3_lgd_fs);
    nexttile
    render_scaling_breakdown(dense_abs, dense_colors, dense_labels, 'CQRRT\_expl', false, '', abs_ylim, label_indices, xtick_labels, dense_lgd_fs);

    sgtitle(sprintf('Runtime Breakdown - Absolute Time (%s, Aspect Ratio %.0f:1)', mode_str, aspect_ratio), ...
        'FontSize', 13, 'FontWeight', 'bold');

end

function [Data_pct, Data_abs, colors, labels, lgd_fontsize] = get_cqrrt_scaling(breakdown)
    % CQRRT\_linop columns 4-14: alloc, saso, qr, trtri, linop_precond, linop_gram, trmm_gram, potrf, finalize, rest, total
    % Stack order: Other(rest+fin+alloc), Chol, TRMM, LinOp(Gram), LinOp(Prec), TRTRI, QR, skop
    raw = [breakdown(:,13)+breakdown(:,12)+breakdown(:,4), breakdown(:,[11,10,9,8,7,6,5])];
    total = breakdown(:, 14);
    Data_pct = 100 * raw ./ total;
    Data_abs = raw / 1e6;
    colors = {[0.6 0.6 0.6], '#7E2F8E', '#FF7F0E', '#EDB120', '#B15928', '#E377C2', '#1F77B4', '#2CA02C'};
    labels = {'Other', 'Chol', 'TRMM', 'LinOp(Gram)', 'LinOp(Prec)', 'TRTRI', 'QR', 'skop'};
    lgd_fontsize = 7;
end

function [Data_pct, Data_abs, colors, labels, lgd_fontsize] = get_cholqr_scaling(breakdown)
    % CholQR columns 15-20: alloc, materialize, gram, potrf, rest, total
    % Stack order: Other(rest+alloc), Chol, Gram, LinOp
    raw = [breakdown(:,19)+breakdown(:,15), breakdown(:,[18, 17, 16])];
    total = breakdown(:, 20);
    Data_pct = 100 * raw ./ total;
    Data_abs = raw / 1e6;
    colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928'};
    labels = {'Other', 'Chol', 'Gram', 'LinOp'};
    lgd_fontsize = 7;
end

function [Data_pct, Data_abs, colors, labels, lgd_fontsize] = get_scholqr3_scaling(breakdown)
    % sCholQR3 columns 21-33: alloc, materialize, gram1, potrf1, trsm1, syrk2, potrf2, update2, syrk3, potrf3, update3, rest, total
    % Stack order: Other(rest+alloc), Upd_3, Chol_3, Syrk_3, Upd_2, Chol_2, Syrk_2, Trsm_1, Chol_1, Gram, LinOp
    raw = [breakdown(:,32)+breakdown(:,21), breakdown(:, [31, 30, 29, 28, 27, 26, 25, 24, 23, 22])];
    total = breakdown(:, 33);
    Data_pct = 100 * raw ./ total;
    Data_abs = raw / 1e6;
    colors = {[0.6 0.6 0.6], '#E31A1C', '#984EA3', '#FF7F0E', '#FB9A99', '#7E2F8E', '#FDBF6F', '#17BECF', '#CAB2D6', '#EDB120', '#B15928'};
    labels = {'Other', 'Upd_3', 'Chol_3', 'Syrk_3', 'Upd_2', 'Chol_2', 'Syrk_2', 'Trsm_1', 'Chol_1', 'Gram', 'LinOp'};
    lgd_fontsize = 6;
end

function [Data_pct, Data_abs, colors, labels, lgd_fontsize] = get_dense_cqrrt_scaling(breakdown)
    % CQRRT\_expl columns 34-44: materialize, saso, qr, trtri(=0), precond, gram, trmm_gram(=0), potrf, finalize, rest, total
    % Stack order: Fin+Other, Chol, Gram, Precond, QR, skop, Materialize
    fin_plus_rest = breakdown(:, 42) + breakdown(:, 43);
    raw = [fin_plus_rest, breakdown(:, [41, 39, 38, 36, 35, 34])];
    total = breakdown(:, 44);
    Data_pct = 100 * raw ./ total;
    Data_abs = raw / 1e6;
    colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928', '#1F77B4', '#2CA02C', '#17BECF'};
    labels = {'Fin+Other', 'Chol', 'Gram', 'Precond', 'QR', 'skop', 'Materialize'};
    lgd_fontsize = 7;
end

function [] = render_scaling_breakdown(Data, colors, labels, title_str, show_ylabel, ylabel_str, ylim_val, label_indices, xtick_labels, lgd_fontsize)
    bplot = bar(Data, 'stacked');
    for i = 1:length(colors)
        bplot(i).FaceColor = colors{i};
        bplot(i).FaceAlpha = 0.9;
    end
    ylim(ylim_val);
    set(gca, 'XTick', label_indices, 'XTickLabel', xtick_labels, 'FontSize', 9);
    xtickangle(45);
    xlabel('m', 'FontSize', 10);
    if show_ylabel
        ylabel(ylabel_str, 'FontSize', 10);
    end
    title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
    lgd = legend(labels{:});
    lgd.FontSize = lgd_fontsize;
    lgd.Location = 'eastoutside';
end

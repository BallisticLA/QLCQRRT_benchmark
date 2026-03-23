function [] = plot_basic_results(data_dir, results_file, breakdown_file, plot_mode)
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

    if nargin < 4, plot_mode = 'best_speed'; end

    % Mode display string for titles
    if strcmp(plot_mode, 'best_speed')
        mode_str = 'Best Speed';
    elseif strcmp(plot_mode, 'worst_ortho')
        mode_str = 'Worst Ortho';
    else
        mode_str = 'Best Ortho';
    end

    %% Load main results using readtable with named columns
    results_path = fullfile(data_dir, results_file);
    n_comments = count_comment_lines(results_path);
    opts = detectImportOptions(results_path, 'NumHeaderLines', n_comments);
    T = readtable(results_path, opts);

    %% Load breakdown data
    breakdown_path = fullfile(data_dir, breakdown_file);
    n_comments_bd = count_comment_lines(breakdown_path);
    bd_opts = detectImportOptions(breakdown_path, 'NumHeaderLines', n_comments_bd);
    T_bd = readtable(breakdown_path, bd_opts);

    %% Select one row per (m, algorithm) based on plot_mode
    % Algorithm columns for time and orthogonality
    algo_names = {'cqrrt', 'cholqr', 'scholqr3', 'dense_cqrrt'};
    algo_time_fields = {'cqrrt_time_us', 'cholqr_time_us', 'scholqr3_time_us', 'dense_cqrrt_time_us'};
    algo_orth_fields = {'cqrrt_orth_error', 'cholqr_orth_error', 'scholqr3_orth_error', 'dense_cqrrt_orth_error'};

    unique_m = unique(T.m);
    num_sizes = length(unique_m);
    num_algos = 4;

    % sel_idx(i, a) = row index into T for size i, algorithm a
    sel_idx = zeros(num_sizes, num_algos);
    for i = 1:num_sizes
        mask = T.m == unique_m(i);
        rows = find(mask);
        for a = 1:num_algos
            time_vals = T.(algo_time_fields{a})(mask);
            orth_vals = T.(algo_orth_fields{a})(mask);
            if strcmp(plot_mode, 'best_speed')
                [~, j] = min(time_vals);
            elseif strcmp(plot_mode, 'worst_ortho')
                [~, j] = max(orth_vals);
            else  % best_ortho
                [~, j] = min(orth_vals);
            end
            sel_idx(i, a) = rows(j);
        end
    end

    m = unique_m;
    n = T.n(sel_idx(:, 1));
    aspect_ratio = T.aspect_ratio(1);

    % Extract selected data using named columns
    cqrrt_orth_error      = T.cqrrt_orth_error(sel_idx(:,1));
    cqrrt_max_orth_cols   = T.cqrrt_max_orth_cols(sel_idx(:,1));
    cqrrt_time            = T.cqrrt_time_us(sel_idx(:,1)) / 1e6;

    cholqr_orth_error     = T.cholqr_orth_error(sel_idx(:,2));
    cholqr_max_orth_cols  = T.cholqr_max_orth_cols(sel_idx(:,2));
    cholqr_time           = T.cholqr_time_us(sel_idx(:,2)) / 1e6;

    scholqr3_orth_error   = T.scholqr3_orth_error(sel_idx(:,3));
    scholqr3_max_orth_cols = T.scholqr3_max_orth_cols(sel_idx(:,3));
    scholqr3_time         = T.scholqr3_time_us(sel_idx(:,3)) / 1e6;

    dense_cqrrt_orth_error     = T.dense_cqrrt_orth_error(sel_idx(:,4));
    dense_cqrrt_max_orth_cols  = T.dense_cqrrt_max_orth_cols(sel_idx(:,4));
    dense_cqrrt_time           = T.dense_cqrrt_time_us(sel_idx(:,4)) / 1e6;

    % Memory (KB -> MB). Use first run of largest m.
    last_m_rows = find(T.m == unique_m(end));
    last = last_m_rows(1);
    mem_peak_rss = [T.cqrrt_peak_rss_kb(last), T.cholqr_peak_rss_kb(last), ...
                    T.scholqr3_peak_rss_kb(last), T.dense_cqrrt_peak_rss_kb(last)] / 1024;
    mem_analytical = [T.cqrrt_analytical_kb(last), T.cholqr_analytical_kb(last), ...
                      T.scholqr3_analytical_kb(last), T.dense_cqrrt_analytical_kb(last)] / 1024;
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
    % Extract breakdown data using named columns
    cqrrt_bd  = get_cqrrt_scaling(T_bd, sel_idx(:,1));
    cholqr_bd = get_cholqr_scaling(T_bd, sel_idx(:,2));
    scholqr3_bd = get_scholqr3_scaling(T_bd, sel_idx(:,3));
    dense_bd  = get_dense_cqrrt_scaling(T_bd, sel_idx(:,4));

    % Create x-tick labels from m values
    num_labels = min(10, num_sizes);
    label_indices = round(linspace(1, num_sizes, num_labels));
    xtick_labels = cell(1, length(label_indices));
    for i = 1:length(label_indices)
        idx = label_indices(i);
        xtick_labels{i} = sprintf('%d', m(idx));
    end

    %% Figure 2: Runtime Breakdown - Percentage (1x4 horizontal layout)
    figure('Position', [150, 150, 2000, 450]);
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    render_scaling_breakdown(cqrrt_bd.pct, cqrrt_bd.colors, cqrrt_bd.labels, 'CQRRT\_linop', true, 'Runtime %', [0 100], label_indices, xtick_labels, cqrrt_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(cholqr_bd.pct, cholqr_bd.colors, cholqr_bd.labels, 'CholQR', false, '', [0 100], label_indices, xtick_labels, cholqr_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(scholqr3_bd.pct, scholqr3_bd.colors, scholqr3_bd.labels, 'sCholQR3', false, '', [0 100], label_indices, xtick_labels, scholqr3_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(dense_bd.pct, dense_bd.colors, dense_bd.labels, 'CQRRT\_expl', false, '', [0 100], label_indices, xtick_labels, dense_bd.lgd_fs);

    sgtitle(sprintf('Runtime Breakdown %% (%s, Aspect Ratio %.0f:1)', mode_str, aspect_ratio), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% Figure 3: Runtime Breakdown - Absolute Time (1x4 horizontal layout)
    max_abs_time = max([max(sum(cqrrt_bd.abs, 2)), max(sum(cholqr_bd.abs, 2)), ...
                        max(sum(scholqr3_bd.abs, 2)), max(sum(dense_bd.abs, 2))]);
    abs_ylim = [0, max_abs_time * 1.05];

    figure('Position', [200, 200, 2000, 450]);
    tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile
    render_scaling_breakdown(cqrrt_bd.abs, cqrrt_bd.colors, cqrrt_bd.labels, 'CQRRT\_linop', true, 'Runtime (seconds)', abs_ylim, label_indices, xtick_labels, cqrrt_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(cholqr_bd.abs, cholqr_bd.colors, cholqr_bd.labels, 'CholQR', false, '', abs_ylim, label_indices, xtick_labels, cholqr_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(scholqr3_bd.abs, scholqr3_bd.colors, scholqr3_bd.labels, 'sCholQR3', false, '', abs_ylim, label_indices, xtick_labels, scholqr3_bd.lgd_fs);
    nexttile
    render_scaling_breakdown(dense_bd.abs, dense_bd.colors, dense_bd.labels, 'CQRRT\_expl', false, '', abs_ylim, label_indices, xtick_labels, dense_bd.lgd_fs);

    sgtitle(sprintf('Runtime Breakdown - Absolute Time (%s, Aspect Ratio %.0f:1)', mode_str, aspect_ratio), ...
        'FontSize', 13, 'FontWeight', 'bold');

end

%% ------------------------------------------------------------------
%  Breakdown extraction functions (using named columns)
%  ------------------------------------------------------------------

function bd = get_cqrrt_scaling(T_bd, sel_rows)
    % Stack order: Other(rest+fin+alloc), Chol, TRMM, LinOp(Gram), LinOp(Prec), TRTRI, QR, skop
    raw = [T_bd.cqrrt_rest(sel_rows) + T_bd.cqrrt_finalize(sel_rows) + T_bd.cqrrt_alloc(sel_rows), ...
           T_bd.cqrrt_potrf(sel_rows), T_bd.cqrrt_trmm_gram(sel_rows), ...
           T_bd.cqrrt_linop_gram(sel_rows), T_bd.cqrrt_linop_precond(sel_rows), ...
           T_bd.cqrrt_trtri(sel_rows), T_bd.cqrrt_qr(sel_rows), T_bd.cqrrt_saso(sel_rows)];
    total = T_bd.cqrrt_total(sel_rows);
    bd.pct = 100 * raw ./ total;
    bd.abs = raw / 1e6;
    bd.colors = {[0.6 0.6 0.6], '#7E2F8E', '#FF7F0E', '#EDB120', '#B15928', '#E377C2', '#1F77B4', '#2CA02C'};
    bd.labels = {'Other', 'Chol', 'TRMM', 'LinOp(Gram)', 'LinOp(Prec)', 'TRTRI', 'QR', 'skop'};
    bd.lgd_fs = 7;
end

function bd = get_cholqr_scaling(T_bd, sel_rows)
    % Stack order: Other(rest+alloc), Chol, Gram, LinOp
    raw = [T_bd.cholqr_rest(sel_rows) + T_bd.cholqr_alloc(sel_rows), ...
           T_bd.cholqr_potrf(sel_rows), T_bd.cholqr_gram(sel_rows), ...
           T_bd.cholqr_materialize(sel_rows)];
    total = T_bd.cholqr_total(sel_rows);
    bd.pct = 100 * raw ./ total;
    bd.abs = raw / 1e6;
    bd.colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928'};
    bd.labels = {'Other', 'Chol', 'Gram', 'LinOp'};
    bd.lgd_fs = 7;
end

function bd = get_scholqr3_scaling(T_bd, sel_rows)
    % Stack order: Other(rest+alloc), Upd_3, Chol_3, Syrk_3, Upd_2, Chol_2, Syrk_2, Trsm_1, Chol_1, Gram, LinOp
    raw = [T_bd.scholqr3_rest(sel_rows) + T_bd.scholqr3_alloc(sel_rows), ...
           T_bd.scholqr3_update3(sel_rows), T_bd.scholqr3_potrf3(sel_rows), ...
           T_bd.scholqr3_syrk3(sel_rows), T_bd.scholqr3_update2(sel_rows), ...
           T_bd.scholqr3_potrf2(sel_rows), T_bd.scholqr3_syrk2(sel_rows), ...
           T_bd.scholqr3_trsm1(sel_rows), T_bd.scholqr3_potrf1(sel_rows), ...
           T_bd.scholqr3_gram1(sel_rows), T_bd.scholqr3_materialize(sel_rows)];
    total = T_bd.scholqr3_total(sel_rows);
    bd.pct = 100 * raw ./ total;
    bd.abs = raw / 1e6;
    bd.colors = {[0.6 0.6 0.6], '#E31A1C', '#984EA3', '#FF7F0E', '#FB9A99', '#7E2F8E', '#FDBF6F', '#17BECF', '#CAB2D6', '#EDB120', '#B15928'};
    bd.labels = {'Other', 'Upd_3', 'Chol_3', 'Syrk_3', 'Upd_2', 'Chol_2', 'Syrk_2', 'Trsm_1', 'Chol_1', 'Gram', 'LinOp'};
    bd.lgd_fs = 6;
end

function bd = get_dense_cqrrt_scaling(T_bd, sel_rows)
    % Stack order: Fin+Other, Chol, Gram, Precond, QR, skop, Materialize
    raw = [T_bd.dense_finalize(sel_rows) + T_bd.dense_rest(sel_rows), ...
           T_bd.dense_potrf(sel_rows), T_bd.dense_gram(sel_rows), ...
           T_bd.dense_precond(sel_rows), T_bd.dense_qr(sel_rows), ...
           T_bd.dense_saso(sel_rows), T_bd.dense_materialize(sel_rows)];
    total = T_bd.dense_total(sel_rows);
    bd.pct = 100 * raw ./ total;
    bd.abs = raw / 1e6;
    bd.colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928', '#1F77B4', '#2CA02C', '#17BECF'};
    bd.labels = {'Fin+Other', 'Chol', 'Gram', 'Precond', 'QR', 'skop', 'Materialize'};
    bd.lgd_fs = 7;
end

%% ------------------------------------------------------------------
%  Rendering helper
%  ------------------------------------------------------------------

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

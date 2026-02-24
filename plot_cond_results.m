function [] = plot_cond_results(data_dir, sparse_file, sparse_breakdown_file, plot_mode, dense_file, dense_breakdown_file)
% PLOT_CONDITIONING_COMPARISON Plot CQRRT\_linop vs CholQR vs sCholQR3 conditioning study results
%
% Usage (sparse only):
%   plot_cond_results(data_dir, sparse_file, sparse_breakdown_file)
%   plot_cond_results(data_dir, sparse_file, sparse_breakdown_file, plot_mode)
%
% Usage (sparse + dense):
%   plot_cond_results(data_dir, sparse_file, sparse_breakdown_file, plot_mode, dense_file, dense_breakdown_file)
%
% plot_mode: 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%   'best_speed'  - for each (cond_num, algorithm), select the run with fastest time
%   'worst_ortho' - for each (cond_num, algorithm), select the run with worst orthogonality error
%   'best_ortho'  - for each (cond_num, algorithm), select the run with best orthogonality error
%
% Example (sparse only):
%   plot_cond_results('C:\Users\mymel\data\', ...
%       '20260122_143052_conditioning_sparse_results.csv', ...
%       '20260122_143052_conditioning_sparse_results_breakdown.csv', ...
%       'worst_ortho')
%
% Example (sparse + dense):
%   plot_cond_results('C:\Users\mymel\data\', ...
%       '20260122_143052_conditioning_sparse_results.csv', ...
%       '20260122_143052_conditioning_sparse_results_breakdown.csv', ...
%       'best_speed', ...
%       '20260122_143052_conditioning_dense_results.csv', ...
%       '20260122_143052_conditioning_dense_results_breakdown.csv')

    % Handle optional arguments
    if nargin < 4, plot_mode = 'best_speed'; end
    has_dense = (nargin >= 6) && ~isempty(dense_file) && ~isempty(dense_breakdown_file);

    % Mode display string for titles
    if strcmp(plot_mode, 'best_speed')
        mode_str = 'Best Speed';
    elseif strcmp(plot_mode, 'worst_ortho')
        mode_str = 'Worst Ortho';
    else
        mode_str = 'Best Ortho';
    end

    % Load comparison data (CQRRT\_linop vs CholQR vs sCholQR3)
    % Main results: 11 comment lines + 1 column header = 12 header lines
    % (includes prepended "# Total benchmark runtime: ..." line)
    data_sparse = readmatrix(fullfile(data_dir, sparse_file), 'NumHeaderLines', 12);
    % Breakdown: 15 comment lines + 1 column header = 16 header lines
    % (includes prepended "# Total benchmark runtime: ..." line
    %  and 4 lines describing CQRRT\_linop/CholQR/sCholQR3/CQRRT\_expl columns)
    breakdown_sparse = readmatrix(fullfile(data_dir, sparse_breakdown_file), 'NumHeaderLines', 16);

    if has_dense
        data_dense = readmatrix(fullfile(data_dir, dense_file), 'NumHeaderLines', 12);
        breakdown_dense = readmatrix(fullfile(data_dir, dense_breakdown_file), 'NumHeaderLines', 16);
    end

    %% Select data based on plot_mode
    % Results columns (26 total):
    % 1:cond_num, 2:run
    % CQRRT\_linop (3-6): orth_error, max_orth_cols, is_orth, time_us
    % CholQR (7-10): orth_error, max_orth_cols, is_orth, time_us
    % sCholQR3 (11-14): orth_error, max_orth_cols, is_orth, time_us
    % CQRRT\_expl (15-18): orth_error, max_orth_cols, is_orth, time_us
    % Memory (19-26): cqrrt_peak_rss_kb, cqrrt_analytical_kb,
    %   cholqr_peak_rss_kb, cholqr_analytical_kb,
    %   scholqr3_peak_rss_kb, scholqr3_analytical_kb,
    %   dense_cqrrt_peak_rss_kb, dense_cqrrt_analytical_kb

    algo_time_cols = [6, 10, 14, 18];   % CQRRT, CholQR, sCholQR3, CQRRT_expl
    algo_orth_cols = [3, 7, 11, 15];

    % Select sparse data
    [cond_num_sparse, sel_idx_sparse] = select_cond_data(data_sparse, plot_mode, algo_time_cols, algo_orth_cols);

    cqrrt_orth_sparse = data_sparse(sel_idx_sparse(:,1), 3);
    cqrrt_max_orth_cols_sparse = data_sparse(sel_idx_sparse(:,1), 4);
    cqrrt_time_sparse = data_sparse(sel_idx_sparse(:,1), 6) / 1e6;  % Convert microseconds to seconds
    cholqr_orth_sparse = data_sparse(sel_idx_sparse(:,2), 7);
    cholqr_max_orth_cols_sparse = data_sparse(sel_idx_sparse(:,2), 8);
    cholqr_time_sparse = data_sparse(sel_idx_sparse(:,2), 10) / 1e6;
    scholqr3_orth_sparse = data_sparse(sel_idx_sparse(:,3), 11);
    scholqr3_max_orth_cols_sparse = data_sparse(sel_idx_sparse(:,3), 12);
    scholqr3_time_sparse = data_sparse(sel_idx_sparse(:,3), 14) / 1e6;
    dense_cqrrt_orth_sparse = data_sparse(sel_idx_sparse(:,4), 15);
    dense_cqrrt_max_orth_cols_sparse = data_sparse(sel_idx_sparse(:,4), 16);
    dense_cqrrt_time_sparse = data_sparse(sel_idx_sparse(:,4), 18) / 1e6;

    % Memory (KB -> MB). Dimensions are constant across conditions and runs; use first row.
    first_sparse = find(data_sparse(:, 1) == cond_num_sparse(1));
    first_sparse = first_sparse(1);
    mem_peak_rss_sparse = [data_sparse(first_sparse,19), data_sparse(first_sparse,21), data_sparse(first_sparse,23), data_sparse(first_sparse,25)] / 1024;
    mem_analytical_sparse = [data_sparse(first_sparse,20), data_sparse(first_sparse,22), data_sparse(first_sparse,24), data_sparse(first_sparse,26)] / 1024;

    %% Extract columns - DENSE (if provided)
    if has_dense
        [cond_num_dense, sel_idx_dense] = select_cond_data(data_dense, plot_mode, algo_time_cols, algo_orth_cols);

        cqrrt_orth_dense = data_dense(sel_idx_dense(:,1), 3);
        cqrrt_max_orth_cols_dense = data_dense(sel_idx_dense(:,1), 4);
        cqrrt_time_dense = data_dense(sel_idx_dense(:,1), 6) / 1e6;
        cholqr_orth_dense = data_dense(sel_idx_dense(:,2), 7);
        cholqr_max_orth_cols_dense = data_dense(sel_idx_dense(:,2), 8);
        cholqr_time_dense = data_dense(sel_idx_dense(:,2), 10) / 1e6;
        scholqr3_orth_dense = data_dense(sel_idx_dense(:,3), 11);
        scholqr3_max_orth_cols_dense = data_dense(sel_idx_dense(:,3), 12);
        scholqr3_time_dense = data_dense(sel_idx_dense(:,3), 14) / 1e6;
        dense_cqrrt_orth_dense = data_dense(sel_idx_dense(:,4), 15);
        dense_cqrrt_max_orth_cols_dense = data_dense(sel_idx_dense(:,4), 16);
        dense_cqrrt_time_dense = data_dense(sel_idx_dense(:,4), 18) / 1e6;

        % Memory (KB -> MB)
        first_dense = find(data_dense(:, 1) == cond_num_dense(1));
        first_dense = first_dense(1);
        mem_peak_rss_dense = [data_dense(first_dense,19), data_dense(first_dense,21), data_dense(first_dense,23), data_dense(first_dense,25)] / 1024;
        mem_analytical_dense = [data_dense(first_dense,20), data_dense(first_dense,22), data_dense(first_dense,24), data_dense(first_dense,26)] / 1024;
    end

    %% Figure 1: Main Comparison
    if has_dense
        % 4x2 grid for sparse + dense
        figure('Position', [100, 100, 1400, 1200]);
        num_cols = 2;
    else
        % 4x1 grid for sparse only
        figure('Position', [100, 100, 700, 1200]);
        num_cols = 1;
    end

    %% SPARSE results (left column or single column)

    % SPARSE: Orthogonality error vs condition number
    subplot(4, num_cols, 1);
    loglog(cond_num_sparse, cqrrt_orth_sparse, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
    hold on;
    loglog(cond_num_sparse, cholqr_orth_sparse, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
    loglog(cond_num_sparse, scholqr3_orth_sparse, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
    loglog(cond_num_sparse, dense_cqrrt_orth_sparse, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
    xlabel('Condition Number \kappa', 'FontSize', 12);
    ylabel('Orthogonality Error ||Q''Q - I||_F / sqrt(n)', 'FontSize', 12);
    title(sprintf('Sparse: Orthogonality Error (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);

    % SPARSE: Max orthonormal columns
    subplot(4, num_cols, 1 + num_cols);
    n_cols = max(cqrrt_max_orth_cols_sparse);  % Get n from data
    loglog(cond_num_sparse, cqrrt_max_orth_cols_sparse, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
    hold on;
    loglog(cond_num_sparse, cholqr_max_orth_cols_sparse, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
    loglog(cond_num_sparse, scholqr3_max_orth_cols_sparse, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
    loglog(cond_num_sparse, dense_cqrrt_max_orth_cols_sparse, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
    xlims = xlim;
    loglog(xlims, [n_cols n_cols], 'k--', 'LineWidth', 1, 'DisplayName', sprintf('Perfect (%d/%d)', n_cols, n_cols));
    xlabel('Condition Number \kappa', 'FontSize', 12);
    ylabel(sprintf('Max Orthonormal Columns (out of %d)', n_cols), 'FontSize', 12);
    title(sprintf('Sparse: Orthonormal Columns (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    ylim([1 n_cols * 1.5]);
    set(gca, 'FontSize', 11);

    % SPARSE: Execution time
    subplot(4, num_cols, 1 + 2*num_cols);
    loglog(cond_num_sparse, cqrrt_time_sparse, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
    hold on;
    loglog(cond_num_sparse, cholqr_time_sparse, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
    loglog(cond_num_sparse, scholqr3_time_sparse, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
    loglog(cond_num_sparse, dense_cqrrt_time_sparse, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
    xlabel('Condition Number \kappa', 'FontSize', 12);
    ylabel('Execution Time (seconds)', 'FontSize', 12);
    title(sprintf('Sparse: Runtime (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);

    % SPARSE: Peak Working Memory
    subplot(4, num_cols, 1 + 3*num_cols);
    mem_data_sparse = [mem_peak_rss_sparse; mem_analytical_sparse]';  % 4x2 matrix
    b = bar(mem_data_sparse);
    b(1).FaceColor = [0.2 0.4 0.8];  % Blue for Peak RSS
    b(2).FaceColor = [1.0 0.5 0.0];  % Orange for Analytical
    set(gca, 'XTickLabel', {'CQRRT\_linop', 'CholQR', 'sCholQR3', 'CQRRT\_expl'});
    ylabel('Memory (MB)', 'FontSize', 12);
    title('Sparse: Peak Working Memory', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Peak RSS', 'Analytical', 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);

    %% Right column: DENSE results (only if provided)
    if has_dense
        % DENSE: Orthogonality error vs condition number
        subplot(4, 2, 2);
        loglog(cond_num_dense, cqrrt_orth_dense, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
        hold on;
        loglog(cond_num_dense, cholqr_orth_dense, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
        loglog(cond_num_dense, scholqr3_orth_dense, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
        loglog(cond_num_dense, dense_cqrrt_orth_dense, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
        xlabel('Condition Number \kappa', 'FontSize', 12);
        ylabel('Orthogonality Error ||Q''Q - I||_F / sqrt(n)', 'FontSize', 12);
        title(sprintf('Dense: Orthogonality Error (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 10);
        grid on;
        set(gca, 'FontSize', 11);

        % DENSE: Max orthonormal columns
        subplot(4, 2, 4);
        n_cols_dense = max(cqrrt_max_orth_cols_dense);
        loglog(cond_num_dense, cqrrt_max_orth_cols_dense, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
        hold on;
        loglog(cond_num_dense, cholqr_max_orth_cols_dense, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
        loglog(cond_num_dense, scholqr3_max_orth_cols_dense, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
        loglog(cond_num_dense, dense_cqrrt_max_orth_cols_dense, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
        xlims = xlim;
        loglog(xlims, [n_cols_dense n_cols_dense], 'k--', 'LineWidth', 1, 'DisplayName', sprintf('Perfect (%d/%d)', n_cols_dense, n_cols_dense));
        xlabel('Condition Number \kappa', 'FontSize', 12);
        ylabel(sprintf('Max Orthonormal Columns (out of %d)', n_cols_dense), 'FontSize', 12);
        title(sprintf('Dense: Orthonormal Columns (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 10);
        grid on;
        ylim([1 n_cols_dense * 1.5]);
        set(gca, 'FontSize', 11);

        % DENSE: Execution time
        subplot(4, 2, 6);
        loglog(cond_num_dense, cqrrt_time_dense, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_linop');
        hold on;
        loglog(cond_num_dense, cholqr_time_dense, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CholQR');
        loglog(cond_num_dense, scholqr3_time_dense, 'g-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'sCholQR3');
        loglog(cond_num_dense, dense_cqrrt_time_dense, 'm-d', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'CQRRT\_expl');
        xlabel('Condition Number \kappa', 'FontSize', 12);
        ylabel('Execution Time (seconds)', 'FontSize', 12);
        title(sprintf('Dense: Runtime (%s)', mode_str), 'FontSize', 14, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 10);
        grid on;
        set(gca, 'FontSize', 11);

        % DENSE: Peak Working Memory
        subplot(4, 2, 8);
        mem_data_dense = [mem_peak_rss_dense; mem_analytical_dense]';  % 4x2 matrix
        b = bar(mem_data_dense);
        b(1).FaceColor = [0.2 0.4 0.8];  % Blue for Peak RSS
        b(2).FaceColor = [1.0 0.5 0.0];  % Orange for Analytical
        set(gca, 'XTickLabel', {'CQRRT\_linop', 'CholQR', 'sCholQR3', 'CQRRT\_expl'});
        ylabel('Memory (MB)', 'FontSize', 12);
        title('Dense: Peak Working Memory', 'FontSize', 14, 'FontWeight', 'bold');
        legend('Peak RSS', 'Analytical', 'Location', 'best', 'FontSize', 10);
        grid on;
        set(gca, 'FontSize', 11);
    end

    %% Overall title for Figure 1
    if has_dense
        sgtitle(sprintf('CQRRT\\_linop vs CholQR vs sCholQR3 (%s) - Sparse (Left) vs Dense (Right)', mode_str), ...
            'FontSize', 16, 'FontWeight', 'bold');
    else
        sgtitle(sprintf('CQRRT\\_linop vs CholQR vs sCholQR3 (%s) - Sparse Input', mode_str), ...
            'FontSize', 16, 'FontWeight', 'bold');
    end

    %% Figure 2: Runtime Breakdown - Percentage (Horizontal Layout)
    % Breakdown columns (51 total):
    % 1:cond_num, 2:run
    % CQRRT\_linop (3-13): alloc, saso, qr, trtri, linop_precond, linop_gram, trmm_gram, potrf, finalize, rest, total
    % CholQR (14-19): alloc, materialize, gram, potrf, rest, total
    % sCholQR3 (20-32): alloc, materialize, gram1, potrf1, trsm1, syrk2, potrf2, update2, syrk3, potrf3, update3, rest, total
    % CQRRT\_expl (33-43): materialize, saso, qr, trtri(=0), precond, gram, trmm_gram(=0), potrf, finalize, rest, total
    % Memory (44-51): cqrrt_peak_rss_kb, cqrrt_analytical_kb, cholqr_peak_rss_kb, cholqr_analytical_kb,
    %   scholqr3_peak_rss_kb, scholqr3_analytical_kb, dense_cqrrt_peak_rss_kb, dense_cqrrt_analytical_kb

    % Select per-algorithm breakdown rows (sparse)
    cqrrt_bd_sparse = breakdown_sparse(sel_idx_sparse(:, 1), :);
    cholqr_bd_sparse = breakdown_sparse(sel_idx_sparse(:, 2), :);
    scholqr3_bd_sparse = breakdown_sparse(sel_idx_sparse(:, 3), :);
    dense_bd_sparse = breakdown_sparse(sel_idx_sparse(:, 4), :);

    if has_dense
        cqrrt_bd_dense = breakdown_dense(sel_idx_dense(:, 1), :);
        cholqr_bd_dense = breakdown_dense(sel_idx_dense(:, 2), :);
        scholqr3_bd_dense = breakdown_dense(sel_idx_dense(:, 3), :);
        dense_bd_dense = breakdown_dense(sel_idx_dense(:, 4), :);
    end

    if has_dense
        % 2x4 layout: rows = sparse/dense, columns = algorithms
        figure('Position', [150, 150, 2000, 700]);
        tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Row 1: Sparse (CQRRT, CholQR, sCholQR3, CQRRT\_expl)
        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_sparse, 'CQRRT\_linop - Sparse', true, 'pct', []);
        nexttile
        plot_cholqr_breakdown(cholqr_bd_sparse, 'CholQR - Sparse', false, 'pct', []);
        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_sparse, 'sCholQR3 - Sparse', false, 'pct', []);
        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_sparse, 'CQRRT\_expl - Sparse', false, 'pct', []);

        % Row 2: Dense (CQRRT, CholQR, sCholQR3, CQRRT\_expl)
        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_dense, 'CQRRT\_linop - Dense', true, 'pct', []);
        nexttile
        plot_cholqr_breakdown(cholqr_bd_dense, 'CholQR - Dense', false, 'pct', []);
        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_dense, 'sCholQR3 - Dense', false, 'pct', []);
        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_dense, 'CQRRT\_expl - Dense', false, 'pct', []);

        sgtitle(sprintf('Runtime Breakdown %% (%s) - Sparse (Top) vs Dense (Bottom)', mode_str), ...
            'FontSize', 14, 'FontWeight', 'bold');
    else
        % 1x4 horizontal layout: sparse only
        figure('Position', [150, 150, 2000, 450]);
        tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_sparse, 'CQRRT\_linop', true, 'pct', []);

        nexttile
        plot_cholqr_breakdown(cholqr_bd_sparse, 'CholQR', false, 'pct', []);

        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_sparse, 'sCholQR3', false, 'pct', []);

        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_sparse, 'CQRRT\_expl', false, 'pct', []);

        sgtitle(sprintf('Runtime Breakdown %% (%s) - Sparse Input', mode_str), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end

    %% Figure 3: Runtime Breakdown - Absolute Time (Horizontal Layout)
    if has_dense
        % Compute common y-axis limit across sparse and dense for each row
        abs_ylim_sparse = compute_common_ylim(cqrrt_bd_sparse, cholqr_bd_sparse, scholqr3_bd_sparse, dense_bd_sparse);
        abs_ylim_dense = compute_common_ylim(cqrrt_bd_dense, cholqr_bd_dense, scholqr3_bd_dense, dense_bd_dense);

        % 2x4 layout: rows = sparse/dense, columns = algorithms
        figure('Position', [200, 200, 2000, 700]);
        tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Row 1: Sparse (CQRRT, CholQR, sCholQR3, CQRRT\_expl)
        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_sparse, 'CQRRT\_linop - Sparse', true, 'abs', abs_ylim_sparse);
        nexttile
        plot_cholqr_breakdown(cholqr_bd_sparse, 'CholQR - Sparse', false, 'abs', abs_ylim_sparse);
        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_sparse, 'sCholQR3 - Sparse', false, 'abs', abs_ylim_sparse);
        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_sparse, 'CQRRT\_expl - Sparse', false, 'abs', abs_ylim_sparse);

        % Row 2: Dense (CQRRT, CholQR, sCholQR3, CQRRT\_expl)
        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_dense, 'CQRRT\_linop - Dense', true, 'abs', abs_ylim_dense);
        nexttile
        plot_cholqr_breakdown(cholqr_bd_dense, 'CholQR - Dense', false, 'abs', abs_ylim_dense);
        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_dense, 'sCholQR3 - Dense', false, 'abs', abs_ylim_dense);
        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_dense, 'CQRRT\_expl - Dense', false, 'abs', abs_ylim_dense);

        sgtitle(sprintf('Runtime Breakdown (seconds) (%s) - Sparse (Top) vs Dense (Bottom)', mode_str), ...
            'FontSize', 14, 'FontWeight', 'bold');
    else
        % Compute common y-axis limit for sparse only
        abs_ylim_sparse = compute_common_ylim(cqrrt_bd_sparse, cholqr_bd_sparse, scholqr3_bd_sparse, dense_bd_sparse);

        % 1x4 horizontal layout: sparse only
        figure('Position', [200, 200, 2000, 450]);
        tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile
        plot_cqrrt_breakdown(cqrrt_bd_sparse, 'CQRRT\_linop', true, 'abs', abs_ylim_sparse);

        nexttile
        plot_cholqr_breakdown(cholqr_bd_sparse, 'CholQR', false, 'abs', abs_ylim_sparse);

        nexttile
        plot_scholqr3_breakdown(scholqr3_bd_sparse, 'sCholQR3', false, 'abs', abs_ylim_sparse);

        nexttile
        plot_dense_cqrrt_breakdown(dense_bd_sparse, 'CQRRT\_expl', false, 'abs', abs_ylim_sparse);

        sgtitle(sprintf('Runtime Breakdown (seconds) (%s) - Sparse Input', mode_str), ...
            'FontSize', 14, 'FontWeight', 'bold');
    end

end

function [unique_conds, sel_idx] = select_cond_data(data, plot_mode, algo_time_cols, algo_orth_cols)
    % Select one row per (cond_num, algorithm) based on plot_mode
    % Returns unique condition numbers and selection indices into data
    unique_conds = unique(data(:, 1));
    num_conds = length(unique_conds);
    num_algos = length(algo_time_cols);
    sel_idx = zeros(num_conds, num_algos);
    for i = 1:num_conds
        mask = data(:, 1) == unique_conds(i);
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
end

function abs_ylim = compute_common_ylim(cqrrt_bd, cholqr_bd, scholqr3_bd, dense_bd)
    % Compute common y-axis limit across all four algorithms
    % Total time columns: CQRRT=13, CholQR=19, sCholQR3=32, CQRRT_expl=43
    cqrrt_total = cqrrt_bd(:, 13) / 1e6;   % Convert to seconds
    cholqr_total = cholqr_bd(:, 19) / 1e6;
    scholqr3_total = scholqr3_bd(:, 32) / 1e6;
    dense_cqrrt_total = dense_bd(:, 43) / 1e6;
    max_time = max([max(cqrrt_total), max(cholqr_total), max(scholqr3_total), max(dense_cqrrt_total)]);
    abs_ylim = [0, max_time * 1.05];  % Add 5% padding
end

function [] = plot_cqrrt_breakdown(breakdown_data, title_str, show_yaxis, mode, common_ylim)
    [Data, colors, labels] = get_cqrrt_breakdown(breakdown_data, mode);
    render_breakdown(Data, colors, labels, title_str, show_yaxis, mode, common_ylim, breakdown_data(:, 1));
end

function [] = plot_cholqr_breakdown(breakdown_data, title_str, show_yaxis, mode, common_ylim)
    [Data, colors, labels] = get_cholqr_breakdown(breakdown_data, mode);
    render_breakdown(Data, colors, labels, title_str, show_yaxis, mode, common_ylim, breakdown_data(:, 1), 7);
end

function [] = plot_scholqr3_breakdown(breakdown_data, title_str, show_yaxis, mode, common_ylim)
    [Data, colors, labels] = get_scholqr3_breakdown(breakdown_data, mode);
    render_breakdown(Data, colors, labels, title_str, show_yaxis, mode, common_ylim, breakdown_data(:, 1));
    xlabel('Condition Number \kappa', 'FontSize', 10);
end

function [] = plot_dense_cqrrt_breakdown(breakdown_data, title_str, show_yaxis, mode, common_ylim)
    [Data, colors, labels] = get_dense_cqrrt_breakdown(breakdown_data, mode);
    render_breakdown(Data, colors, labels, title_str, show_yaxis, mode, common_ylim, breakdown_data(:, 1));
end

function [Data, colors, labels] = get_cqrrt_breakdown(breakdown_data, mode)
    % CQRRT\_linop columns 3-13: alloc, saso, qr, trtri, linop_precond, linop_gram, trmm_gram, potrf, finalize, rest, total
    % Stack order (bottom to top): Other(rest+alloc), Fin, Chol, TRMM, LinOp(Gram), LinOp(Prec), TRTRI, QR, skop
    raw = [breakdown_data(:,12)+breakdown_data(:,3), breakdown_data(:, [11, 10, 9, 8, 7, 6, 5, 4])];
    total = breakdown_data(:, 13);
    if strcmp(mode, 'pct')
        Data = 100 * raw ./ total;
    else
        Data = raw / 1e6;
    end
    colors = {[0.6 0.6 0.6], [0.5 0.5 0.5], '#7E2F8E', '#FF7F0E', '#EDB120', '#B15928', '#E377C2', '#1F77B4', '#2CA02C'};
    labels = {'Other', 'Fin', 'Chol', 'TRMM', 'LinOp(Gram)', 'LinOp(Prec)', 'TRTRI', 'QR', 'skop'};
end

function [Data, colors, labels] = get_cholqr_breakdown(breakdown_data, mode)
    % CholQR columns 14-19: alloc, materialize, gram, potrf, rest, total
    % Stack order: Other(rest+alloc), Chol, Gram, LinOp
    raw = [breakdown_data(:,18)+breakdown_data(:,14), breakdown_data(:, [17, 16, 15])];
    total = breakdown_data(:, 19);
    if strcmp(mode, 'pct')
        Data = 100 * raw ./ total;
    else
        Data = raw / 1e6;
    end
    colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928'};
    labels = {'Other', 'Chol', 'Gram', 'LinOp'};
end

function [Data, colors, labels] = get_scholqr3_breakdown(breakdown_data, mode)
    % sCholQR3 columns 20-32: alloc, materialize, gram1, potrf1, trsm1, syrk2, potrf2, update2, syrk3, potrf3, update3, rest, total
    % Stack order: Other(rest+alloc), Upd_3, Chol_3, Syrk_3, Upd_2, Chol_2, Syrk_2, Trsm_1, Chol_1, Gram, LinOp
    raw = [breakdown_data(:,31)+breakdown_data(:,20), breakdown_data(:, [30, 29, 28, 27, 26, 25, 24, 23, 22, 21])];
    total = breakdown_data(:, 32);
    if strcmp(mode, 'pct')
        Data = 100 * raw ./ total;
    else
        Data = raw / 1e6;
    end
    colors = {[0.6 0.6 0.6], '#E31A1C', '#984EA3', '#FF7F0E', '#FB9A99', '#7E2F8E', '#FDBF6F', '#17BECF', '#CAB2D6', '#EDB120', '#B15928'};
    labels = {'Other', 'Upd_3', 'Chol_3', 'Syrk_3', 'Upd_2', 'Chol_2', 'Syrk_2', 'Trsm_1', 'Chol_1', 'Gram', 'LinOp'};
end

function [Data, colors, labels] = get_dense_cqrrt_breakdown(breakdown_data, mode)
    % CQRRT\_expl columns 33-43: materialize, saso, qr, trtri(=0), precond, gram, trmm_gram(=0), potrf, finalize, rest, total
    % Stack order: Fin+Other, Chol, Gram, Precond, QR, skop, Materialize
    fin_plus_rest = breakdown_data(:, 41) + breakdown_data(:, 42);
    raw = [fin_plus_rest, breakdown_data(:, [40, 38, 37, 35, 34, 33])];
    total = breakdown_data(:, 43);
    if strcmp(mode, 'pct')
        Data = 100 * raw ./ total;
    else
        Data = raw / 1e6;
    end
    colors = {[0.6 0.6 0.6], '#7E2F8E', '#EDB120', '#B15928', '#1F77B4', '#2CA02C', '#17BECF'};
    labels = {'Fin+Other', 'Chol', 'Gram', 'Precond', 'QR', 'skop', 'Materialize'};
end

function [] = render_breakdown(Data, colors, labels, title_str, show_yaxis, mode, common_ylim, cond_nums, lgd_fontsize)
    if nargin < 9, lgd_fontsize = 6; end
    bplot = bar(Data, 'stacked');
    for i = 1:length(colors)
        bplot(i).FaceColor = colors{i};
        bplot(i).FaceAlpha = 0.9;
    end
    if strcmp(mode, 'pct')
        ylim([0 100]);
    elseif ~isempty(common_ylim)
        ylim(common_ylim);
    end
    [label_indices, xtick_labels] = get_cond_xticks(cond_nums);
    set(gca, 'XTick', label_indices, 'XTickLabel', xtick_labels, 'FontSize', 9);
    xtickangle(45);
    title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
    if show_yaxis
        if strcmp(mode, 'pct')
            ylabel('Runtime %', 'FontSize', 10);
        else
            ylabel('Time (s)', 'FontSize', 10);
        end
    else
        set(gca, 'YTickLabel', []);
    end
    lgd = legend(labels{:}, 'Location', 'eastoutside');
    lgd.FontSize = lgd_fontsize;
end

function [label_indices, xtick_labels] = get_cond_xticks(cond_nums)
    % Create x-tick labels from condition numbers (log scale labels)
    num_conds = length(cond_nums);
    num_labels = min(8, num_conds);
    label_indices = round(linspace(1, num_conds, num_labels));

    xtick_labels = cell(1, length(label_indices));
    for i = 1:length(label_indices)
        idx = label_indices(i);
        % Format as 10^x for cleaner display
        exp_val = log10(cond_nums(idx));
        if abs(exp_val - round(exp_val)) < 0.1
            xtick_labels{i} = sprintf('10^{%d}', round(exp_val));
        else
            xtick_labels{i} = sprintf('%.1e', cond_nums(idx));
        end
    end
end

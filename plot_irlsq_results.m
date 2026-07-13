%% plot_irlsq_results.m — Plot CQRRT_linop_irlsq benchmark results
%
% The sparse-input IR-LSQ benchmark runs every Q-less QR variant selected by
% method_mask through the IR-LSQ pipeline (Algorithm 1, Epperly–Meier–Nakatsukasa
% 2025) on a tall sparse matrix loaded from a .mtx file. Each (algorithm, run)
% does Q-less QR → 2-step IR (x_0 = 0) with inner CG.
%
% CSV schema (results):
%   algorithm, run, m, n, qr_status, qr_time_us, peak_rss_kb, analytical_kb,
%   orth_error, ir_total_us, ir_outer_iters, ir_inner_iters_total,
%   ls_residual_norm, ls_solution_error
%
% ls_residual_norm uses the Higham normwise backward-error metric
%   ||A x - b|| / (||A||_2 * ||x|| + ||b||),
% drivable to machine epsilon for a backward-stable LS solver.
%
% ls_solution_error is -1 / NaN when no ground-truth x_true exists (FEM b = L^{-1} r).
%
% CSV schema (breakdown):
%   algorithm, run, phase, t0..t10
%   phase in {QR, IR}.  IR layout (6): outer_total, inner_cg_total, trsm, fwd, adj, other.
%
% Usage:
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv)
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix)
%   plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab)

function plot_irlsq_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab)

if nargin < 4, title_suffix = ''; end
if nargin < 5, main_tab = []; end
if nargin < 6, bd_tab = []; end

% =========================================================================
%  Wong colorblind-friendly palette (matches plot_applications_results.m)
% =========================================================================
w_blue      = [  0 114 178] / 255;
w_orange    = [230 159   0] / 255;
w_skyblue   = [ 86 180 233] / 255;
w_green     = [  0 158 115] / 255;
w_vermilion = [213  94   0] / 255;
w_purple    = [204 121 167] / 255;
w_gray      = [0.65 0.65 0.65];
w_ltgray    = [0.85 0.85 0.85];

% =========================================================================
%  Algorithm display order.  GEQP3-stabilized variant (CQRRT_linop_stb) is
%  intentionally absent: the benchmark dropped that path.
% =========================================================================
alg_csv_order  = {'CQRRT_linop', 'CholQR', 'CholQR2', 'sCholQR3_basic', 'sCholQR3', 'Blendenpik'};
alg_disp_names = {'CQRRT\_linop', 'CholQR', 'CholQR2', 'sCholQR3', 'sCholQR3\_se', 'Blendenpik'};

% =========================================================================
%  Load CSVs.  count_comment_lines is a helper in this directory; we rely on
%  it to skip the leading "#"-prefixed metadata block.
% =========================================================================
results_path   = fullfile(data_dir, results_csv);
breakdown_path = fullfile(data_dir, breakdown_csv);
if ~isfile(results_path)
    error('plot_irlsq_results: file not found: %s', results_path);
end
n_skip = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_skip);
opts = setvartype(opts, 'algorithm', 'char');
T = readtable(results_path, opts);

algorithms       = T.algorithm;
qr_status        = T.qr_status;
qr_time          = T.qr_time_us;
peak_rss_kb      = T.peak_rss_kb;
analytical_kb    = T.analytical_kb;
ir_total_us      = T.ir_total_us;
ir_inner_total   = T.ir_inner_iters_total;
ls_residual_norm = T.ls_residual_norm;
orth_error       = T.orth_error;
m_val            = T.m(1);
n_val            = T.n(1);
% chol_retries (adaptive-shift retries) is present in irlsq_reg CSVs from 2026-07;
% default to 0 for older CSVs that predate the column.
if ismember('chol_retries', T.Properties.VariableNames)
    chol_retries = T.chol_retries;
else
    chol_retries = zeros(height(T), 1);
end

% =========================================================================
%  Sort algorithms by display order
% =========================================================================
unique_algs_csv = unique(algorithms, 'stable');
unique_algs = {};
disp_labels = {};
for k = 1:numel(alg_csv_order)
    if any(strcmp(unique_algs_csv, alg_csv_order{k}))
        unique_algs{end+1} = alg_csv_order{k}; %#ok<AGROW>
        disp_labels{end+1} = alg_disp_names{k}; %#ok<AGROW>
    end
end
% Catch any extras (defensive)
for k = 1:numel(unique_algs_csv)
    if ~any(strcmp(unique_algs, unique_algs_csv{k}))
        unique_algs{end+1} = unique_algs_csv{k}; %#ok<AGROW>
        disp_labels{end+1} = strrep(unique_algs_csv{k}, '_', '\_'); %#ok<AGROW>
    end
end
n_algs = numel(unique_algs);

% =========================================================================
%  Per-algorithm: pick the run with smallest total wall time (QR + IR).
%  If all runs of an algorithm failed (qr_status != 0), keep the first row
%  but flag it visually so failures aren't silently dropped.
% =========================================================================
sel_idx     = zeros(n_algs, 1);
sel_failed  = false(n_algs, 1);
for a = 1:n_algs
    mask     = strcmp(algorithms, unique_algs{a});
    indices  = find(mask);
    succ     = qr_status(mask) == 0;
    if any(succ)
        succ_idx = indices(succ);
        total = qr_time(succ_idx) + ir_total_us(succ_idx);
        [~, best] = min(total);
        sel_idx(a) = succ_idx(best);
    else
        sel_idx(a) = indices(1);
        sel_failed(a) = true;
    end
end

% =========================================================================
%  Main figure: 2x3 layout (5 panels; bottom-right tile blank)
%    (1) Stacked timing: QR + IR
%    (2) Residual: Higham normwise backward error
%    (3) Memory: peak RSS vs analytical prediction
%    (4) Inner CG iterations (total)
%    (5) Q-factor orthogonality loss ||Q^T Q - I||_F / sqrt(n)
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 1100, 900]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
tl_main = tiledlayout(parent_main, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title_label = sprintf('Sparse IR-LSQ Benchmark — %d \\times %d', m_val, n_val);
if ~isempty(title_suffix)
    title_label = sprintf('%s — %s', title_label, title_suffix);
end
title(tl_main, title_label, 'FontWeight', 'bold');

x_pos = 1:n_algs;

% ---- (1) Stacked timing: QR + IR ----
nexttile(tl_main);
qr_ms  = arrayfun(@(i) qr_time(i)/1000,     sel_idx);
ir_ms  = arrayfun(@(i) ir_total_us(i)/1000, sel_idx);
qr_ms(sel_failed) = 0;  ir_ms(sel_failed) = 0;
b = bar(x_pos, [qr_ms, ir_ms], 'stacked');
b(1).FaceColor = w_blue;     b(1).DisplayName = 'QR';
b(2).FaceColor = w_orange;   b(2).DisplayName = 'IR-LSQ (x_0 = 0)';
ylabel('Time (ms)'); title('Wall-time per algorithm');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;
for a = 1:n_algs
    if sel_failed(a)
        text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (2) Residual: Higham normwise backward error ----
nexttile(tl_main);
resid = arrayfun(@(i) ls_residual_norm(i), sel_idx);
resid(sel_failed) = NaN;
resid(resid < 0) = NaN;
bar(x_pos, resid, 'FaceColor', w_skyblue); set(gca, 'YScale', 'log');
ylim([1e-16, 1e0]);
ylabel('||Ax - b|| / (||A||\cdot||x|| + ||b||)'); title('Normwise backward error');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(resid(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (3) Memory: peak RSS vs analytical prediction ----
nexttile(tl_main);
mem_peak = arrayfun(@(i) peak_rss_kb(i),   sel_idx) / 1024;     % MB
mem_pred = arrayfun(@(i) analytical_kb(i), sel_idx) / 1024;     % MB
mem_peak(sel_failed) = 0;  mem_pred(sel_failed) = 0;
b = bar(x_pos, [mem_peak, mem_pred], 'grouped');
b(1).FaceColor = w_purple;     b(1).DisplayName = 'Peak RSS';
b(2).FaceColor = w_ltgray;     b(2).DisplayName = 'Analytical';
ylabel('Memory (MB)'); title('Peak vs predicted working memory');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;

% ---- (4) Inner CG iterations (total across both outer IR steps) ----
% Algorithmic signal: lower = R is a better preconditioner for A^T A.
% Outer iters are fixed at n_refine_steps = 2 by construction; only inner CG
% varies across algorithms, so we plot only inner totals.
nexttile(tl_main);
inner_iters = arrayfun(@(i) ir_inner_total(i), sel_idx);
inner_iters(sel_failed) = NaN;
bar(x_pos, inner_iters, 'FaceColor', w_orange);
ylabel('Inner CG iterations (total)');
title('Inner CG iterations to convergence');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(inner_iters(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    else
        % Annotate the actual integer count on top of each bar.
        text(x_pos(a), inner_iters(a), sprintf('%d', inner_iters(a)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontWeight', 'bold');
    end
end

% ---- (5) Orthogonality loss in Q-factor: ||Q^T Q - I||_F / sqrt(n) ----
nexttile(tl_main);
orth_vals = arrayfun(@(i) orth_error(i), sel_idx);
orth_vals(sel_failed) = NaN;
orth_vals(orth_vals < 0) = NaN;
bar(x_pos, orth_vals, 'FaceColor', w_vermilion); set(gca, 'YScale', 'log');
ylim([1e-16, 1e0]);
ylabel('||Q^T Q - I||_F / \surd n'); title('Q-factor orthogonality loss');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(orth_vals(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (6) Cholesky adaptive-shift retries (0 = clean; N/A for Blendenpik) ----
% How many times each CholeskyQR method had to grow the diagonal shift and retry
% potrf. 0 = the unshifted first attempt succeeded; higher = a more ill-conditioned
% (e.g. single-precision) Gram that needed regularization to factor.
nexttile(tl_main);
retries = arrayfun(@(i) chol_retries(i), sel_idx);
bar(x_pos, retries, 'FaceColor', w_gray);
ylabel('Cholesky shift retries'); title('Adaptive-shift retries');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
ylim([0, max(1, max(retries) * 1.25 + 1)]);
for a = 1:n_algs
    if strcmp(unique_algs{a}, 'Blendenpik')
        lbl = 'N/A';   % Blendenpik is not a CholeskyQR method
    else
        lbl = sprintf('%d', retries(a));
    end
    text(x_pos(a), retries(a), lbl, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

% =========================================================================
%  Breakdown figure (optional): IR-LSQ phase breakdown stacked bar
%  Layout (6 fields): outer_total, inner_cg_total, trsm, fwd, adj, other
% =========================================================================
if isfile(breakdown_path)
    n_skip_b = count_comment_lines(breakdown_path);
    Tb = readtable(breakdown_path, 'NumHeaderLines', n_skip_b);

    if isempty(bd_tab)
        figure('Position', [200, 100, 900, 500]);
        parent_bd = gcf;
    else
        parent_bd = bd_tab;
    end
    tl_bd = tiledlayout(parent_bd, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl_bd, sprintf('IR-LSQ runtime breakdown — %s', title_label), 'FontWeight', 'bold');

    % Build a per-algorithm matrix [inner_cg | trsm | fwd | adj | other] in ms.
    % All 5 segments come from IterRefineLSQ's populate_times() and sum to
    % outer_total = t0. (x_0 = 0 now: there is no sketch-and-solve initial guess,
    % so ir_total_us matches outer_total up to loop overhead — no x0 segment.)
    ir_breakdown = zeros(n_algs, 5);
    for a = 1:n_algs
        if sel_failed(a), continue; end
        run_idx = T.run(sel_idx(a));
        match = strcmp(Tb.algorithm, unique_algs{a}) & Tb.run == run_idx & strcmp(Tb.phase, 'IR');
        idx = find(match, 1);
        if isempty(idx), continue; end
        ir_breakdown(a, 1) = Tb.t1(idx)  / 1000;  % inner_cg ex. fwd/adj/trsm
        ir_breakdown(a, 2) = Tb.t2(idx)  / 1000;  % trsm
        ir_breakdown(a, 3) = Tb.t3(idx)  / 1000;  % fwd
        ir_breakdown(a, 4) = Tb.t4(idx)  / 1000;  % adj
        ir_breakdown(a, 5) = Tb.t5(idx)  / 1000;  % other
    end

    nexttile(tl_bd);
    b = bar(x_pos, ir_breakdown, 'stacked');
    b(1).FaceColor = w_orange;    b(1).DisplayName = 'inner CG control (axpy/dot)';
    b(2).FaceColor = w_skyblue;   b(2).DisplayName = 'TRSM';
    b(3).FaceColor = w_green;     b(3).DisplayName = 'J fwd';
    b(4).FaceColor = w_vermilion; b(4).DisplayName = 'J^T adj';
    b(5).FaceColor = w_gray;      b(5).DisplayName = 'axpy/copy/nrm2';
    ylabel('Time (ms)'); title('IR-LSQ phase breakdown');
    xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
    legend('Location', 'northeastoutside'); grid on; box on;
    for a = 1:n_algs
        if sel_failed(a)
            text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', w_vermilion);
        end
    end
end

end

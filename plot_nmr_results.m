%% plot_nmr_results.m — Plot CQRRT_linop_nmr benchmark results
%
% The NMR benchmark runs every Q-less QR variant selected by method_mask through
% the IR-LSQ pipeline (Algorithm 1, Epperly–Meier–Nakatsukasa 2025) on the
% Kronecker-structured 2D NMR relaxometry operator A = A2 ⊗ A1.
%
% CSV schema (results):
%   algorithm, run, m, n, qr_status, qr_time_us, peak_rss_kb, analytical_kb,
%   ir_total_us, ir_outer_iters, ir_inner_iters_total, ls_residual_norm, ls_solution_error
%
% CSV schema (breakdown):
%   algorithm, run, phase, t0..t10
%   phase ∈ {QR, IR}.  IR layout (6): outer_total, inner_cg_total, trsm, fwd, adj, other.
%
% Usage:
%   plot_nmr_results(data_dir, results_csv, breakdown_csv)
%   plot_nmr_results(data_dir, results_csv, breakdown_csv, title_suffix)
%   plot_nmr_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab)

function plot_nmr_results(data_dir, results_csv, breakdown_csv, title_suffix, main_tab, bd_tab)

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
%  Algorithm display order — same convention as plot_applications_results.m
% =========================================================================
alg_csv_order  = {'CQRRT_linop', 'CQRRT_linop_stb', 'CQRRT_linop_stb_bqrrp', ...
                  'CholQR', 'sCholQR3_basic', 'sCholQR3'};
alg_disp_names = {'CQRRT\_linop', 'CQRRT\_linop\_stb', 'CQRRT\_linop\_stb\_bqrrp', ...
                  'CholQR', 'sCholQR3', 'sCholQR3\_se'};

% =========================================================================
%  Load CSVs.  count_comment_lines is a helper in this directory; we rely on
%  it to skip the leading "#"-prefixed metadata block.
% =========================================================================
results_path   = fullfile(data_dir, results_csv);
breakdown_path = fullfile(data_dir, breakdown_csv);
if ~isfile(results_path)
    error('plot_nmr_results: file not found: %s', results_path);
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
ir_outer         = T.ir_outer_iters;
ir_inner_total   = T.ir_inner_iters_total;
ls_residual_norm = T.ls_residual_norm;
ls_solution_err  = T.ls_solution_error;
m_val            = T.m(1);
n_val            = T.n(1);

% =========================================================================
%  Sort algorithms by display order
% =========================================================================
unique_algs_csv = unique(algorithms, 'stable');
unique_algs = {};
disp_labels = {};
for k = 1:numel(alg_csv_order)
    if any(strcmp(unique_algs_csv, alg_csv_order{k}))
        unique_algs{end+1} = alg_csv_order{k};
        disp_labels{end+1} = alg_disp_names{k};
    end
end
% Catch any extras (defensive: an algorithm we didn't expect)
for k = 1:numel(unique_algs_csv)
    if ~any(strcmp(unique_algs, unique_algs_csv{k}))
        unique_algs{end+1} = unique_algs_csv{k};
        disp_labels{end+1} = strrep(unique_algs_csv{k}, '_', '\_');
    end
end
n_algs = numel(unique_algs);

% =========================================================================
%  Per-algorithm: pick the run with smallest total wall time (QR + IR).
%  If all runs of an algorithm failed (qr_status != 0), keep the first row
%  but flag it visually so failures aren't silently dropped from the chart.
% =========================================================================
sel_idx     = zeros(n_algs, 1);
sel_failed  = false(n_algs, 1);
for a = 1:n_algs
    mask     = strcmp(algorithms, unique_algs{a});
    indices  = find(mask);
    succ     = qr_status(mask) == 0;
    if any(succ)
        % Pick fastest successful run
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
%  Main figure: 2x2 layout
%    (1) Stacked timing: QR + IR
%    (2) Solution error  ||x - x_true|| / ||x_true||
%    (3) Residual error  ||b - A x|| / ||b||
%    (4) Memory: peak RSS vs analytical prediction
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 1100, 900]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
tl_main = tiledlayout(parent_main, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

title_label = sprintf('NMR Benchmark — %d × %d', m_val, n_val);
if ~isempty(title_suffix)
    title_label = sprintf('%s — %s', title_label, title_suffix);
end
title(tl_main, title_label, 'FontWeight', 'bold');

x_pos = 1:n_algs;
fail_text_y_offset = 1.05;  % vertical placement of FAIL annotations (relative)

% ---- (1) Stacked timing: QR + IR ----
nexttile(tl_main);
qr_ms  = arrayfun(@(i) qr_time(i)/1000,         sel_idx);
ir_ms  = arrayfun(@(i) ir_total_us(i)/1000,     sel_idx);
qr_ms(sel_failed) = 0;  ir_ms(sel_failed) = 0;
b = bar(x_pos, [qr_ms, ir_ms], 'stacked');
b(1).FaceColor = w_blue;     b(1).DisplayName = 'QR';
b(2).FaceColor = w_orange;   b(2).DisplayName = 'IR-LSQ';
ylabel('Time (ms)'); title('Wall-time per algorithm');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;
% Mark failures
for a = 1:n_algs
    if sel_failed(a)
        text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (2) Solution error ----
nexttile(tl_main);
sol_err = arrayfun(@(i) ls_solution_err(i), sel_idx);
sol_err(sel_failed) = NaN;
% Replace any sentinel (-1) with NaN for plotting
sol_err(sol_err < 0) = NaN;
bar(x_pos, sol_err, 'FaceColor', w_green); set(gca, 'YScale', 'log');
ylabel('||x - x_{true}|| / ||x_{true}||'); title('Solution error');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;  % preserve auto-y-range before annotating failures
for a = 1:n_algs
    if sel_failed(a) || isnan(sol_err(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (3) Residual ||b - A x|| / ||b|| ----
nexttile(tl_main);
resid = arrayfun(@(i) ls_residual_norm(i), sel_idx);
resid(sel_failed) = NaN;
resid(resid < 0) = NaN;
bar(x_pos, resid, 'FaceColor', w_skyblue); set(gca, 'YScale', 'log');
ylabel('||b - A x|| / ||b||'); title('Relative residual');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
yl = ylim;
for a = 1:n_algs
    if sel_failed(a) || isnan(resid(a))
        text(x_pos(a), yl(2)*0.9, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (4) Memory: peak RSS vs analytical prediction ----
nexttile(tl_main);
mem_peak = arrayfun(@(i) peak_rss_kb(i),    sel_idx) / 1024;     % MB
mem_pred = arrayfun(@(i) analytical_kb(i),  sel_idx) / 1024;     % MB
mem_peak(sel_failed) = 0;  mem_pred(sel_failed) = 0;
b = bar(x_pos, [mem_peak, mem_pred], 'grouped');
b(1).FaceColor = w_purple;     b(1).DisplayName = 'Peak RSS';
b(2).FaceColor = w_ltgray;     b(2).DisplayName = 'Analytical';
ylabel('Memory (MB)'); title('Peak vs predicted working memory');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;

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

    % Build a per-algorithm matrix [inner_cg | trsm | fwd | adj | other] in ms
    % (skipping outer_total = sum of the others, to avoid double counting)
    ir_breakdown = zeros(n_algs, 5);
    for a = 1:n_algs
        if sel_failed(a), continue; end
        run_idx = T.run(sel_idx(a));
        match = strcmp(Tb.algorithm, unique_algs{a}) & Tb.run == run_idx & strcmp(Tb.phase, 'IR');
        idx = find(match, 1);
        if isempty(idx), continue; end
        ir_breakdown(a, 1) = Tb.t1(idx) / 1000;  % inner_cg
        ir_breakdown(a, 2) = Tb.t2(idx) / 1000;  % trsm
        ir_breakdown(a, 3) = Tb.t3(idx) / 1000;  % fwd
        ir_breakdown(a, 4) = Tb.t4(idx) / 1000;  % adj
        ir_breakdown(a, 5) = Tb.t5(idx) / 1000;  % other
    end

    nexttile(tl_bd);
    b = bar(x_pos, ir_breakdown, 'stacked');
    b(1).FaceColor = w_orange;    b(1).DisplayName = 'inner CG (ex. fwd/adj/trsm)';
    b(2).FaceColor = w_skyblue;   b(2).DisplayName = 'TRSM';
    b(3).FaceColor = w_green;     b(3).DisplayName = 'J fwd';
    b(4).FaceColor = w_vermilion; b(4).DisplayName = 'J^T adj';
    b(5).FaceColor = w_gray;      b(5).DisplayName = 'other';
    ylabel('Time (ms)'); title('IR-LSQ phase breakdown');
    xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
    legend('Location', 'northeastoutside'); grid on; box on;
    % Mark failures
    for a = 1:n_algs
        if sel_failed(a)
            text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', w_vermilion);
        end
    end
end

end

%% plot_toeplitz_results.m — plot the Toeplitz LS benchmark (Oleg's 2nd experiment).
%
% One row per method in the CSV (num_runs=1). 6-panel layout matching the App-1
% plotter: wall-time (build+solve), accuracy (data + recovery rel error), memory
% (peak vs analytical), solver iterations, Q-orthogonality loss, Cholesky retries.
%
% CSV schema (see toeplitz_ls_benchmark.cc):
%   algorithm,m,n,qr_status,qr_time_us,solve_time_us,peak_rss_kb,analytical_kb,
%   orth_error,iterations,solver_flag,solver_relres,aug_relres,normal_relres,
%   data_relres,recovery_error,cond_estimate,chol_retries
% (solver_relres / normal_relres / cond_estimate were added to the benchmark on
%  2026-07-15 but omitted from this comment until 2026-07-27; readtable is name-based,
%  so the stale comment was misleading rather than breaking.)
%
% Usage:
%   plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab)
%   plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab, timing_agg)
%
% timing_agg ('best' default | 'mean') controls how multi-run CSVs (num_runs > 1,
% 2026-08-05) collapse to one row per method; see aggregate_runs.m.

function plot_toeplitz_results(data_dir, results_csv, title_suffix, main_tab, timing_agg)

if nargin < 3, title_suffix = ''; end
if nargin < 4, main_tab = []; end
if nargin < 5 || isempty(timing_agg), timing_agg = 'best'; end

% Wong colorblind-friendly palette (matches the App-1 plotter).
w_blue = [0 114 178]/255;  w_orange = [230 159 0]/255;  w_skyblue = [86 180 233]/255;
w_green = [0 158 115]/255; w_vermilion = [213 94 0]/255; w_purple = [204 121 167]/255;
w_gray = [0.65 0.65 0.65]; w_ltgray = [0.85 0.85 0.85];

% "Blendenpik_cold" appears in [NEW 08-05]+ CSVs (warm x_0 is Blendenpik-only and
% both variants run). "Blendenpik" keeps its bare label because its x_0 policy is
% era-dependent (old cold campaigns disabled its warm start); the era note carries it.
alg_csv_order  = {'CQRRT_linop','CholQR','CholQR2','sCholQR3_basic','sCholQR3','Blendenpik','Blendenpik_cold','unpreconditioned'};
alg_disp_names = {'CQRRT\_linop','CholQR','CholQR2','sCholQR3 (basic)','sCholQR3 (blocked)','Blendenpik','Blendenpik (cold x_0)','unprec'};

results_path = fullfile(data_dir, results_csv);
if ~isfile(results_path), error('plot_toeplitz_results: file not found: %s', results_path); end
n_skip = count_comment_lines(results_path);
opts = detectImportOptions(results_path, 'NumHeaderLines', n_skip);
opts = setvartype(opts, 'algorithm', 'char');
T = readtable(results_path, opts);

% Multi-run CSVs (num_runs > 1, 2026-08-05): collapse to one row per method
% before anything indexes the table. agg_note lands in the figure title.
[T, agg_note] = aggregate_runs(T, timing_agg);

algs = T.algorithm;
m_val = T.m(1); n_val = T.n(1);

% Order rows by display order (skip methods not present).
unique_algs = {}; disp_labels = {}; sel = [];
for k = 1:numel(alg_csv_order)
    idx = find(strcmp(algs, alg_csv_order{k}), 1);
    if ~isempty(idx)
        unique_algs{end+1} = alg_csv_order{k}; %#ok<AGROW>
        disp_labels{end+1} = alg_disp_names{k}; %#ok<AGROW>
        sel(end+1) = idx; %#ok<AGROW>
    end
end
sel = sel(:);   % column, so metric vectors are Nx1 and [a,b] is Nx2 for grouped/stacked bars
na = numel(sel);
x = 1:na;
failed = (T.qr_status(sel) ~= 0);

% --- Figure / tab ---
if isempty(main_tab), figure('Position',[100 100 1200 900]); parent = gcf; else parent = main_tab; end
tl = tiledlayout(parent, 2, 4, 'TileSpacing','compact', 'Padding','compact');
ttl = sprintf('Toeplitz LS Benchmark — %d \\times %d', m_val, n_val);
if ~isempty(title_suffix), ttl = sprintf('%s — %s', ttl, title_suffix); end
if ~isempty(agg_note), ttl = sprintf('%s — %s', ttl, agg_note); end
title(tl, ttl, 'FontWeight','bold');

% (1) Wall-time: build (QR/sketch) + warm start + solve (LSQR), stacked, ms.
% The warm-start segment (setup_us column, 2026-07-31) is the auxiliary
% sketch-and-solve x0 build; orange is reserved for it figure-wide. CSVs
% without the column (all campaigns before 07-31) render as before.
nexttile(tl, 1);
build_ms = T.qr_time_us(sel)/1000;  solve_ms = T.solve_time_us(sel)/1000;
ws_ms = zeros(na, 1);
if ismember('setup_us', T.Properties.VariableNames)
    ws_ms = max(T.setup_us(sel), 0)/1000;
end
build_ms(failed) = 0; solve_ms(failed) = 0; ws_ms(failed) = 0;
b = bar(x, [build_ms, ws_ms, solve_ms], 'stacked');
b(1).FaceColor = w_blue;      b(1).DisplayName = 'QR (Q-less)';
b(2).FaceColor = w_vermilion; b(2).DisplayName = 'warm start (x_0 build)';
b(3).FaceColor = w_orange;    b(3).DisplayName = 'LSQR solve';
if ~any(ws_ms > 0), delete(b(2)); end
% Blendenpik's warm start reuses its own sketch factors (setup_us = 0 by design),
% so it never gets an orange segment -- which reads as "not warm-started"
% (Max, 2026-07-31). Say it explicitly when the CSV provenance header says so.
hdr = fileread(results_path);
if ~isempty(regexp(hdr, 'bp_warm_start=1', 'once'))
    bp_i = find(strcmp(unique_algs, 'Blendenpik'), 1);
    if ~isempty(bp_i) && ~failed(bp_i)
        text(x(bp_i), build_ms(bp_i) + ws_ms(bp_i) + solve_ms(bp_i), 'warm x_0 (embedded)', ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 8, 'Color', w_vermilion, 'FontWeight', 'bold');
    end
end
ylabel('Time (ms)'); title('Wall-time per algorithm');   % linear scale, starts at 0
ylim([0 max(1, 1.15*max(build_ms + ws_ms + solve_ms))]);
xticks(x); xticklabels(disp_labels); xtickangle(35); legend('Location','northeast'); grid on; box on;

% (2) Accuracy: data rel error (||Tx-b||/||b||) + recovery error (||x-xtrue||/||xtrue||).
nexttile(tl, 2);
dre = T.data_relres(sel); rec = T.recovery_error(sel);
dre(failed) = NaN; rec(failed) = NaN; dre(dre<0) = NaN; rec(rec<0) = NaN;
b = bar(x, [dre, rec], 'grouped'); set(gca,'YScale','log');
b(1).FaceColor = w_skyblue; b(1).DisplayName = '||Tx-b||/||b||';
b(2).FaceColor = w_green;   b(2).DisplayName = '||x-x_{true}||/||x_{true}||';
ylabel('relative error'); title('Data & recovery error');
xticks(x); xticklabels(disp_labels); xtickangle(35); legend('Location','northwest'); grid on; box on;

% (3) Memory: peak RSS vs analytical (MB).
nexttile(tl, 3);
peak_mb = T.peak_rss_kb(sel)/1024; ana_mb = T.analytical_kb(sel)/1024;
peak_mb(failed) = 0; ana_mb(failed) = 0;
b = bar(x, [peak_mb, ana_mb], 'grouped');
b(1).FaceColor = w_purple; b(1).DisplayName = 'Peak RSS';
b(2).FaceColor = w_ltgray; b(2).DisplayName = 'Analytical';
ylabel('Memory (MB)'); title('Peak vs predicted memory');
xticks(x); xticklabels(disp_labels); xtickangle(35); legend('Location','northwest'); grid on; box on;

% (4) Solver iterations (LSQR / CG), with the warm-start x0 build shown as a
% base segment in TIME-EQUIVALENT iterations (setup time / per-iteration solve
% time), so the stacked total stays proportional to the wall-clock bar.
nexttile(tl, 5);
iters = T.iterations(sel); iters(failed) = NaN;
setup_iters = zeros(na, 1);
if ismember('setup_us', T.Properties.VariableNames)
    su = max(T.setup_us(sel), 0); sv = max(T.solve_time_us(sel), 1);
    setup_iters = su ./ sv .* iters;
    setup_iters(~isfinite(setup_iters)) = 0;
end
hb = bar(x, [setup_iters, iters], 'stacked');
hb(1).FaceColor = w_vermilion; hb(1).DisplayName = 'warm start (time-equiv iters)';
hb(2).FaceColor = w_orange;    hb(2).DisplayName = 'LSQR';
if any(setup_iters > 0), legend('Location', 'northwest'); end
ylabel('solver iterations'); title('LSQR iterations to convergence');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
for a = 1:na
    if failed(a), continue; end
    text(x(a), setup_iters(a) + iters(a), sprintf('%d', iters(a)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
end

% (4b, tile 6) QR build: canonical GEQRF rate. KB / dissertation convention:
% flops of standard Householder QR for this (m,n) over the tested method's
% measured build time -- method-independent numerator, so different QRs compare
% fairly. unprec has no QR build => N/A.
nexttile(tl, 6);
fl_canon   = 2*m_val*n_val^2 - (2/3)*n_val^3;
qr_s       = T.qr_time_us(sel) * 1e-6;
canon_rate = (fl_canon ./ qr_s) / 1e9;
canon_rate(qr_s <= 0 | failed) = NaN;
bar(x, canon_rate, 'FaceColor', w_blue);
ylabel('canonical GFLOP/s'); title('QR build: canonical GEQRF rate');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
for a = 1:na
    if isnan(canon_rate(a))
        text(x(a), 0, 'N/A', 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','Color',w_gray);
    else
        text(x(a), canon_rate(a), sprintf('%.0f', canon_rate(a)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
    end
end

% (5) Q-orthogonality loss ||R^-T A'A R^-1 - I||_F/sqrt(n). -1 => no preconditioner.
nexttile(tl, 4);
orth = T.orth_error(sel); orth(failed) = NaN; orth(orth<0) = NaN;
% w_green, NOT vermilion: orange is reserved for warm-start segments (2026-07-31).
bar(x, orth, 'FaceColor', w_green); set(gca,'YScale','log'); ylim([1e-16 1e1]);
ylabel('||R^{-T}A^TA R^{-1}-I||_F/\surd n'); title('Preconditioner orthogonality loss');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
yl = ylim;
for a = 1:na
    if strcmp(unique_algs{a}, 'unpreconditioned') || isnan(orth(a))
        text(x(a), yl(2)*0.5, 'N/A', 'HorizontalAlignment','center','FontWeight','bold','Color',w_gray);
    end
end

% (6) Cholesky adaptive-shift retries (N/A for Blendenpik/unpreconditioned).
nexttile(tl, 7);
retr = T.chol_retries(sel);
bar(x, retr, 'FaceColor', w_gray);
ylabel('Cholesky shift retries'); title('Adaptive-shift retries');
xticks(x); xticklabels(disp_labels); xtickangle(35); grid on; box on;
ylim([0, max(1, max(retr)*1.25 + 1)]);
for a = 1:na
    if startsWith(unique_algs{a}, 'Blendenpik') || strcmp(unique_algs{a}, 'unpreconditioned')
        lbl = 'N/A';
    else
        lbl = sprintf('%d', retr(a));
    end
    text(x(a), retr(a), lbl, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
end

end

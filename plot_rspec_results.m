function plot_rspec_results(data_dir, results_csv, title_suffix, main_tab)
% plot_rspec_results — visualize one RSPEC (Algorithm 4, reduced spectral
% approximation) CSV from CQRRT_linop_applications.
%
% CSV schema (post-2026-06-05 binary):
%   algorithm, run, m, n, omega, power_j, qr_status, qr_time_us,
%   peak_rss_kb, analytical_kb, factor_time_us, rspec_total_us,
%   eig_0 .. eig_9, resid_0 .. resid_9
%
% Layout (2x3): wall-time stacked (PCholQR + post-processing),
%   memory (peak vs analytical), Ritz residuals (log scale, grouped by
%   eigenvalue rank), top-k eigenvalues (signed, grouped), eigenvalue
%   agreement table, blank.
%
% Usage:
%   plot_rspec_results(data_dir, results_csv)
%   plot_rspec_results(data_dir, results_csv, title_suffix)
%   plot_rspec_results(data_dir, results_csv, title_suffix, main_tab)

if nargin < 3, title_suffix = ''; end
if nargin < 4, main_tab = []; end

% ---- Colors (Wong / Okabe-Ito palette) ------------------------------------
w_blue      = [  0 114 178] / 255;
w_orange    = [230 159   0] / 255;
w_skyblue   = [ 86 180 233] / 255;
w_green     = [  0 158 115] / 255;  %#ok<NASGU>
w_vermilion = [213  94   0] / 255;
w_purple    = [204 121 167] / 255;
w_gray      = [0.65 0.65 0.65];     %#ok<NASGU>
w_ltgray    = [0.85 0.85 0.85];

% ---- Algorithm display order (matches plot_irlsq_results.m) ---------------
alg_csv_order  = {'CQRRT_linop', 'CholQR', 'CholQR2', 'sCholQR3_basic', 'sCholQR3'};
alg_disp_names = {'CQRRT\_linop', 'CholQR', 'CholQR2', 'sCholQR3', 'sCholQR3\_se'};

% =========================================================================
%  Load CSV (skip "#" header lines)
% =========================================================================
results_path = fullfile(data_dir, results_csv);
n_skip = count_comment_lines(results_path);
T = readtable(results_path, 'NumHeaderLines', n_skip);

algorithms    = T.algorithm;
qr_status     = T.qr_status;
qr_time       = T.qr_time_us;
peak_rss_kb   = T.peak_rss_kb;
analytical_kb = T.analytical_kb;
rspec_total   = T.rspec_total_us;
m_val         = T.m(1);
n_val         = T.n(1);
omega_val     = T.omega(1);
power_j_val   = T.power_j(1);

% Top-k eigenvalues and residuals: detect how many columns exist.
top_k = 0;
while ismember(sprintf('eig_%d', top_k), T.Properties.VariableNames) && top_k < 32
    top_k = top_k + 1;
end
eigvals = zeros(height(T), top_k);
resids  = zeros(height(T), top_k);
for k = 0:top_k-1
    eigvals(:, k+1) = T.(sprintf('eig_%d',   k));
    resids (:, k+1) = T.(sprintf('resid_%d', k));
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
for k = 1:numel(unique_algs_csv)
    if ~any(strcmp(unique_algs, unique_algs_csv{k}))
        unique_algs{end+1} = unique_algs_csv{k}; %#ok<AGROW>
        disp_labels{end+1} = strrep(unique_algs_csv{k}, '_', '\_'); %#ok<AGROW>
    end
end
n_algs = numel(unique_algs);

sel_idx    = zeros(n_algs, 1);
sel_failed = false(n_algs, 1);
for a = 1:n_algs
    mask    = strcmp(algorithms, unique_algs{a});
    indices = find(mask);
    succ    = qr_status(mask) == 0;
    if any(succ)
        succ_idx = indices(succ);
        [~, best] = min(rspec_total(succ_idx));
        sel_idx(a) = succ_idx(best);
    else
        sel_idx(a) = indices(1);
        sel_failed(a) = true;
    end
end

% =========================================================================
%  Figure
% =========================================================================
if isempty(main_tab)
    figure('Position', [100, 100, 1200, 900]);
    parent_main = gcf;
else
    parent_main = main_tab;
end
tl_main = tiledlayout(parent_main, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title_label = sprintf('RSPEC (Algorithm 4) — %d \\times %d, \\omega=%g, j=%d', ...
                       m_val, n_val, omega_val, power_j_val);
if ~isempty(title_suffix)
    title_label = sprintf('%s — %s', title_label, title_suffix);
end
title(tl_main, title_label, 'FontWeight', 'bold');

x_pos = 1:n_algs;

% ---- (1) Stacked wall-time: PCholQR + (RR + syevd + residuals) ----
nexttile(tl_main);
qr_ms    = arrayfun(@(i) qr_time(i)    / 1000, sel_idx);
total_ms = arrayfun(@(i) rspec_total(i) / 1000, sel_idx);
post_ms  = max(0, total_ms - qr_ms);
qr_ms  (sel_failed) = 0;
post_ms(sel_failed) = 0;
b = bar(x_pos, [qr_ms, post_ms], 'stacked');
b(1).FaceColor = w_blue;    b(1).DisplayName = 'PCholQR (QR)';
b(2).FaceColor = w_orange;  b(2).DisplayName = 'RR + syevd + Ritz residuals';
ylabel('Time (ms)'); title('Wall-time per algorithm');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;
for a = 1:n_algs
    if sel_failed(a)
        text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

% ---- (2) Memory: peak vs analytical ----
nexttile(tl_main);
mem_peak = arrayfun(@(i) peak_rss_kb(i),   sel_idx) / 1024;
mem_pred = arrayfun(@(i) analytical_kb(i), sel_idx) / 1024;
mem_peak(sel_failed) = 0;  mem_pred(sel_failed) = 0;
b = bar(x_pos, [mem_peak, mem_pred], 'grouped');
b(1).FaceColor = w_purple;   b(1).DisplayName = 'Peak RSS';
b(2).FaceColor = w_ltgray;   b(2).DisplayName = 'Analytical';
ylabel('Memory (MB)'); title('Peak vs predicted working memory');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
legend('Location', 'northwest'); grid on; box on;

% ---- (3) Top-k Ritz residuals (log scale, grouped by eigenvalue rank) ----
nexttile(tl_main);
res_per_alg = zeros(n_algs, top_k);
for a = 1:n_algs
    if sel_failed(a), continue; end
    res_per_alg(a, :) = max(resids(sel_idx(a), :), 1e-16);
end
b = bar(1:top_k, res_per_alg', 'grouped');
for a = 1:n_algs
    b(a).DisplayName = disp_labels{a};
end
set(gca, 'YScale', 'log');
ylim([1e-12, 1e3]);
ylabel('Ritz residual  ||K v - \lambda M v|| / (||K v|| + |\lambda| ||M v||)');
xlabel('Ritz pair rank (by |\lambda|)');
title('Top-k Ritz residuals');
legend('Location', 'northeastoutside'); grid on; box on;

% ---- (4) Top-k eigenvalues (signed, grouped by rank) ----
nexttile(tl_main);
eig_per_alg = zeros(n_algs, top_k);
for a = 1:n_algs
    if sel_failed(a), continue; end
    eig_per_alg(a, :) = eigvals(sel_idx(a), :);
end
b = bar(1:top_k, eig_per_alg', 'grouped');
for a = 1:n_algs
    b(a).DisplayName = disp_labels{a};
end
ylabel('\lambda  (signed)');
xlabel('Ritz pair rank (by |\lambda|)');
title('Top-k Ritz eigenvalues');
legend('Location', 'northeastoutside'); grid on; box on;

% ---- (5) Eigenvalue agreement: relative spread across methods ----
nexttile(tl_main);
abs_eig = abs(eig_per_alg);
abs_eig(sel_failed, :) = NaN;
spread = (max(abs_eig, [], 1, 'omitnan') - min(abs_eig, [], 1, 'omitnan')) ...
       ./ max(1e-30, max(abs_eig, [], 1, 'omitnan'));
bar(1:top_k, spread, 'FaceColor', w_skyblue);
set(gca, 'YScale', 'log');
ylim([1e-16, 1e0]);
ylabel('(max - min) / max  across surviving methods');
xlabel('Ritz pair rank (by |\lambda|)');
title('Eigenvalue agreement across methods');
grid on; box on;

% ---- (6) Total wall time per algorithm (single bar; complements panel 1) ----
nexttile(tl_main);
total_ms_clean = arrayfun(@(i) rspec_total(i) / 1000, sel_idx);
total_ms_clean(sel_failed) = 0;
bar(x_pos, total_ms_clean / 1000, 'FaceColor', w_blue);   % seconds
ylabel('Total rspec wall time (s)');
title('Total wall time');
xticks(x_pos); xticklabels(disp_labels); xtickangle(35);
grid on; box on;
for a = 1:n_algs
    if sel_failed(a)
        text(x_pos(a), 1, 'FAIL', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', w_vermilion);
    end
end

end

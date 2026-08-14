%% plot_noise_sweep.m
% Log-log plot of LS solution-error vs noise level for the sparse IR-LSQ
% benchmark. One marker per (algorithm, noise level); algorithms overlay
% because they all converge to the same LS optimum.
%
% datasets: N x 3 cell, columns = {timestamp, noise_level, label}.
%   Rows with noise_level == 0 are plotted as a horizontal dashed line
%   (machine-epsilon floor) rather than as points on the log axis.
%
% Usage:
%   plot_noise_sweep(data_dir, datasets, title_suffix)
%   plot_noise_sweep(data_dir, datasets, title_suffix, parent)
%   parent may be a figure / uitab handle; defaults to a new figure.

function fig = plot_noise_sweep(data_dir, datasets, title_suffix, parent)

if nargin < 3, title_suffix = ''; end
if nargin < 4, parent = []; end

% Wong colorblind-friendly palette
w_blue      = [  0 114 178] / 255;
w_orange    = [230 159   0] / 255;
w_green     = [  0 158 115] / 255;
w_vermilion = [213  94   0] / 255;
w_purple    = [204 121 167] / 255;
w_sky       = [ 86 180 233] / 255;   % Wong sky-blue, for the newly-added CholQR2 series
w_gray      = [0.5 0.5 0.5];

% Algorithm rendering order + colors + markers
% CholQR2 was silently missing from this list, so it was dropped from every noise-sweep
% figure even when present in the CSV (all other plotters include it). Added 2026-07-27.
% Display names for the two sCholQR3 variants are spelled out rather than relying on the
% CSV suffix: CSV 'sCholQR3_basic' is the NON-blocked one and CSV 'sCholQR3' is the
% blocked/storage-efficient one, which reads backwards if left implicit.
alg_csv   = {'CQRRT_linop',  'CQRRT_linop_bqrrp',   'CholQR',  'CholQR2', 'sCholQR3_basic',   'sCholQR3'};
alg_disp  = {'CQRRT\_linop', 'CQRRT\_linop\_bqrrp', 'CholQR',  'CholQR2', 'sCholQR3 (basic)', 'sCholQR3 (blocked)'};
alg_color = {w_blue,         w_orange,              w_green,   w_sky,     w_vermilion,        w_purple};
alg_mark  = {'o',            's',                   '^',       'x',       'd',                'v'};

n_algs = numel(alg_csv);
n_ds   = size(datasets, 1);

% Allocate: (algorithm, dataset) -> sol_err (NaN if missing/failed)
sol_err_mat = nan(n_algs, n_ds);
resid_mat   = nan(n_algs, n_ds);
noise_vec   = zeros(n_ds, 1);

for j = 1:n_ds
    ts             = datasets{j, 1};
    noise_vec(j)   = datasets{j, 2};
    res_csv        = fullfile(data_dir, sprintf('%s_irlsq_results.csv', ts));
    if ~isfile(res_csv)
        warning('plot_noise_sweep: missing %s', res_csv);
        continue;
    end
    n_skip = count_comment_lines(res_csv);
    opts   = detectImportOptions(res_csv, 'NumHeaderLines', n_skip);
    opts   = setvartype(opts, 'algorithm', 'char');
    T      = readtable(res_csv, opts);

    for a = 1:n_algs
        mask  = strcmp(T.algorithm, alg_csv{a}) & T.qr_status == 0;
        if ~any(mask), continue; end
        % All runs of a given algorithm hit the same LS optimum to >= 6
        % significant digits, so the median across runs is uncontroversial.
        sol_err_mat(a, j) = median(T.ls_solution_error(mask));
        resid_mat(a, j)   = median(T.ls_residual_norm(mask));
    end
end

% Split noise=0 row(s) from the log-plottable rows
zero_idx    = noise_vec == 0;
nonzero_idx = ~zero_idx;
noise_pos   = noise_vec(nonzero_idx);
sol_err_pos = sol_err_mat(:, nonzero_idx);

% Empirical amplification factor (median over non-zero rows, median over algos)
amp_factor = median(sol_err_pos(:) ./ ...
                    repmat(noise_pos(:)', n_algs, 1), 'all', 'omitnan');

if isempty(parent)
    fig = figure('Position', [100, 100, 850, 600]);
    parent_fig = fig;
else
    parent_fig = parent;
    fig = ancestor(parent_fig, 'figure');
end
ax = axes(parent_fig); hold(ax, 'on');

% Reference line: y = amp_factor * x, drawn across the plotted noise range
x_lo = min(noise_pos) / 2;
x_hi = max(noise_pos) * 2;
plot(ax, [x_lo, x_hi], amp_factor * [x_lo, x_hi], '--', ...
     'Color', w_gray, 'LineWidth', 1.3, ...
     'DisplayName', sprintf('y = %.2f \\cdot x  (LS noise amplification)', amp_factor));

% One series per algorithm (lines + markers overlay because all algorithms
% converge to the same LS optimum)
for a = 1:n_algs
    y = sol_err_pos(a, :);
    if all(isnan(y)), continue; end
    plot(ax, noise_pos, y, '-', 'Color', alg_color{a}, 'LineWidth', 1.2, ...
         'Marker', alg_mark{a}, 'MarkerSize', 9, 'MarkerFaceColor', alg_color{a}, ...
         'DisplayName', alg_disp{a});
end

% Machine-eps floor from noise=0 rows (if any)
if any(zero_idx)
    floor_y = max(sol_err_mat(:, zero_idx), [], 'all', 'omitnan');
    plot(ax, [x_lo, x_hi], [floor_y, floor_y], ':', ...
         'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
         'DisplayName', sprintf('noise = 0 floor (\\approx %.1e)', floor_y));
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, [x_lo, x_hi]);
xlabel(ax, 'noise level  (||noise|| / ||A x_{true}||)');
ylabel(ax, '||x - x_{true}|| / ||x_{true}||');
title(ax, sprintf('LS solution error vs noise — %s', title_suffix), 'Interpreter', 'tex');
grid(ax, 'on'); box(ax, 'on');
legend(ax, 'Location', 'southeast');
end

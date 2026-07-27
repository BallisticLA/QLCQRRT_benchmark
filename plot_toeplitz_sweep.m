%% plot_toeplitz_sweep.m — cross-size summary of the Toeplitz LS 2:1 sweep.
%
% Reads the newest *_toeplitz_ls_results.csv in every size subfolder (n%05d,
% one ISAAC experiment each, m = 2n) and plots each metric as a line per method
% versus n. Panel 1 (total wall-time) is the preconditioned-vs-unpreconditioned
% crossover picture the sweep was run for.
%
% Empty CSVs (failed/capped jobs) and qr_status~=0 rows are skipped as NaN
% points, so a missing size leaves a gap instead of breaking the figure.
%
% Usage:
%   plot_toeplitz_sweep(toep_dir, size_dirs, main_tab)

function plot_toeplitz_sweep(toep_dir, size_dirs, main_tab)

if nargin < 3, main_tab = []; end

% Wong colorblind-friendly palette (matches plot_toeplitz_results).
w_blue = [0 114 178]/255;  w_orange = [230 159 0]/255;  w_skyblue = [86 180 233]/255;
w_green = [0 158 115]/255; w_vermilion = [213 94 0]/255; w_purple = [204 121 167]/255;
w_dkgray = [0.3 0.3 0.3];

% Same CSV-name -> display-name mapping as the per-size plotter.
alg_csv_order  = {'CQRRT_linop','CholQR','CholQR2','sCholQR3_basic','sCholQR3','Blendenpik','unpreconditioned'};
alg_disp_names = {'CQRRT\_linop','CholQR','CholQR2','sCholQR3 (basic)','sCholQR3 (blocked)','Blendenpik','unprec'};
alg_colors     = {w_blue, w_orange, w_skyblue, w_green, w_vermilion, w_purple, w_dkgray};
alg_markers    = {'o','s','d','^','v','>','x'};
na = numel(alg_csv_order);

% --- Gather: one row per size, one column per method, NaN where absent. ---
ns = numel(size_dirs);
n_vals   = nan(ns, 1);
build_s  = nan(ns, na);  solve_s = nan(ns, na);
iters    = nan(ns, na);  rec_err = nan(ns, na);  orth = nan(ns, na);
for s = 1:ns
    size_dir = fullfile(toep_dir, size_dirs{s});
    f = dir(fullfile(size_dir, '*_toeplitz_ls_results.csv'));
    f = f([f.bytes] > 0);                      % empty CSV = job died before writing
    if isempty(f)
        fprintf('(sweep: skipped %s — no non-empty CSV)\n', size_dirs{s});
        continue;
    end
    [~, ord] = sort({f.name});                 % newest last (timestamped names)
    results_path = fullfile(size_dir, f(ord(end)).name);
    n_skip = count_comment_lines(results_path);
    opts = detectImportOptions(results_path, 'NumHeaderLines', n_skip);
    opts = setvartype(opts, 'algorithm', 'char');
    T = readtable(results_path, opts);
    n_vals(s) = T.n(1);
    for k = 1:na
        idx = find(strcmp(T.algorithm, alg_csv_order{k}), 1);
        if isempty(idx) || T.qr_status(idx) ~= 0, continue; end
        build_s(s, k) = T.qr_time_us(idx) / 1e6;
        solve_s(s, k) = T.solve_time_us(idx) / 1e6;
        iters(s, k)   = T.iterations(idx);
        rec_err(s, k) = T.recovery_error(idx);
        orth(s, k)    = T.orth_error(idx);
    end
end
rec_err(rec_err < 0) = NaN;
orth(orth < 0) = NaN;                          % -1 => no preconditioner
have = ~isnan(n_vals);
n_vals = n_vals(have);
build_s = build_s(have, :); solve_s = solve_s(have, :);
iters = iters(have, :); rec_err = rec_err(have, :); orth = orth(have, :);
total_s = build_s + solve_s;
total_s(:, na) = solve_s(:, na);               % unprec: build = 0 by definition

% --- Figure / tab ---
if isempty(main_tab), figure('Position',[100 100 1250 900]); parent = gcf; else parent = main_tab; end
tl = tiledlayout(parent, 2, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('Toeplitz LS sweep — m = 2n, n = %d ... %d', min(n_vals), max(n_vals)), ...
      'FontWeight','bold');

lineargs = @(k) {'Color', alg_colors{k}, 'Marker', alg_markers{k}, ...
                 'MarkerFaceColor', alg_colors{k}, 'LineWidth', 1.3, ...
                 'DisplayName', alg_disp_names{k}};

    function sweep_panel(vals, ylab, ttl, yscale)
        nexttile(tl); hold on;
        for kk = 1:na
            args = lineargs(kk);
            if strcmp(alg_csv_order{kk}, 'unpreconditioned'), args = [args, {'LineStyle','--'}]; end %#ok<AGROW>
            plot(n_vals, vals(:, kk), args{:});
        end
        set(gca, 'XScale','log', 'YScale', yscale);
        xticks(n_vals); xtickangle(35);
        xlabel('n'); ylabel(ylab); title(ttl); grid on; box on;
    end

% (1) Total wall-time — the crossover panel.
sweep_panel(total_s, 'time (s)', 'Total wall-time (build + solve)', 'log');
legend('Location','northwest', 'NumColumns', 2);

% (2) Preconditioner build (QR) time. unprec has none -> absent.
sweep_panel(build_s, 'time (s)', 'Preconditioner build (QR) time', 'log');

% (3) LSQR solve time.
sweep_panel(solve_s, 'time (s)', 'LSQR solve time', 'log');

% (4) LSQR iterations to convergence.
sweep_panel(iters, 'iterations', 'LSQR iterations', 'log');

% (5) Recovery error ||x - x_true|| / ||x_true||.
sweep_panel(rec_err, '||x-x_{true}||/||x_{true}||', 'Recovery error', 'log');

% (6) Preconditioner orthogonality loss (N/A for unprec).
sweep_panel(orth, '||R^{-T}A^TA R^{-1}-I||_F/\surd n', 'Orthogonality loss', 'log');

end

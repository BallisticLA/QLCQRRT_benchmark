function plot_irlsq_reg_sweep(data_dir, method_name)
% plot_irlsq_reg_sweep — visualize the irlsq_reg kappa x precision sweep.
%
% Reads every *_irlsq_reg_results.csv under data_dir (recursively), then for each
% matrix size produces a 2x2 figure of the four key metrics vs kappa_target, with
% one line per precision combo (precond_prec x solve_prec):
%
%   (1) Forward error   ||x - x_true|| / ||x_true||   (~ u*kappa; precision-sensitive)
%   (2) Backward error  ||Ax - b|| / (||A|| ||x|| + ||b||)  (~ u(solve); kappa-flat)
%   (3) Inner CG iters to convergence
%   (4) Q-factor orthogonality loss ||Q^T Q - I||_F / sqrt(n)
%
% The headline result: single-precond / double-solve (sd) tracks the all-double
% baseline (dd) in forward error even at large kappa, while single/single (ss)
% degrades as ~u_single*kappa.
%
% Usage:
%   plot_irlsq_reg_sweep(data_dir)               % method = 'CholQR'
%   plot_irlsq_reg_sweep(data_dir, 'CQRRT_linop')

if nargin < 2 || isempty(method_name), method_name = 'CholQR'; end

% ---- Colors (Wong / Okabe-Ito) keyed by precision combo ----
combo_color = containers.Map();
combo_color('double/double') = [0 0 0];                 % baseline: black
combo_color('single/double') = [0 114 178] / 255;       % scenario 1: blue
combo_color('single/single') = [213 94 0] / 255;        % scenario 2: vermilion
combo_color('double/single') = [0 158 115] / 255;       % completeness: green

% ---- Collect all reg CSVs ----
files = dir(fullfile(data_dir, '**', '*_irlsq_reg_results.csv'));
if isempty(files)
    error('plot_irlsq_reg_sweep: no *_irlsq_reg_results.csv under %s', data_dir);
end
T = table();
klab_all = strings(0, 1);     % kappa-variant label per row, parsed from the cell folder
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    n_skip = count_comment_lines(p);
    t = readtable(p, 'NumHeaderLines', n_skip, 'TextType', 'string');
    [~, foldername] = fileparts(files(k).folder);   % e.g. small_k1e7_dd
    kl = regexp(foldername, 'k1e\d+', 'match', 'once');
    if isempty(kl), kl = "k?"; end
    klab_all = [klab_all; repmat(string(kl), height(t), 1)]; %#ok<AGROW>
    T = [T; t]; %#ok<AGROW>
end

% String-ify precision columns and form the combo label "precond/solve".
precond = string(T.precond_prec);
solve   = string(T.solve_prec);
combo   = precond + "/" + solve;

% Keep only the requested method and successful QR runs.
keep = (string(T.algorithm) == string(method_name)) & (T.qr_status == 0);
T = T(keep, :); precond = precond(keep); solve = solve(keep); combo = combo(keep); %#ok<NASGU>
klab = klab_all(keep);

if isempty(T)
    error('plot_irlsq_reg_sweep: no successful rows for method %s', method_name);
end

% kappa_target is now 1 for every cell (kappa lives in the data, not runtime
% scaling), so the x-axis must come from kappa_measured. Use one canonical kappa
% per variant -- the max over that variant's rows, i.e. the double/double value
% (single-precision sCholQR3 distorts its own kappa_measured downward).
ulab = unique(klab);
kap_of_lab = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(ulab)
    kap_of_lab(char(ulab(i))) = max(T.kappa_measured(klab == ulab(i)));
end
kappa_plot = zeros(height(T), 1);
for r = 1:height(T)
    kappa_plot(r) = kap_of_lab(char(klab(r)));
end

sizes = unique(T.m, 'stable');
combos_order = ["double/double", "single/double", "single/single", "double/single"];

for s = 1:numel(sizes)
    m_s = sizes(s);
    rows_s = (T.m == m_s);
    n_s = T.n(find(rows_s, 1));

    fig = figure('Name', sprintf('irlsq_reg sweep — %s — %d x %d', method_name, m_s, n_s), ...
                 'Position', [80 + 30*s, 80, 1200, 850]);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Regularized IR-LSQ  (%s)   %d \\times %d', ...
          strrep(method_name,'_','\_'), m_s, n_s), 'FontWeight', 'bold');

    % Panel definitions: {field, title, ylabel, yscale}
    panels = {
        'ls_solution_error', 'Forward error vs \kappa', '||x - x_{true}|| / ||x_{true}||', 'log';
        'ls_residual_norm',  'Backward error vs \kappa', '||Ax-b|| / (||A|| ||x|| + ||b||)', 'log';
        'ir_inner_iters_total', 'Inner CG iterations vs \kappa', 'inner CG iters (total)', 'linear';
        'orth_error',        'Q-factor orthogonality loss vs \kappa', '||Q^T Q - I||_F / \surd n', 'log';
    };

    for pp = 1:size(panels, 1)
        ax = nexttile(tl); hold(ax, 'on');
        fld = panels{pp, 1};
        for c = 1:numel(combos_order)
            cb = combos_order(c);
            sel = rows_s & (combo == cb);
            if ~any(sel), continue; end
            kap = kappa_plot(sel);
            val = T.(fld)(sel);
            [kap, idx] = sort(kap); val = val(idx);
            % Average duplicate kappa (multiple runs) if present.
            [ukap, ~, g] = unique(kap);
            uval = accumarray(g, val, [], @mean);
            col = combo_color(char(cb));
            plot(ax, ukap, max(uval, realmin), '-o', 'Color', col, ...
                 'MarkerFaceColor', col, 'LineWidth', 1.5, 'DisplayName', char(cb));
        end
        set(ax, 'XScale', 'log');
        if strcmp(panels{pp, 4}, 'log'), set(ax, 'YScale', 'log'); end
        xlabel(ax, '\kappa (measured)'); ylabel(ax, panels{pp, 3});
        title(ax, panels{pp, 2}); grid(ax, 'on'); box(ax, 'on');
        if pp == 1, legend(ax, 'Location', 'northwest'); end
    end
end
end

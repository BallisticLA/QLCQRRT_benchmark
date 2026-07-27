% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));

tab_groups = {};   % {tabgroup handle, export prefix}

%% ========================================================================
%  Toeplitz least-squares benchmark (Oleg's 2nd experiment) — ACTIVE
%
%  Prolate-Toeplitz regularized LS  min ||T x - b||^2 + lambda ||x||^2, augmented
%  A = [T; sqrt(lambda) I], solved by LSQR with Q-less QR right preconditioners.
%  6-panel plot_toeplitz_results layout: wall-time (build+solve), data & recovery
%  error, peak-vs-analytical memory, LSQR iterations, preconditioner orthogonality
%  loss, Cholesky shift retries. Methods: CQRRT / CholQR / CholQR2 / sCholQR3 /
%  sCholQR3_basic / Blendenpik / unpreconditioned. Data from the C++ benchmark
%  (toeplitz_ls_benchmark), one CSV per run.
%
%  First tab = plot_toeplitz_sweep cross-size summary (metric-vs-n lines per
%  method; the build-vs-solve crossover picture), then one tab per size.
% =========================================================================

toep_dir = fullfile(script_dir, 'benchmark-output', 'toeplitz_ls');
if exist(toep_dir, 'dir')
    % Discover every size subfolder (one ISAAC experiment each) that holds a CSV,
    % sorted by name (subfolders are zero-padded n%05d so this is size order).
    d = dir(toep_dir);
    size_dirs = sort({d([d.isdir] & ~ismember({d.name}, {'.','..'})).name});
    fig = figure('Name', 'Toeplitz LS — Q-less QR right preconditioners', ...
                 'Position', [80, 80, 1250, 900]);
    tg = uitabgroup(fig);
    st = uitab(tg, 'Title', 'sweep');            % cross-size summary, first tab
    plot_toeplitz_sweep(toep_dir, size_dirs, st);
    any_size = false;
    for s = 1:numel(size_dirs)
        size_dir = fullfile(toep_dir, size_dirs{s});
        f = dir(fullfile(size_dir, '*_toeplitz_ls_results.csv'));
        f = f([f.bytes] > 0);                    % empty CSV = job died before writing
        if isempty(f)
            fprintf('(skipped Toeplitz LS %s — no non-empty CSV)\n', size_dirs{s});
            continue;
        end
        [~, ord] = sort({f.name});          % newest last (timestamped names)
        res_csv = f(ord(end)).name;
        % Label the tab with the FULL m x n, not the folder name. The folder is
        % named after n only (n%05d), so a tab reading "n16000" sits above a plot
        % titled "32000 x 16000" -- an easy misread of the tab as the matrix shape
        % (hit 2026-07-27). Read the shape out of the CSV so the two always agree.
        shape = toeplitz_shape_label(fullfile(size_dir, res_csv), size_dirs{s});
        mt = uitab(tg, 'Title', shape);
        plot_toeplitz_results(size_dir, res_csv, ...
            sprintf('prolate Toeplitz, LSQR — %s', shape), mt);
        any_size = true;
    end
    if any_size
        tab_groups(end+1, :) = {tg, 'toeplitz_ls'}; %#ok<SAGROW>
    else
        close(fig);
    end
else
    fprintf('(skipped Toeplitz LS — no data dir %s)\n', toep_dir);
end

%% ========================================================================
%  FEM_Problem_2 App-1 IR-LSQ (regularized) + Blendenpik, single hard
%  kappa^colnorm=1e6, two sizes, three precision combos (ISAAC 2026-07-11) — ACTIVE.
%
%  min ||A x - b||^2 solved by LSQR with Q-less QR right preconditioners, plus the
%  Blendenpik competitor (SASO sketch + Householder QR + LSQR). Per combo: one
%  6-panel results tab per size. Runtime-breakdown tabs are intentionally skipped
%  (pass an empty breakdown CSV name to suppress them).
%  Data: benchmark-output/irlsq_reg_1e6/<size>_kc1e6_<combo>/.
% =========================================================================
if true
    clean_dir = fullfile(script_dir, 'benchmark-output', 'irlsq_reg_1e6');
    sizes  = { 'small', 'small (75466x8256)';
               'large', 'large (301864x33024)' };
    kaps   = { 'kc1e6', '\kappa^{colnorm}\approx1e6' };
    combos = { 'dd', 'double precond / double solve';
               'sd', 'single precond / double solve';
               'ss', 'single precond / single solve' };
    if exist(clean_dir, 'dir')
        for c = 1:size(combos, 1)
            fig = figure('Name', sprintf('FEM2 App 1 — IR-LSQ (reg) — %s', combos{c, 2}), ...
                         'Position', [60 + 25*c, 60, 1150, 900]);
            tg = uitabgroup(fig);
            for s = 1:size(sizes, 1)
                for k = 1:size(kaps, 1)
                    cell_name = sprintf('%s_%s_%s', sizes{s, 1}, kaps{k, 1}, combos{c, 1});
                    cell_dir  = fullfile(clean_dir, cell_name);
                    f = dir(fullfile(cell_dir, '*_irlsq_reg_results.csv'));
                    if isempty(f), fprintf('(missing cell: %s)\n', cell_name); continue; end
                    tab_title = sprintf('%s  k~%s', sizes{s, 1}, kaps{k, 2});
                    suffix    = sprintf('\\kappa \\approx %s   [%s]', kaps{k, 2}, combos{c, 1});
                    mt = uitab(tg, 'Title', tab_title);
                    plot_irlsq_results(cell_dir, f(1).name, '', suffix, mt);
                end
            end
            tab_groups(end+1, :) = {tg, sprintf('fem2_irlsq_reg_%s', combos{c, 1})}; %#ok<SAGROW>
        end
    end
end

%% ========================================================================
%  Export — one vector PDF per tab.
% =========================================================================
export_dir = fullfile(script_dir, 'figures-export');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
% Wipe ALL prior PDFs so figures-export/ holds only the current run's tabs.
old = dir(fullfile(export_dir, '*.pdf'));
for i = 1:numel(old), delete(fullfile(export_dir, old(i).name)); end
expfig = @(h, name) exportgraphics(h, fullfile(export_dir, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

for g = 1:size(tab_groups, 1)
    tg     = tab_groups{g, 1};
    prefix = tab_groups{g, 2};
    for k = 1:numel(tg.Children)
        tab  = tg.Children(k);
        slug = lower(regexprep(tab.Title, '[^a-zA-Z0-9]+', '_'));
        slug = regexprep(slug, '^_+|_+$', '');
        fname = [prefix, '_', slug];
        expfig(tab, fname);
        fprintf('  Exported: %s.pdf\n', fname);
    end
end
fprintf('All figures exported to %s\n', export_dir);

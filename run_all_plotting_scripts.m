% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));                % apply_butterfly, compute_kappa_variants

%% ========================================================================
%  FEM_Problem_2 — Application 1 (IR-LSQ, regularized augmented operator)
%                  SINGLE hard kappa^colnorm = 1e6 campaign (ISAAC 2026-07-11)
%
%  6-panel plot_irlsq_results layout: wall-time, normwise backward error, peak-vs-
%  analytical memory, inner-CG (or LSQR) iters, Q-orth loss, + Cholesky shift retries.
%  Data = geometric-scale + sparse-butterfly matrices (gen_fem2_hard_variants),
%  kappa^colnorm=1e6 (the regime where single-precond WORKS for the robust methods;
%  kappa>=1e7 was uniformly hopeless). Run with kappa_target=1, mu_factor=100, mask=63
%  (6 methods incl. Blendenpik = sketch + Householder-QR + LSQR, an independent solver).
%
%  One tabbed figure per precision combo, 2 tabs each (small + large):
%    dd = double/double  -- all Q-less methods machine-accurate; Blendenpik converges
%         to its LSQR tol (looser, more iters).
%    sd = single precond / double solve -- CQRRT is the unique winner (bwd 1e-16 from a
%         single-precision sketch preconditioner); the CholeskyQR-of-Gram methods fail.
%    ss = single/single -- all fail (pure single at kappa~2.5e6), expected.
%  Sizes: small 75466x8256, large 301864x33024. Cell dirs carry a breakdown CSV.
% =========================================================================

clean_dir = fullfile(script_dir, 'benchmark-output', 'irlsq_reg_1e6');
sizes  = { 'small', 'small (75466x8256)';
           'large', 'large (301864x33024)' };
kaps   = { 'kc1e6', '\kappa^{colnorm}\approx1e6' };
combos = { 'dd', 'double precond / double solve';
           'sd', 'single precond / double solve';
           'ss', 'single precond / single solve' };

tab_groups = {};   % {tabgroup handle, export prefix}
if exist(clean_dir, 'dir')
    for c = 1:size(combos, 1)
        fig = figure('Name', sprintf('FEM2 App 1 — IR-LSQ (reg) — %s', combos{c, 2}), ...
                     'Position', [60 + 25*c, 60, 1150, 900]);
        tg = uitabgroup(fig);
        fig_bd = figure('Name', sprintf('FEM2 App 1 — runtime breakdown — %s', combos{c, 2}), ...
                        'Position', [80 + 25*c, 80, 1100, 600]);
        tg_bd = uitabgroup(fig_bd);
        for s = 1:size(sizes, 1)
            for k = 1:size(kaps, 1)
                cell_name = sprintf('%s_%s_%s', sizes{s, 1}, kaps{k, 1}, combos{c, 1});
                cell_dir  = fullfile(clean_dir, cell_name);
                f = dir(fullfile(cell_dir, '*_irlsq_reg_results.csv'));
                if isempty(f)
                    fprintf('(missing cell: %s)\n', cell_name);
                    continue;
                end
                fb = dir(fullfile(cell_dir, '*_irlsq_reg_breakdown.csv'));
                bd_csv = '';  if ~isempty(fb), bd_csv = fb(1).name; end
                tab_title = sprintf('%s  k~%s', sizes{s, 1}, kaps{k, 2});
                suffix    = sprintf('\\kappa \\approx %s   [%s]', kaps{k, 2}, combos{c, 1});
                mt = uitab(tg,    'Title', tab_title);
                bt = uitab(tg_bd, 'Title', tab_title);
                plot_irlsq_results(cell_dir, f(1).name, bd_csv, suffix, mt, bt);
            end
        end
        tab_groups(end+1, :) = {tg,    sprintf('fem2_irlsq_reg_%s',    combos{c, 1})}; %#ok<SAGROW>
        tab_groups(end+1, :) = {tg_bd, sprintf('fem2_irlsq_reg_%s_bd', combos{c, 1})}; %#ok<SAGROW>
    end
else
    fprintf('(skipped App 1 IR-LSQ — no data in %s)\n', clean_dir);
end

%% ========================================================================
%  OLD — superseded. Kept (guarded) for reference: the ISAAC 2026-06-12
%  single-kappa App-1 (IR-LSQ) + App-2 (RSPEC) tabs on the pre-contrast FEM2
%  (the moot kappa~2e13 _corrected operator). Flip the guard to `if true` to revive.
% =========================================================================
if false
    data_dir_isaac_jun12_fem2 = fullfile(script_dir, 'benchmark-output', 'isaac-jun12-fem2');

    datasets_fem2_jun12_irlsq = {
        '20260612_034608',  'FEM_Problem_2 small  (75824x8304)';
        '20260612_035744',  'FEM_Problem_2 medium (151648x16608)';
        '20260612_043556',  'FEM_Problem_2 large  (303296x33216)';
    };
    fig_fem2    = figure('Name', 'FEM_Problem_2 App 1 — Sparse IR-LSQ (ISAAC 2026-06-12)',  'Position', [ 50, 50, 1200, 900]);
    tg_fem2     = uitabgroup(fig_fem2);
    fig_fem2_bd = figure('Name', 'FEM_Problem_2 App 1 — Runtime Breakdowns (ISAAC 2026-06-12)', 'Position', [100, 50, 1100, 600]);
    tg_fem2_bd  = uitabgroup(fig_fem2_bd);
    if exist(data_dir_isaac_jun12_fem2, 'dir')
        for i = 1:size(datasets_fem2_jun12_irlsq, 1)
            ts    = datasets_fem2_jun12_irlsq{i, 1};
            label = datasets_fem2_jun12_irlsq{i, 2};
            mt = uitab(tg_fem2,    'Title', label);
            bt = uitab(tg_fem2_bd, 'Title', label);
            res_csv = sprintf('%s_irlsq_results.csv',   ts);
            bd_csv  = sprintf('%s_irlsq_breakdown.csv', ts);
            plot_irlsq_results(data_dir_isaac_jun12_fem2, res_csv, bd_csv, label, mt, bt);
        end
    else
        fprintf('(skipped FEM_Problem_2 IR-LSQ tabs — no data in %s)\n', data_dir_isaac_jun12_fem2);
    end

    datasets_fem2_jun12_rspec = {
        '20260612_040233',  'FEM_Problem_2 small  (75824x8304), \omega=0, j=1';
    };
    fig_fem2_rspec = figure('Name', 'FEM_Problem_2 App 2 — RSPEC (ISAAC 2026-06-12)',  'Position', [150, 50, 1200, 900]);
    tg_fem2_rspec  = uitabgroup(fig_fem2_rspec);
    if exist(data_dir_isaac_jun12_fem2, 'dir')
        for i = 1:size(datasets_fem2_jun12_rspec, 1)
            ts    = datasets_fem2_jun12_rspec{i, 1};
            label = datasets_fem2_jun12_rspec{i, 2};
            res_csv = sprintf('%s_rspec_results.csv', ts);
            csv_path = fullfile(data_dir_isaac_jun12_fem2, res_csv);
            if isfile(csv_path)
                mt = uitab(tg_fem2_rspec, 'Title', label);
                plot_rspec_results(data_dir_isaac_jun12_fem2, res_csv, label, mt);
            else
                fprintf('(missing rspec CSV: %s)\n', csv_path);
            end
        end
    else
        fprintf('(skipped FEM_Problem_2 RSPEC tabs — no data in %s)\n', data_dir_isaac_jun12_fem2);
    end
end  % if false (old jun12 App1/App2)

%% ========================================================================
%  Export — one vector PDF per tab.
% =========================================================================
export_dir = fullfile(script_dir, 'figures-export');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
% Wipe this script's own prior outputs so figures-export/ holds ONLY the current
% campaign's figures (stale kc1e7/1e9/1e11 etc. from earlier runs are removed).
% Scoped to the fem2_irlsq_reg_* prefix so other projects' figures are untouched.
old = dir(fullfile(export_dir, 'fem2_irlsq_reg_*.pdf'));
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

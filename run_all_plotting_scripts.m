% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));                % resize_fem_pair, resize_sparse_matrix, etc.

%% ========================================================================
%  FEM_Problem_2 — Application 1 (IR-LSQ, regularized augmented operator)
%                  ISAAC 2026-06-20 kappa x precision campaign
%
%  Same 5-panel plot_irlsq_results layout as always (wall-time, normwise
%  backward error, peak-vs-analytical memory, inner-CG iters, Q-orth loss) —
%  just driven by the new irlsq_reg data. Each cell CSV is one config
%  (5 Q-less QR methods, one matrix, one precision combo), so it drops straight
%  into plot_irlsq_results; the extra reg columns (kappa_measured, mu,
%  precond/solve_prec) are ignored by that plotter.
%
%  One tabbed figure per precision combo, 9 tabs each (3 sizes x 3 kappa):
%    dd = double/double (baseline)
%    sd = single precond / double solve   (collaborator scenario 1)
%    ss = single precond / single solve   (collaborator scenario 2)
%  Achieved kappa (kappa_measured, dd cells): k1e7~8.8e5, k1e9~2.3e8, k1e11~6.0e10.
%  No breakdown CSV is emitted for irlsq_reg, so the phase-breakdown figure is
%  skipped automatically.
% =========================================================================

clean_dir = fullfile(script_dir, 'benchmark-output', 'irlsq_reg_3kappa_clean_v2');
sizes  = { 'small',  'small (75824x8304)';
           'medium', 'medium (151648x16608)';
           'large',  'large (303296x33216)' };
kaps   = { 'k1e7',  '8.8e5';
           'k1e9',  '2.3e8';
           'k1e11', '6.0e10' };
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

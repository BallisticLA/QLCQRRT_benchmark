% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));                % resize_fem_pair, resize_sparse_matrix, etc.

data_dir_irlsq            = fullfile(script_dir, 'benchmark-output', 'irlsq');
data_dir_isaac_jun5_fem2  = fullfile(script_dir, 'benchmark-output', 'isaac-jun5-fem2');
data_dir_isaac_apr7       = fullfile(script_dir, 'benchmark-output', 'isaac-apr7');     %#ok<NASGU>
data_dir_isaac_apr9       = fullfile(script_dir, 'benchmark-output', 'isaac-apr9');     %#ok<NASGU>
data_dir_photogrammetry   = fullfile(script_dir, 'benchmark-output', 'photogrammetry2'); %#ok<NASGU>
data_dir_nmr              = fullfile(script_dir, 'benchmark-output', 'nmr');             %#ok<NASGU>

%% ========================================================================
%  FEM_Problem_2 — ISAAC 2026-06-06 campaign (Application 1: IR-LSQ +
%                                              Application 2: RSPEC)
%
%  Binary HEAD = 98627b5 (CholQR family rework: adaptive shift on
%  cholqr_primitive + pcholqr_primitive; new CholQR2 driver;
%  sCholQR3_linops_basic refactored onto the shared primitives;
%  CQRRT_linop_bqrrp dropped from dispatch). method_mask = 31
%  (CQRRT_linop, CholQR, sCholQR3, sCholQR3_basic, CholQR2).
%  d_factor=2, sketch_nnz=4, block_size=256.
%  ISAAC SPR ai-tenn partition, 128 OpenMP threads.
%
%  IR-LSQ status: all 3 sizes complete. With the new adaptive shift,
%    sCholQR3 + sCholQR3_basic + CholQR2 all survive at medium+large
%    (previously failed). Plain CholQR still fails at medium (intentional;
%    no retries on the unshifted baseline).
%  RSPEC status: small completed; medium + large were cancelled at
%    ~2h05 by the ai-tenn QoS cap (my 12h / 24h walltime requests were
%    overridden). Only small is plotted; medium + large need a per-method
%    split or QoS bump to land.
% =========================================================================

% --- Application 1: IR-LSQ (small / medium / large) ----------------------
datasets_fem2_jun6_irlsq = {
    '20260606_054624',  'FEM_Problem_2 small  (75824x8304)';
    '20260606_055818',  'FEM_Problem_2 medium (151648x16608)';
    '20260606_065439',  'FEM_Problem_2 large  (303296x33216)';
};

fig_fem2    = figure('Name', 'FEM_Problem_2 App 1 — Sparse IR-LSQ (ISAAC 2026-06-06)',  'Position', [ 50, 50, 1200, 900]);
tg_fem2     = uitabgroup(fig_fem2);
fig_fem2_bd = figure('Name', 'FEM_Problem_2 App 1 — Runtime Breakdowns (ISAAC 2026-06-06)', 'Position', [100, 50, 1100, 600]);
tg_fem2_bd  = uitabgroup(fig_fem2_bd);

if exist(data_dir_isaac_jun5_fem2, 'dir')
    for i = 1:size(datasets_fem2_jun6_irlsq, 1)
        ts    = datasets_fem2_jun6_irlsq{i, 1};
        label = datasets_fem2_jun6_irlsq{i, 2};
        mt = uitab(tg_fem2,    'Title', label);
        bt = uitab(tg_fem2_bd, 'Title', label);
        res_csv = sprintf('%s_irlsq_results.csv',   ts);
        bd_csv  = sprintf('%s_irlsq_breakdown.csv', ts);
        plot_irlsq_results(data_dir_isaac_jun5_fem2, res_csv, bd_csv, label, mt, bt);
    end
else
    fprintf('(skipped FEM_Problem_2 IR-LSQ tabs — no data in %s)\n', data_dir_isaac_jun5_fem2);
end

% --- Application 2: RSPEC (small only; medium / large pending QoS bump) --
datasets_fem2_jun6_rspec = {
    '20260606_074538',  'FEM_Problem_2 small  (75824x8304), \omega=0, j=1';
};

fig_fem2_rspec = figure('Name', 'FEM_Problem_2 App 2 — RSPEC (ISAAC 2026-06-06)',  'Position', [150, 50, 1200, 900]);
tg_fem2_rspec  = uitabgroup(fig_fem2_rspec);

if exist(data_dir_isaac_jun5_fem2, 'dir')
    for i = 1:size(datasets_fem2_jun6_rspec, 1)
        ts    = datasets_fem2_jun6_rspec{i, 1};
        label = datasets_fem2_jun6_rspec{i, 2};
        res_csv = sprintf('%s_rspec_results.csv', ts);
        csv_path = fullfile(data_dir_isaac_jun5_fem2, res_csv);
        if isfile(csv_path)
            mt = uitab(tg_fem2_rspec, 'Title', label);
            plot_rspec_results(data_dir_isaac_jun5_fem2, res_csv, label, mt);
        else
            fprintf('(missing rspec CSV: %s)\n', csv_path);
        end
    end
else
    fprintf('(skipped FEM_Problem_2 RSPEC tabs — no data in %s)\n', data_dir_isaac_jun5_fem2);
end

%% ========================================================================
%  Export
% =========================================================================

export_dir = fullfile(script_dir, 'figures-export');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
expfig = @(h, name) exportgraphics(h, fullfile(export_dir, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

% Tabbed figures → one PDF per tab.
tab_exports = {
    tg_fem2,        'fem2_irlsq_isaac_jun6';
    tg_fem2_bd,     'fem2_irlsq_isaac_jun6_bd';
    tg_fem2_rspec,  'fem2_rspec_isaac_jun6';
};
for g = 1:size(tab_exports, 1)
    tg     = tab_exports{g, 1};
    prefix = tab_exports{g, 2};
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

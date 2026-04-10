
% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));                % resize_fem_pair, resize_sparse_matrix, etc.

data_dir_isaac_apr7     = fullfile(script_dir, 'benchmark-output', 'isaac-apr7');
data_dir_isaac_apr9     = fullfile(script_dir, 'benchmark-output', 'isaac-apr9');
data_dir_photogrammetry = fullfile(script_dir, 'benchmark-output', 'photogrammetry2');

%% ========================================================================
%  FEM_Problem_2 — GEQP3 preconditioner path
%  Methods: CQRRT_linop, CholQR, sCholQR3, sCholQR3_basic, CQRRT_linop_stb
%  method_mask=47, upcast_orth=0, run_expl=0, 128 threads, d=2, nnz=4, b=256
%  Double/Float x Small/Medium/Large; double large merged from apr9 parts
%  ========================================================================

fig_geqp3    = figure('Name', 'GEQP3 Path — Applications Benchmarks', 'Position', [50, 50, 1400, 950]);
tg_geqp3     = uitabgroup(fig_geqp3);
fig_geqp3_bd = figure('Name', 'GEQP3 Path — Runtime Breakdowns',      'Position', [100, 50, 1800, 500]);
tg_geqp3_bd  = uitabgroup(fig_geqp3_bd);

% {timestamp_prefix, tab_label, data_dir}
datasets_geqp3 = {
    '20260407_181505', 'Double, Small (75860x4812)',       data_dir_isaac_apr7;
    '20260407_192748', 'Double, Medium (450000x28767)',    data_dir_isaac_apr7;
    'dbl_large_merged','Double, Large (600000x38286)',     data_dir_isaac_apr9;
    '20260407_161653', 'Float, Small (75860x4812)',        data_dir_isaac_apr7;
    '20260407_171702', 'Float, Medium (450000x28767)',     data_dir_isaac_apr7;
    '20260407_190904', 'Float, Large (600000x38286)',      data_dir_isaac_apr7;
};

for i = 1:size(datasets_geqp3, 1)
    ts    = datasets_geqp3{i, 1};
    label = datasets_geqp3{i, 2};
    ddir  = datasets_geqp3{i, 3};

    mt = uitab(tg_geqp3,    'Title', label);
    bt = uitab(tg_geqp3_bd, 'Title', label);

    res_csv = sprintf('%s_gsvd_results.csv',   ts);
    bd_csv  = sprintf('%s_gsvd_breakdown.csv', ts);
    plot_applications_results(ddir, res_csv, bd_csv, 'worst_ortho', label, mt, bt);
end

%% ========================================================================
%  FEM_Problem_2 — BQRRP preconditioner path
%  Methods: CQRRT_linop, CholQR, sCholQR3, sCholQR3_basic, CQRRT_linop_stb_bqrrp
%  Base mask=47 runs (excl. linop_stb) + mask=64 stb_bqrrp runs; merged in isaac-apr9
%  ========================================================================

fig_bqrrp    = figure('Name', 'BQRRP Path — Applications Benchmarks', 'Position', [50, 50, 1400, 950]);
tg_bqrrp     = uitabgroup(fig_bqrrp);
fig_bqrrp_bd = figure('Name', 'BQRRP Path — Runtime Breakdowns',      'Position', [100, 50, 1800, 500]);
tg_bqrrp_bd  = uitabgroup(fig_bqrrp_bd);

datasets_bqrrp = {
    'bqrrp_dbl_small',  'Double, Small (75860x4812)';
    'bqrrp_dbl_medium', 'Double, Medium (450000x28767)';
    'bqrrp_dbl_large',  'Double, Large (600000x38286)';
    'bqrrp_flt_small',  'Float, Small (75860x4812)';
    'bqrrp_flt_medium', 'Float, Medium (450000x28767)';
    'bqrrp_flt_large',  'Float, Large (600000x38286)';
};

for i = 1:size(datasets_bqrrp, 1)
    ts    = datasets_bqrrp{i, 1};
    label = datasets_bqrrp{i, 2};

    mt = uitab(tg_bqrrp,    'Title', label);
    bt = uitab(tg_bqrrp_bd, 'Title', label);

    res_csv = sprintf('%s_gsvd_results.csv',   ts);
    bd_csv  = sprintf('%s_gsvd_breakdown.csv', ts);
    plot_applications_results(data_dir_isaac_apr9, res_csv, bd_csv, 'worst_ortho', label, mt, bt);
end

%% ========================================================================
%  photogrammetry2 (4472 x 936) — local runs, all 6 methods + diagnostic
%  ========================================================================

fig_photo = figure('Name', 'photogrammetry2 Benchmarks', 'Position', [50, 50, 1400, 950]);
tg_photo  = uitabgroup(fig_photo);

fig_photo_bd = figure('Name', 'photogrammetry2 Breakdown', 'Position', [100, 50, 1800, 500]);
tg_photo_bd  = uitabgroup(fig_photo_bd);

% Applications benchmark (all 6 methods)
photo_mt = uitab(tg_photo,    'Title', 'photogrammetry2 (double)');
photo_bt = uitab(tg_photo_bd, 'Title', 'photogrammetry2 (double)');
plot_applications_results(data_dir_photogrammetry, ...
    '20260407_152645_gsvd_results.csv', ...
    '20260407_152645_gsvd_breakdown.csv', ...
    'worst_ortho', 'photogrammetry2 (double, all methods)', photo_mt, photo_bt);

% Diagnostic benchmark — two standalone figures (orthogonality + step-by-step)
diag_files = dir(fullfile(data_dir_photogrammetry, 'diagnostic_*.csv'));
if isempty(diag_files)
    error('No diagnostic CSV found in %s', data_dir_photogrammetry);
end
[~, diag_idx] = max([diag_files.datenum]);
diag_csv = fullfile(data_dir_photogrammetry, diag_files(diag_idx).name);
[diag_fig_orth, diag_fig_step] = plot_diagnostic_results(diag_csv, 'photogrammetry2 (double)');

%% ========================================================================
%  CQRRT linop vs explicit orthogonality gap (synthetic kappas + photogrammetry2)
%  Data: CQRRT_diagnostic generate mode, kappa=1e2,1e4,1e6,1e8 + photogrammetry2
%  ========================================================================

data_dir_orth_gap = fullfile(script_dir, 'benchmark-output', 'orth-gap');
fig_orth_gap = figure('Name', 'Orth Gap: linop vs expl', 'Position', [100 100 1000 580]);
plot_orth_gap(data_dir_orth_gap, fig_orth_gap);

%% ========================================================================
%  Export all figures
%  - Tabbed figures: each tab exported as <prefix>_<slug>.pdf
%  - Paper figures (matching \includegraphics in overleaf): explicit names
%  Output: figures-export/  (upload contents to overleaf figures/)
%  ========================================================================

export_dir = fullfile(script_dir, 'figures-export');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end

expfig = @(h, name) exportgraphics(h, fullfile(export_dir, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

% --- Export all tabs of each tabbed figure ---
% {uitabgroup handle, filename prefix}
tab_exports = {
    tg_geqp3,    'geqp3';
    tg_geqp3_bd, 'geqp3_bd';
    tg_bqrrp,    'bqrrp';
    tg_bqrrp_bd, 'bqrrp_bd';
    tg_photo,    'photo2';
    tg_photo_bd, 'photo2_bd';
};
for g = 1:size(tab_exports, 1)
    tg     = tab_exports{g, 1};
    prefix = tab_exports{g, 2};
    for k = 1:numel(tg.Children)
        tab = tg.Children(k);
        slug = lower(regexprep(tab.Title, '[^a-zA-Z0-9]+', '_'));
        slug = regexprep(slug, '^_+|_+$', '');
        fname = [prefix, '_', slug];
        expfig(tab, fname);
        fprintf('  Exported: %s.pdf\n', fname);
    end
end

% --- Diagnostic figures (paper names) ---
expfig(diag_fig_orth, 'cqrrt_orth_diagnostic');
fprintf('  Exported: cqrrt_orth_diagnostic.pdf\n');
expfig(diag_fig_step, 'cqrrt_stepwise_diagnostic');
fprintf('  Exported: cqrrt_stepwise_diagnostic.pdf\n');

% --- Orth gap (paper name, standalone figure) ---
expfig(fig_orth_gap, 'cqrrt_orth_gap');
fprintf('  Exported: cqrrt_orth_gap.pdf\n');

fprintf('All figures exported to %s\n', export_dir);

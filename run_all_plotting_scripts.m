% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));                % resize_fem_pair, resize_sparse_matrix, etc.

data_dir_irlsq          = fullfile(script_dir, 'benchmark-output', 'irlsq');
data_dir_isaac_apr7     = fullfile(script_dir, 'benchmark-output', 'isaac-apr7');     %#ok<NASGU>
data_dir_isaac_apr9     = fullfile(script_dir, 'benchmark-output', 'isaac-apr9');     %#ok<NASGU>
data_dir_photogrammetry = fullfile(script_dir, 'benchmark-output', 'photogrammetry2'); %#ok<NASGU>
data_dir_nmr            = fullfile(script_dir, 'benchmark-output', 'nmr');             %#ok<NASGU>

%% ========================================================================
%  PRtomo n=96  (parallel-beam X-ray CT, IRtools PRtomo, sm=true / sparse)
%
%  Active dataset is the larger n=96 run.  The n=64 noise-sweep data is
%  preserved in the LEGACY block below; uncomment it if you want the
%  log-log noise-sweep figure or a direct size comparison.
% =========================================================================

% {timestamp,  noise_level,  tab_label}
datasets_prtomo = {
    '20260525_122918',  0.05,   'PRtomo n=96 (24480x9216), noise=0.05';
};

% --- Tabbed figures (main + breakdown) ----------------------------------
fig_prtomo    = figure('Name', 'PRtomo n=96 — Sparse IR-LSQ',         'Position', [ 50, 50, 1200, 900]);
tg_prtomo     = uitabgroup(fig_prtomo);
fig_prtomo_bd = figure('Name', 'PRtomo n=96 — Runtime Breakdowns',    'Position', [100, 50, 1100, 600]);
tg_prtomo_bd  = uitabgroup(fig_prtomo_bd);

if exist(data_dir_irlsq, 'dir')
    for i = 1:size(datasets_prtomo, 1)
        ts    = datasets_prtomo{i, 1};
        label = datasets_prtomo{i, 3};
        if startsWith(ts, '__')
            fprintf('(skipping placeholder tab "%s")\n', label);
            continue;
        end
        mt = uitab(tg_prtomo,    'Title', label);
        bt = uitab(tg_prtomo_bd, 'Title', label);
        res_csv = sprintf('%s_irlsq_results.csv',   ts);
        bd_csv  = sprintf('%s_irlsq_breakdown.csv', ts);
        plot_irlsq_results(data_dir_irlsq, res_csv, bd_csv, label, mt, bt);
    end
else
    fprintf('(skipped PRtomo tabs — no data in %s yet)\n', data_dir_irlsq);
end

% Noise-sweep figure intentionally omitted at this size — we only have a
% single noise level (0.05) at n=96.  To restore the log-log noise-sweep
% plot, uncomment the n=64 noise-sweep block under "LEGACY" below; or
% generate a 5-point noise sweep at n=96 and reinstate the
% datasets_prtomo_noise_sweep / plot_noise_sweep call.
fig_prtomo_noise = [];

%% ========================================================================
%  Export
% =========================================================================

export_dir = fullfile(script_dir, 'figures-export');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
expfig = @(h, name) exportgraphics(h, fullfile(export_dir, [name, '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');

% Tabbed figures → one PDF per tab.
tab_exports = {
    tg_prtomo,    'prtomo_irlsq';
    tg_prtomo_bd, 'prtomo_irlsq_bd';
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

% Standalone noise-sweep figure
if ~isempty(fig_prtomo_noise) && isgraphics(fig_prtomo_noise)
    expfig(fig_prtomo_noise, 'prtomo_irlsq_noise_sweep');
    fprintf('  Exported: prtomo_irlsq_noise_sweep.pdf\n');
end

fprintf('All figures exported to %s\n', export_dir);


%% ========================================================================
%  LEGACY (commented out) — uncomment a block to re-render historical data.
% =========================================================================

% %% --- PRtomo n=64 noise sweep (5 noise levels, feeds plot_noise_sweep) --
% datasets_prtomo_noise_sweep = {
%     '20260525_122141',  0.0,    'PRtomo n=64, noise=0.0';
%     '20260525_121850',  0.001,  'PRtomo n=64, noise=0.001';
%     '20260525_123430',  0.01,   'PRtomo n=64, noise=0.01';
%     '20260525_120336',  0.05,   'PRtomo n=64, noise=0.05';
%     '20260525_123653',  0.1,    'PRtomo n=64, noise=0.1';
% };
% fig_prtomo_n64    = figure('Name', 'PRtomo n=64 — Sparse IR-LSQ',      'Position', [ 50, 50, 1200, 900]);
% tg_prtomo_n64     = uitabgroup(fig_prtomo_n64);
% fig_prtomo_n64_bd = figure('Name', 'PRtomo n=64 — Runtime Breakdowns', 'Position', [100, 50, 1100, 600]);
% tg_prtomo_n64_bd  = uitabgroup(fig_prtomo_n64_bd);
% for i = 1:size(datasets_prtomo_noise_sweep, 1)
%     ts    = datasets_prtomo_noise_sweep{i, 1};
%     label = datasets_prtomo_noise_sweep{i, 3};
%     mt = uitab(tg_prtomo_n64,    'Title', label);
%     bt = uitab(tg_prtomo_n64_bd, 'Title', label);
%     res_csv = sprintf('%s_irlsq_results.csv',   ts);
%     bd_csv  = sprintf('%s_irlsq_breakdown.csv', ts);
%     plot_irlsq_results(data_dir_irlsq, res_csv, bd_csv, label, mt, bt);
% end
% fig_prtomo_n64_noise = plot_noise_sweep(data_dir_irlsq, ...
%     datasets_prtomo_noise_sweep, 'PRtomo n=64 (5 algorithms overlay)');

% %% --- Photogrammetry2 IR-LSQ (sparse SpLinOp, retired pending PRtomo focus) -----
% datasets_photogrammetry_irlsq = {
%     '20260522_121857', 'photogrammetry2 (4472x936), 5 runs';
% };
% fig_photo_ir    = figure('Name', 'photogrammetry2 — Sparse IR-LSQ',      'Position', [ 50, 50, 1200, 900]);
% tg_photo_ir     = uitabgroup(fig_photo_ir);
% fig_photo_ir_bd = figure('Name', 'photogrammetry2 — Runtime Breakdowns', 'Position', [100, 50, 1100, 600]);
% tg_photo_ir_bd  = uitabgroup(fig_photo_ir_bd);
% for i = 1:size(datasets_photogrammetry_irlsq, 1)
%     ts    = datasets_photogrammetry_irlsq{i, 1};
%     label = datasets_photogrammetry_irlsq{i, 2};
%     mt = uitab(tg_photo_ir,    'Title', label);
%     bt = uitab(tg_photo_ir_bd, 'Title', label);
%     res_csv = sprintf('%s_irlsq_results.csv',   ts);
%     bd_csv  = sprintf('%s_irlsq_breakdown.csv', ts);
%     plot_irlsq_results(data_dir_irlsq, res_csv, bd_csv, label, mt, bt);
% end

% %% --- FEM_Problem_2 — GEQP3 preconditioner path -------------------------
% fig_geqp3    = figure('Name', 'GEQP3 Path — Applications Benchmarks', 'Position', [50, 50, 1400, 950]);
% tg_geqp3     = uitabgroup(fig_geqp3);
% fig_geqp3_bd = figure('Name', 'GEQP3 Path — Runtime Breakdowns',      'Position', [100, 50, 1400, 800]);
% tg_geqp3_bd  = uitabgroup(fig_geqp3_bd);
%
% datasets_geqp3 = {
%     '20260407_181505', 'Double, Small (75860x4812)',       data_dir_isaac_apr7;
%     '20260407_192748', 'Double, Medium (450000x28767)',    data_dir_isaac_apr7;
%     'dbl_large_merged','Double, Large (600000x38286)',     data_dir_isaac_apr9;
%     '20260407_161653', 'Float, Small (75860x4812)',        data_dir_isaac_apr7;
%     '20260407_171702', 'Float, Medium (450000x28767)',     data_dir_isaac_apr7;
%     '20260407_190904', 'Float, Large (600000x38286)',      data_dir_isaac_apr7;
% };
%
% for i = 1:size(datasets_geqp3, 1)
%     ts    = datasets_geqp3{i, 1};
%     label = datasets_geqp3{i, 2};
%     ddir  = datasets_geqp3{i, 3};
%     mt = uitab(tg_geqp3,    'Title', label);
%     bt = uitab(tg_geqp3_bd, 'Title', label);
%     res_csv = sprintf('%s_gsvd_results.csv',   ts);
%     bd_csv  = sprintf('%s_gsvd_breakdown.csv', ts);
%     plot_applications_results(ddir, res_csv, bd_csv, 'worst_ortho', label, mt, bt);
% end

% %% --- FEM_Problem_2 — BQRRP preconditioner path -------------------------
% fig_bqrrp    = figure('Name', 'BQRRP Path — Applications Benchmarks', 'Position', [50, 50, 1400, 950]);
% tg_bqrrp     = uitabgroup(fig_bqrrp);
% fig_bqrrp_bd = figure('Name', 'BQRRP Path — Runtime Breakdowns',      'Position', [100, 50, 1400, 800]);
% tg_bqrrp_bd  = uitabgroup(fig_bqrrp_bd);
%
% datasets_bqrrp = {
%     'bqrrp_dbl_small',  'Double, Small (75860x4812)';
%     'bqrrp_dbl_medium', 'Double, Medium (450000x28767)';
%     'bqrrp_dbl_large',  'Double, Large (600000x38286)';
%     'bqrrp_flt_small',  'Float, Small (75860x4812)';
%     'bqrrp_flt_medium', 'Float, Medium (450000x28767)';
%     'bqrrp_flt_large',  'Float, Large (600000x38286)';
% };
%
% for i = 1:size(datasets_bqrrp, 1)
%     ts    = datasets_bqrrp{i, 1};
%     label = datasets_bqrrp{i, 2};
%     mt = uitab(tg_bqrrp,    'Title', label);
%     bt = uitab(tg_bqrrp_bd, 'Title', label);
%     res_csv = sprintf('%s_gsvd_results.csv',   ts);
%     bd_csv  = sprintf('%s_gsvd_breakdown.csv', ts);
%     plot_applications_results(data_dir_isaac_apr9, res_csv, bd_csv, 'worst_ortho', label, mt, bt);
% end

% %% --- old photogrammetry2 (FEM applications benchmark + diagnostic) -----
% fig_photo = figure('Name', 'photogrammetry2 Benchmarks', 'Position', [50, 50, 1400, 950]);
% tg_photo  = uitabgroup(fig_photo);
% fig_photo_bd = figure('Name', 'photogrammetry2 Breakdown', 'Position', [100, 50, 1400, 800]);
% tg_photo_bd  = uitabgroup(fig_photo_bd);
%
% photo_mt = uitab(tg_photo,    'Title', 'photogrammetry2 (double)');
% photo_bt = uitab(tg_photo_bd, 'Title', 'photogrammetry2 (double)');
% plot_applications_results(data_dir_photogrammetry, ...
%     '20260407_152645_gsvd_results.csv', ...
%     '20260407_152645_gsvd_breakdown.csv', ...
%     'worst_ortho', 'photogrammetry2 (double, all methods)', photo_mt, photo_bt);
%
% diag_files = dir(fullfile(data_dir_photogrammetry, 'diagnostic_*.csv'));
% if isempty(diag_files)
%     error('No diagnostic CSV found in %s', data_dir_photogrammetry);
% end
% [~, diag_idx] = max([diag_files.datenum]);
% diag_csv = fullfile(data_dir_photogrammetry, diag_files(diag_idx).name);
% [diag_fig_orth, diag_fig_step] = plot_diagnostic_results(diag_csv, 'photogrammetry2 (double)');

% %% --- NMR Relaxometry (Kronecker operator, retired) ---------------------
% fig_nmr    = figure('Name', 'NMR Benchmark — IR-LSQ', 'Position', [50, 50, 1200, 900]);
% tg_nmr     = uitabgroup(fig_nmr);
% fig_nmr_bd = figure('Name', 'NMR Benchmark — Runtime Breakdowns', 'Position', [100, 50, 1100, 600]);
% tg_nmr_bd  = uitabgroup(fig_nmr_bd);
%
% datasets_nmr = {
%     '20260508_020045', 'Double, Small (16384x4096), λ=1.0';
%     '20260508_021844', 'Double, Medium (65536x16384), λ=1.0';
%     '20260508_022756', 'Double, Large (102400x25600), λ=1.0';
%     '20260507_171053', 'Double, Smoke n=16 (1024x256), λ=1.0';
% };
%
% if exist(data_dir_nmr, 'dir') && ~isempty(datasets_nmr)
%     for i = 1:size(datasets_nmr, 1)
%         ts    = datasets_nmr{i, 1};
%         label = datasets_nmr{i, 2};
%         mt = uitab(tg_nmr,    'Title', label);
%         bt = uitab(tg_nmr_bd, 'Title', label);
%         res_csv = sprintf('%s_nmr_results.csv',   ts);
%         bd_csv  = sprintf('%s_nmr_breakdown.csv', ts);
%         plot_nmr_results(data_dir_nmr, res_csv, bd_csv, label, mt, bt);
%     end
% end

% %% --- CQRRT linop-vs-explicit orthogonality gap -------------------------
% data_dir_orth_gap = fullfile(script_dir, 'benchmark-output', 'orth-gap');
% fig_orth_gap = figure('Name', 'Orth Gap: linop vs expl', 'Position', [100 100 1000 580]);
% plot_orth_gap(data_dir_orth_gap, fig_orth_gap);

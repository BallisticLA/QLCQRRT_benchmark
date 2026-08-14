% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % plotting helpers + count_comment_lines
addpath(fullfile(script_dir, 'utils'));

tab_groups = {};   % {tabgroup handle, export prefix}

%% ========================================================================
%  CAMPAIGN ERAS (rebased 2026-08-05, Max)
%
%  The era ladder rolled forward: the pre-07-29 campaigns and the [DIAG 07-29]
%  diagnostic cells are RETIRED from the figure set (their data stays on disk
%  under benchmark-output/{toeplitz_ls, toeplitz_ls_3case, irlsq_reg_1e6,
%  irlsq_reg_3kappa*, irlsq_diag}); the 07-31 warm/cold pair is now [OLD ...];
%  [NEW 08-05] is the rerun with the 08-05 policy changes.
%
%  Every era is marked in THREE places:
%    1. the figure Name,
%    2. the tab title (prefixed [OLD ...] / [NEW ...]),
%    3. the exported PDF filename prefix.
%  2026-07-30 (Max): era tags REMOVED from the in-plot super-titles -- the tab
%  strip and PDF filename carry the era; the panel title stays clean.
%
%  WHAT CHANGED in the [NEW 08-05] era (read before comparing to OLD):
%    * Warm-start policy: warm x_0 is Blendenpik-ONLY, and Blendenpik appears
%      TWICE per experiment ("Blendenpik" warm / "Blendenpik_cold"); every
%      Q-less QR method is cold. The OLD pair is the last of the
%      everything-warm / everything-cold design.
%    * Inner CG: single restart from the returned iterate against the TRUE
%      residual (IterRefineLSQ inner_restarts = 1); per-step cap 500 (was 1000).
%    * Timing hygiene: untimed CPU warmups before any measurement; num_runs = 5
%      recorded per method (run column); plotters aggregate via TIMING_AGG.
%    * Memory: PeakRSSTracker now malloc_trim's at start(), so peak-RSS bars
%      compare across methods regardless of execution order. OLD-era peak-RSS
%      bars over-report whichever method ran FIRST (CQRRT) and under-report the
%      rest -- do not read storage conclusions off the OLD figures.
%    * Wall-clock is NOT comparable across eras (policy + repeat changes).
% =========================================================================

ERA_OLD_WARM = '[OLD 07-31 warm]';   % warm x_0 in EVERY method, cap 1000
ERA_OLD_COLD = '[OLD 07-31 cold]';   % cold x_0 in EVERY method (Blendenpik too)
ERA_NEW      = '[NEW 08-05]';        % Blendenpik-only warm (both variants), CG restart, 5 runs

% Toeplitz campaigns: {data subdir, era tag, note for the plot title, ordered
%                      size list ({} = discover+sort), draw cross-size sweep tab}
TOEP_CAMPAIGNS = {
    'toeplitz_ls_warm', ERA_OLD_WARM, '4 fixed cases, warm x_0 all methods (superseded policy)', {'small','fixedm','middle','large'},  false
    'toeplitz_ls_cold', ERA_OLD_COLD, '4 fixed cases, cold x_0 all methods (Blendenpik too)',    {'small','fixedm','middle','large'},  false
    'toeplitz_ls_0805', ERA_NEW,      '4 fixed cases; Blendenpik warm+cold, 5 runs, CPU warmup', {'small','fixedm','middle','large'},  false
};

% FEM2 IR-LSQ campaigns: {data subdir, era tag, note, combos to look for}
FEM2_CAMPAIGNS = {
    'irlsq_reg_1e10',      ERA_OLD_WARM, 'warm x_0 all methods, cap 1000, stagnation exit; constructed: \kappa^{colnorm}\approx2e11 measured (1e10 target); ill = native, \kappa^{colnorm}\approx3e13 measured', {'dd'}
    'irlsq_reg_1e10_cold', ERA_OLD_COLD, 'cold x_0 all methods (Blendenpik too), cap 1000, stagnation exit -- the full-cold companion to the warm row', {'dd'}
    'irlsq_reg_0805',      ERA_NEW,      'Blendenpik warm+cold, all else cold; cap 500, CG single restart, 5 runs', {'dd'}
};

% 2026-08-05 (Max): the benchmarks now record num_runs repetitions per method
% (single-run timings at 4-7 solver iterations are at the noise floor).
% TIMING_AGG picks how multi-run CSVs collapse: 'best' (min-total run, default)
% or 'mean'. Single-run CSVs (every campaign before 08-05) are unaffected, and
% the aggregation is stamped into each figure title by aggregate_runs.
TIMING_AGG = 'best';

%% ========================================================================
%  Toeplitz least-squares benchmark (Oleg's 2nd experiment)
%
%  Prolate-Toeplitz regularized LS  min ||T x - b||^2 + lambda ||x||^2, augmented
%  A = [T; sqrt(lambda) I], solved by LSQR with Q-less QR right preconditioners.
%  6-panel plot_toeplitz_results layout: wall-time (build+solve), data & recovery
%  error, peak-vs-analytical memory, LSQR iterations, preconditioner orthogonality
%  loss, Cholesky shift retries. Methods: CQRRT / CholQR / CholQR2 / sCholQR3 /
%  sCholQR3_basic / Blendenpik / unpreconditioned.
% =========================================================================
for cc = 1:size(TOEP_CAMPAIGNS, 1)
    [sub, era, note, order, want_sweep] = deal(TOEP_CAMPAIGNS{cc, :});
    toep_dir = fullfile(script_dir, 'benchmark-output', sub);
    if ~exist(toep_dir, 'dir')
        fprintf('(skipped Toeplitz %s %s -- no data dir %s)\n', era, sub, toep_dir);
        continue;
    end
    if isempty(order)
        % Discover size subfolders, sorted by name (old dirs are zero-padded
        % n%05d, so name order == size order).
        d = dir(toep_dir);
        size_dirs = sort({d([d.isdir] & ~ismember({d.name}, {'.','..'})).name});
    else
        % Explicit order: the new tags (small/middle/large) do NOT sort into size
        % order alphabetically, and a mis-ordered x-axis reads as a real trend.
        size_dirs = order(cellfun(@(s) exist(fullfile(toep_dir, s), 'dir') == 7, order));
    end
    if isempty(size_dirs)
        fprintf('(skipped Toeplitz %s -- no size subfolders)\n', era); continue;
    end

    fig = figure('Name', sprintf('%s Toeplitz LS -- Q-less QR right preconditioners  (%s)', era, note), ...
                 'Position', [80 + 25*cc, 80, 1250, 900]);
    tg = uitabgroup(fig);
    if want_sweep
        st = uitab(tg, 'Title', sprintf('%s sweep', era));
        plot_toeplitz_sweep(toep_dir, size_dirs, st, TIMING_AGG);
    end
    any_size = false;
    for s = 1:numel(size_dirs)
        size_dir = fullfile(toep_dir, size_dirs{s});
        f = dir(fullfile(size_dir, '*_toeplitz_ls_results.csv'));
        f = f([f.bytes] > 0);                    % empty CSV = job died before writing
        if isempty(f)
            fprintf('(skipped Toeplitz %s %s -- no non-empty CSV)\n', era, size_dirs{s});
            continue;
        end
        [~, ord] = sort({f.name});               % timestamped names -> newest last
        res_csv = f(ord(end)).name;
        % Label the tab with the FULL m x n, not the folder name. The old folders are
        % named after n only (n%05d), so a tab reading "n16000" sits above a plot
        % titled "32000 x 16000" -- an easy misread of the tab as the matrix shape
        % (hit 2026-07-27). Read the shape out of the CSV so the two always agree.
        shape = toeplitz_shape_label(fullfile(size_dir, res_csv), size_dirs{s});
        mt = uitab(tg, 'Title', sprintf('%s %s', era, shape));
        plot_toeplitz_results(size_dir, res_csv, ...
            sprintf('prolate Toeplitz, LSQR -- %s  %s', shape, note), mt, TIMING_AGG);
        any_size = true;
    end
    if any_size
        tab_groups(end+1, :) = {tg, sprintf('toeplitz_%s', era_slug(era))}; %#ok<SAGROW>
    else
        close(fig);
    end
end

%% ========================================================================
%  FEM_Problem_2 App-1 IR-LSQ (regularized) + Blendenpik
%
%  min ||A x - b||^2 solved by IterRefineLSQ (2 outer refinement steps, inner CG on
%  the preconditioned normal equations) with Q-less QR right preconditioners, plus
%  the Blendenpik competitor (SASO sketch + Householder QR + LSQR). One 6-panel
%  results tab per cell. Runtime-breakdown tabs intentionally suppressed (empty
%  breakdown CSV name).
%
%  Cells are DISCOVERED from the directory rather than built from a
%  sizes x kappas x combos product: the 1e10 campaign includes `native_ill_dd`,
%  which has no kappa label and no size ladder (it is Oleg's unmodified generator
%  at exactly one size), so the old triple loop silently skipped it.
% =========================================================================
for cc = 1:size(FEM2_CAMPAIGNS, 1)
    [sub, era, note, combos] = deal(FEM2_CAMPAIGNS{cc, :});
    camp_dir = fullfile(script_dir, 'benchmark-output', sub);
    if ~exist(camp_dir, 'dir')
        fprintf('(skipped FEM2 %s %s -- no data dir)\n', era, sub); continue;
    end
    d = dir(camp_dir);
    cells = sort({d([d.isdir] & ~ismember({d.name}, {'.','..'})).name});
    for c = 1:numel(combos)
        combo = combos{c};
        % Cells for this precision combo, i.e. dirs ending _<combo>.
        mine = cells(endsWith(cells, ['_' combo]));
        if isempty(mine), continue; end
        fig = figure('Name', sprintf('%s FEM2 App 1 -- IR-LSQ (reg) -- %s  (%s)', era, combo, note), ...
                     'Position', [60 + 25*c + 40*cc, 60, 1150, 900]);
        tg = uitabgroup(fig);
        any_cell = false;
        for k = 1:numel(mine)
            cell_dir = fullfile(camp_dir, mine{k});
            f = dir(fullfile(cell_dir, '*_irlsq_reg_results.csv'));
            f = f([f.bytes] > 0);
            if isempty(f), fprintf('(missing cell: %s/%s)\n', sub, mine{k}); continue; end
            [~, ord] = sort({f.name});
            % NEWEST, consistently. This used to take f(1) (the OLDEST) for FEM2
            % while the Toeplitz path took the newest -- so a re-run silently did
            % not show up on one of the two plots.
            res_csv = f(ord(end)).name;
            label = strrep(mine{k}, '_', '\_');
            mt = uitab(tg, 'Title', sprintf('%s %s', era, mine{k}));
            plot_irlsq_results(cell_dir, res_csv, '', ...
                sprintf('%s   [%s]   %s', label, combo, note), mt, [], TIMING_AGG);
            any_cell = true;
        end
        if any_cell
            tab_groups(end+1, :) = {tg, sprintf('fem2_irlsq_%s_%s', era_slug(era), combo)}; %#ok<SAGROW>
        else
            close(fig);
        end
    end
end

% (The [DIAG 07-29] inner-CG diagnostic section was removed 2026-08-05 (Max):
%  the stagnation-detection fix it motivated has landed, its question is
%  answered in the 07-29/07-30 session logs, and its CSVs remain under
%  benchmark-output/irlsq_diag/ if ever needed again.)

%% ========================================================================
%  Export -- one vector PDF per tab. The era tag is already in the tab title, so
%  it survives into the slug; the prefix carries it too, so OLD and NEW PDFs can
%  never collide or be confused once separated from the figure window.
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
        % Strip the leading [ERA ...] tag before slugging: the prefix already
        % carries the era, and keeping both gave names like
        % toeplitz_OLD_old_07_11_07_15_sweep.pdf.
        title_bare = regexprep(tab.Title, '^\s*\[[^\]]*\]\s*', '');
        slug = lower(regexprep(title_bare, '[^a-zA-Z0-9]+', '_'));
        slug = regexprep(slug, '^_+|_+$', '');
        fname = [prefix, '_', slug];
        expfig(tab, fname);
        fprintf('  Exported: %s.pdf\n', fname);
    end
end
fprintf('All figures exported to %s\n', export_dir);

% ---- helpers -------------------------------------------------------------
function s = era_slug(era)
% '[OLD 07-11/07-15]' -> 'OLD'.  Keeps PDF prefixes short but unambiguous.
% warm/cold eras keep their qualifier so the two figure sets export to
% DISTINCT filenames instead of overwriting each other (2026-07-31).
    t = regexprep(era, '[\[\]]', '');
    s = regexp(t, '^[A-Z]+', 'match', 'once');
    if isempty(s), s = 'ERA'; end
    if contains(t, 'warm'), s = [s '_warm'];
    elseif contains(t, 'cold'), s = [s '_cold']; end
end

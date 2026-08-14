% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);                                  % repo root
addpath(fullfile(script_dir, 'utils'));
addpath(fullfile(script_dir, 'plotting'));                            % plot drivers

tab_groups = {};   % {tabgroup handle, export prefix}

%% ========================================================================
%  CAMPAIGN ERAS (rebased 2026-08-07, Max)
%
%  The era ladder rolled forward: every pre-08-06 campaign (the 07-31
%  warm/cold pairs, the 08-05 lsqr rows, the 3kappa/1e6/hard_kcolnorm sets,
%  irlsq_diag, the June isaac-* dirs) is RETIRED and its local data DELETED
%  (2026-08-07 cleanup; the raw output survives on ISAAC under
%  ~/slurm/CQRRT/benchmark-out/). The 08-06 campaign is now [OLD 08-06];
%  [NEW 08-07] is the rerun on the unified restarted_pcg_ne engine.
%
%  Every era is marked in THREE places:
%    1. the figure Name,
%    2. the tab title (prefixed [OLD ...] / [NEW ...]),
%    3. the exported PDF filename prefix.
%  2026-07-30 (Max): era tags REMOVED from the in-plot super-titles -- the tab
%  strip and PDF filename carry the era; the panel title stays clean.
%
%  WHAT CHANGED in the [NEW 08-10] era (RandLAPACK commits 4c32aae -> 60f366a
%  -> d148e73 -> b5a2886 -> 6d026d5 -> a73d03d; read before comparing to OLD):
%    * ONE solver engine: both benchmarks now run restarted_pcg_ne.
%      IterRefineLSQ is a thin adapter over it (bitwise-identical solutions
%      pinned by test); inner_restarts is GONE -- every round is a restart
%      from the returned iterate against the TRUE residual.
%    * Per-round restart pacing: restart_drop default 1e-2 -> 1e-4 (Oleg's
%      pacing), plus an inner_abs_tol guard so rounds stop once below the
%      absolute target. FEM2 CLI slot 13 is ir_round_drop now (>= 1 rejected).
%    * Round caps raised to 20 in BOTH benchmarks (FEM2 ir_n_steps, Toeplitz
%      pcg_max_restarts). The OLD 3-4 round caps truncated the
%      unpreconditioned baseline ~5 orders above tol -- its OLD-era accuracy
%      bars are cap artifacts, not solver quality.
%    * Outer stagnation exit: 2 rounds without true-residual improvement
%      stop the method, so noise-floor methods no longer grind the cap.
%    * BLAS-2 THREAD GUARD (d148e73), the reason 08-07 was rerun as 08-08:
%      MKL threads the preconditioner's triangular solves badly at 64 threads.
%      Calibrated on the benchmark node (Gold 6430): dtrsv is fastest at 8-16
%      threads and degrades past that, so the solves are now capped (8 for
%      n <= 4000, 16 above). The preconditioned methods apply the factor twice
%      per inner iteration and the unpreconditioned baseline not at all, so
%      the overhead taxed exactly the methods that converge fastest. Measured
%      0807 -> 0808 on solve wall-clock: preconditioned methods 1.3-2.0x
%      faster, unpreconditioned ~1.0x (unchanged, as expected).
%      => 08-07 SOLVE TIMINGS ARE RETIRED as thread-contaminated; its accuracy
%      and iteration counts were never affected.
%    * FFT THREAD CAP (b5a2886) -- THE ACTUAL FIX for the wall-clock inversion,
%      and the reason 08-08 was rerun as 08-09. MKL's threaded FFT
%      INTERMITTENTLY STALLS at 64 threads: single DftiComputeForward calls
%      were caught taking 16-32 ms instead of ~0.1 ms, ~99% of the stall inside
%      the FFT call (fill/multiply/backward stay at 80-100 us). A method making
%      276 applies averages those away; one making 11 cannot, so the artifact
%      scaled INVERSELY with preconditioner quality and inverted the ranking.
%      Thread sweep, 3 repeats: at 8 threads CQRRT solves in 19 ms against
%      unpreconditioned's 200 ms (10.6x, matching 5-vs-267 iterations); at 64
%      both read ~344 ms. The transform is now capped (default 16). 08-08 SOLVE
%      TIMINGS ARE THEREFORE ALSO RETIRED; 08-09 is the first trustworthy
%      wall-clock. Accuracy and iteration counts were never affected in any era.
%    * FOUR BLENDENPIK ROWS (6d026d5). As published Blendenpik is sketch + QR +
%      LSQR with no refinement, so comparing it against Q-less methods that all
%      run through refinement conflated preconditioner quality with solver
%      structure. The published warm/cold rows are unchanged; two rows now hand
%      Blendenpik's OWN R and answer to our restarted-PCG engine. Measured:
%      cold Blendenpik's recovery error 3.008e+03 becomes 4.107e-04 with
%      refinement at IDENTICAL preconditioner conditioning, so its failure is
%      solver structure, not a bad preconditioner.
%      (a73d03d, why 08-09 was rerun as 08-10: 6d026d5 added those rows only to
%      the sparse-mode selector, but run_irlsq_reg has its OWN selector, so the
%      08-09 FEM2 cells produced 7 rows instead of 9 despite method_mask=127.
%      The 08-09 Toeplitz data was complete and correct; both families were
%      rerun together so one era means one commit.)
%    * Wall-clock is NOT comparable across eras (round policy changed).
% =========================================================================

ERA_OLD = '[OLD 08-06]';   % 4 outer rounds + outer_tol early exit, 3 runs (see below)
ERA_NEW = '[NEW 08-10]';   % unified engine, FFT + BLAS-2 thread caps, 4 Blendenpik rows (see above)
% [OLD 08-06] era (commits 17a60e2 + 3b2ee61, campaign 2026-08-06):
%   * FEM2: ir_n_steps=4 outer refinement steps (was 2), outer_tol=10*eps early
%     exit, inner_restarts=1 -- healthy dd methods reach relres ~1.2e-16 in 3
%     steps (machine precision; the 08-05 era's 2-step runs plateaued ~2-3e-12).
%   * Toeplitz: solver = restarted_pcg_ne (Oleg's second solver) with the STABLE
%     round residual; all preconditioned methods sit at the 1.005e-10 noise floor
%     (the data's 1e-11 noise) with flag=1 at tol 1e-12 -- flag=1 here is NORMAL,
%     read relres. (Both eras ran Xeon Gold 6430 / 64 threads.)
TOEP_CAMPAIGNS = {
    'toeplitz_ls_0806_pcg_ne',  ERA_OLD, '4 fixed cases; restarted PCG-NE solver (stable residual), 3 runs', {'small','fixedm','middle','large'},  false
    'toeplitz_ls_0810_pcg_ne',  ERA_NEW, '4 fixed cases; unified PCG-NE engine, restart\_drop 1e-4, cap 50, stagnation exit, thread caps, 4 Blendenpik rows, 3 runs', {'small','fixedm','middle','large'},  false
};

% FEM2 IR-LSQ campaigns: {data subdir, era tag, note, combos to look for}
FEM2_CAMPAIGNS = {
    'irlsq_reg_0806',  ERA_OLD, '4 outer steps + 10\epsilon early exit, cap 500, 1 CG restart, 3 runs', {'dd'}
    'irlsq_reg_0810',  ERA_NEW, 'unified PCG-NE engine: round\_drop 1e-4, cap 50 rounds, stagnation exit, thread caps, 4 Blendenpik rows, 3 runs', {'dd'}
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
    toep_dir = fullfile(script_dir, 'results', sub);
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
    camp_dir = fullfile(script_dir, 'results', sub);
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
%  the stagnation-detection fix it motivated has landed and its question is
%  answered in the 07-29/07-30 session logs. Its CSVs were deleted locally in
%  the 2026-08-07 cleanup; ISAAC benchmark-out still has them if ever needed.)

%% ========================================================================
%  Export -- one vector PDF per tab. The era tag is already in the tab title, so
%  it survives into the slug; the prefix carries it too, so OLD and NEW PDFs can
%  never collide or be confused once separated from the figure window.
% =========================================================================
export_dir = fullfile(script_dir, 'figures');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end
% Wipe ALL prior PDFs so figures/ holds only the current run's tabs.
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

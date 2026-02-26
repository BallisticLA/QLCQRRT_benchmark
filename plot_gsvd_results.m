%% plot_gsvd_results.m — Plot GSVD / Generalized LS benchmark results
%
% Usage:
%   plot_gsvd_results(data_dir, results_csv, svals_csv)
%   plot_gsvd_results(data_dir, results_csv, svals_csv, plot_mode)
%
% Arguments:
%   data_dir    — directory containing the CSV files
%   results_csv — filename of *_gsvd_results.csv
%   svals_csv   — filename of *_gsvd_svals.csv
%   plot_mode   — 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%
% Produces three figures:
%   1. Timing bar chart — QR time + 3 application totals per algorithm
%   2. Orthogonality comparison — bar chart (log scale)
%   3. Singular value spectrum — overlay from each algorithm

function plot_gsvd_results(data_dir, results_csv, svals_csv, plot_mode)

if nargin < 4
    plot_mode = 'best_speed';
end

results_csv = fullfile(data_dir, results_csv);
svals_csv   = fullfile(data_dir, svals_csv);

%% ------------------------------------------------------------------
%  Load results CSV
%  ------------------------------------------------------------------
% Count comment lines (lines starting with #)
fid = fopen(results_csv, 'r');
n_comments = 0;
while true
    line = fgetl(fid);
    if line(1) == '#'
        n_comments = n_comments + 1;
    else
        break;
    end
end
fclose(fid);

opts = detectImportOptions(results_csv, 'NumHeaderLines', n_comments);
T = readtable(results_csv, opts);

algorithms = T.algorithm;
qr_time    = T.qr_time_us;
orth_error = T.orth_error;
max_orth   = T.max_orth_cols;
total_a    = T.total_a_time_us;
total_b    = T.total_b_time_us;
total_c    = T.total_c_time_us;
ls_err     = T.ls_rel_error;
m_val      = T.m(1);
n_val      = T.n(1);

unique_algs = unique(algorithms, 'stable');
n_algs = numel(unique_algs);

%% ------------------------------------------------------------------
%  Select representative run per algorithm based on plot_mode
%  ------------------------------------------------------------------
sel_idx = zeros(n_algs, 1);  % index into the full data arrays

for a = 1:n_algs
    mask = strcmp(algorithms, unique_algs{a});
    indices = find(mask);

    switch plot_mode
        case 'best_speed'
            [~, best] = min(total_a(mask));
            sel_idx(a) = indices(best);
        case 'worst_ortho'
            [~, worst] = max(orth_error(mask));
            sel_idx(a) = indices(worst);
        case 'best_ortho'
            [~, best] = min(orth_error(mask));
            sel_idx(a) = indices(best);
        otherwise
            error('Unknown plot_mode: %s', plot_mode);
    end
end

%% ------------------------------------------------------------------
%  Figure 1: Timing bar chart
%  ------------------------------------------------------------------
figure('Position', [100 100 800 500]);

bar_data = zeros(n_algs, 4);
for a = 1:n_algs
    i = sel_idx(a);
    bar_data(a, 1) = qr_time(i) / 1000;    % ms
    bar_data(a, 2) = total_a(i) / 1000;     % ms
    bar_data(a, 3) = total_b(i) / 1000;     % ms
    bar_data(a, 4) = total_c(i) / 1000;     % ms
end

b = bar(bar_data);
set(gca, 'XTickLabel', unique_algs);
ylabel('Time (ms)');
title(sprintf('GSVD Benchmark Timing (%d x %d, mode: %s)', m_val, n_val, plot_mode));
legend({'Q-less QR', 'App (a): Gen. LS', 'App (b): Gen. svals', 'App (c): Gen. svecs'}, ...
       'Location', 'northwest');
grid on;

% Color scheme
colors = [0.2 0.4 0.8;   % blue - QR
          0.8 0.3 0.2;   % red - App a
          0.3 0.7 0.3;   % green - App b
          0.7 0.4 0.8];  % purple - App c
for k = 1:4
    b(k).FaceColor = colors(k, :);
end

%% ------------------------------------------------------------------
%  Figure 2: Orthogonality comparison
%  ------------------------------------------------------------------
figure('Position', [100 650 600 400]);

orth_vals = zeros(n_algs, 1);
for a = 1:n_algs
    orth_vals(a) = orth_error(sel_idx(a));
end

bar(orth_vals);
set(gca, 'YScale', 'log');
set(gca, 'XTickLabel', unique_algs);
ylabel('||Q^T Q - I||_F / \surd{n}');
title(sprintf('Orthogonality (%d x %d, mode: %s)', m_val, n_val, plot_mode));
grid on;

% Add text labels on bars
for a = 1:n_algs
    text(a, orth_vals(a) * 2, sprintf('%.1e', orth_vals(a)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

%% ------------------------------------------------------------------
%  Figure 3: Singular value spectrum
%  ------------------------------------------------------------------
if nargin >= 3 && ~isempty(svals_csv)
    figure('Position', [750 650 600 400]);

    % Count comment lines in svals CSV
    fid2 = fopen(svals_csv, 'r');
    n_comments_sv = 0;
    while true
        line = fgetl(fid2);
        if line(1) == '#'
            n_comments_sv = n_comments_sv + 1;
        else
            break;
        end
    end
    fclose(fid2);

    sv_opts = detectImportOptions(svals_csv, 'NumHeaderLines', n_comments_sv);
    T_sv = readtable(svals_csv, sv_opts);

    sv_algs = T_sv.algorithm;
    % Singular value columns are sigma_0, sigma_1, ...
    sv_col_names = T_sv.Properties.VariableNames;
    sv_mask = startsWith(sv_col_names, 'sigma_');
    n_svals = sum(sv_mask);
    sv_col_idx = find(sv_mask);

    % Plot singular values for the selected run of each algorithm
    colors_line = lines(n_algs);
    markers = {'o', 's', '^', 'd'};
    hold on;
    for a = 1:n_algs
        mask = strcmp(sv_algs, unique_algs{a});
        sv_indices = find(mask);

        % Use the first matching run (all runs should give similar svals)
        if ~isempty(sv_indices)
            row = sv_indices(1);
            svals = zeros(1, n_svals);
            for j = 1:n_svals
                svals(j) = T_sv{row, sv_col_idx(j)};
            end
            semilogy(1:n_svals, svals, ['-', markers{mod(a-1, numel(markers))+1}], ...
                     'Color', colors_line(a, :), 'MarkerSize', 4, ...
                     'DisplayName', unique_algs{a});
        end
    end
    hold off;
    xlabel('Index i');
    ylabel('\sigma_i');
    title(sprintf('Generalized Singular Values (%d x %d)', m_val, n_val));
    legend('Location', 'northeast');
    grid on;
end

fprintf('Plots generated for %s\n', results_csv);
fprintf('  Mode: %s\n', plot_mode);
fprintf('  Algorithms: %s\n', strjoin(unique_algs, ', '));
fprintf('  Matrix size: %d x %d\n', m_val, n_val);

end

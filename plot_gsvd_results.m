%% plot_gsvd_results.m — Plot GSVD / Generalized LS benchmark results
%
% Usage:
%   plot_gsvd_results(results_csv, svals_csv)
%   plot_gsvd_results(results_csv, svals_csv, plot_mode)
%
% Arguments:
%   results_csv — path to *_gsvd_results.csv
%   svals_csv   — path to *_gsvd_svals.csv
%   plot_mode   — 'best_speed' (default), 'worst_ortho', or 'best_ortho'
%
% Produces three figures:
%   1. Timing bar chart — QR time + 3 application totals per algorithm
%   2. Orthogonality comparison — bar chart (log scale)
%   3. Singular value spectrum — overlay from each algorithm

function plot_gsvd_results(results_csv, svals_csv, plot_mode)

if nargin < 3
    plot_mode = 'best_speed';
end

%% ------------------------------------------------------------------
%  Load results CSV
%  ------------------------------------------------------------------
fid = fopen(results_csv, 'r');
% Skip comment lines
while true
    pos = ftell(fid);
    line = fgetl(fid);
    if line(1) ~= '#'
        fseek(fid, pos, 'bof');
        break;
    end
end
header = fgetl(fid);
cols = strsplit(header, ',');
n_cols = numel(cols);

% Read data
fmt = repmat('%f', 1, n_cols);
% Replace %f for algorithm column with %s
alg_col = find(strcmp(cols, 'algorithm'));
fmt_parts = cell(1, n_cols);
for i = 1:n_cols
    if i == alg_col
        fmt_parts{i} = '%s';
    else
        fmt_parts{i} = '%f';
    end
end
fseek(fid, 0, 'bof');
% Skip comments again
while true
    pos = ftell(fid);
    line = fgetl(fid);
    if line(1) ~= '#'
        fseek(fid, pos, 'bof');
        break;
    end
end
fgetl(fid);  % skip header
data = textscan(fid, strjoin(fmt_parts, ','), 'Delimiter', ',');
fclose(fid);

% Build column map
col_map = containers.Map();
for i = 1:n_cols
    col_map(cols{i}) = i;
end

algorithms = data{col_map('algorithm')};
qr_time    = data{col_map('qr_time_us')};
orth_error = data{col_map('orth_error')};
max_orth   = data{col_map('max_orth_cols')};
total_a    = data{col_map('total_a_time_us')};
total_b    = data{col_map('total_b_time_us')};
total_c    = data{col_map('total_c_time_us')};
ls_err     = data{col_map('ls_rel_error')};
m_val      = data{col_map('m')}(1);
n_val      = data{col_map('n')}(1);

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
if nargin >= 2 && ~isempty(svals_csv)
    figure('Position', [750 650 600 400]);

    fid2 = fopen(svals_csv, 'r');
    % Skip comment line
    while true
        pos = ftell(fid2);
        line = fgetl(fid2);
        if line(1) ~= '#'
            fseek(fid2, pos, 'bof');
            break;
        end
    end
    sv_header = fgetl(fid2);
    sv_cols = strsplit(sv_header, ',');
    n_sv_cols = numel(sv_cols);

    % Read: run, algorithm, sigma_0, sigma_1, ...
    sv_fmt = ['%f%s', repmat('%f', 1, n_sv_cols - 2)];
    sv_data = textscan(fid2, sv_fmt, 'Delimiter', ',');
    fclose(fid2);

    sv_algs = sv_data{2};
    n_svals = n_sv_cols - 2;

    % Plot singular values for the selected run of each algorithm
    colors_line = lines(n_algs);
    markers = {'o', 's', '^', 'd'};
    hold on;
    for a = 1:n_algs
        % Find the selected run for this algorithm in the svals data
        mask = strcmp(sv_algs, unique_algs{a});
        sv_indices = find(mask);

        % Use the first matching run (all runs should give similar svals)
        if ~isempty(sv_indices)
            row = sv_indices(1);
            svals = zeros(1, n_svals);
            for j = 1:n_svals
                svals(j) = sv_data{2 + j}(row);
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

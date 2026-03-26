%% plot_orth_gap.m — CQRRT linop vs explicit orthogonality gap visualization
%
% Reads CSV files from benchmark-output/orth-gap/ and produces two figures:
%   Figure 1: Orthogonality error comparison (linop vs expl) across matrices
%   Figure 2: Step-by-step diagnostic — relative difference at each algorithm step
%
% Usage: run from /home/mymel/data/QLCQRRT_benchmark/

clear; close all;
data_dir = fullfile(pwd, 'benchmark-output', 'orth-gap');

%% ===== Figure 1: Orthogonality gap across matrices =====
% Collect all orth_gap_*.csv files (excluding diag files)
files = dir(fullfile(data_dir, 'orth_gap_*.csv'));
files = files(~contains({files.name}, 'diag'));

matrix_names = {};
linop_orth   = [];
expl_orth    = [];
gap_ratios   = [];

for i = 1:numel(files)
    fpath = fullfile(data_dir, files(i).name);

    % Read header comment to extract matrix info
    fid = fopen(fpath, 'r');
    header1 = fgetl(fid); % # comment line 1
    header2 = fgetl(fid); % # comment line 2 with parameters
    fclose(fid);

    % Extract a display name
    fname = files(i).name;
    if contains(fname, '_gen_')
        % Synthetic: extract kappa from header
        kappa_match = regexp(header2, 'kappa=([0-9.e+]+)', 'tokens');
        if ~isempty(kappa_match)
            matrix_names{end+1} = sprintf('synthetic k=%.0e', str2double(kappa_match{1}{1}));
        else
            matrix_names{end+1} = 'synthetic';
        end
    else
        % File mode: extract matrix name from filename
        parts = regexp(fname, 'orth_gap_(.+?)_\d{8}', 'tokens');
        if ~isempty(parts)
            matrix_names{end+1} = strrep(parts{1}{1}, '_', '\_');
        else
            matrix_names{end+1} = fname;
        end
    end

    % Read data
    T = readtable(fpath, 'CommentStyle', '#');
    % Average over runs
    linop_orth(end+1)  = mean(T.linop_orth);
    expl_orth(end+1)   = mean(T.expl_orth);
    gap_ratios(end+1)  = mean(T.gap_ratio);
end

% Sort by gap ratio for better visualization
[~, idx] = sort(gap_ratios);
matrix_names = matrix_names(idx);
linop_orth   = linop_orth(idx);
expl_orth    = expl_orth(idx);
gap_ratios   = gap_ratios(idx);

n_matrices = numel(matrix_names);

figure('Position', [100 100 900 500]);

% Grouped bar chart on log scale
x = 1:n_matrices;
hold on;
b1 = bar(x - 0.15, log10(linop_orth), 0.3, 'FaceColor', [0.85 0.33 0.10]);
b2 = bar(x + 0.15, log10(expl_orth),  0.3, 'FaceColor', [0.00 0.45 0.74]);

% Add gap ratio labels above bars
for i = 1:n_matrices
    y_top = max(log10(linop_orth(i)), log10(expl_orth(i)));
    text(x(i), y_top + 0.5, sprintf('%.0fx', gap_ratios(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

hold off;
set(gca, 'XTick', x, 'XTickLabel', matrix_names, 'XTickLabelRotation', 30);
ylabel('log_{10}( ||Q^TQ - I||_F )');
title('CQRRT: Linop vs Explicit Orthogonality Error');
legend([b1 b2], {'CQRRT\_linop', 'CQRRT\_expl'}, 'Location', 'northwest');
grid on;

% Machine epsilon reference line
hold on;
yline(log10(eps), '--k', '\epsilon_{mach}', 'LabelHorizontalAlignment', 'left');
hold off;

%% ===== Figure 2: Step-by-step diagnostic =====
diag_files = dir(fullfile(data_dir, 'orth_gap_diag_*.csv'));
if isempty(diag_files)
    fprintf('No diagnostic CSV files found in %s\n', data_dir);
    return;
end

% Read the most recent diagnostic file
diag_path = fullfile(data_dir, diag_files(end).name);
fid = fopen(diag_path, 'r');
header1 = fgetl(fid); % # comment
header2 = fgetl(fid); % # parameters
fclose(fid);

T = readtable(diag_path, 'CommentStyle', '#');

% Filter to steps 1-7 (exclude orth_linop, orth_expl, orth_gap rows)
mask = T.step < 8;
steps = T.step(mask);
descriptions = T.description(mask);
rel_diffs = T.rel_diff(mask);

figure('Position', [100 650 800 400]);
bar_data = log10(rel_diffs);
b = barh(bar_data, 'FaceColor', [0.47 0.67 0.19]);

set(gca, 'YTick', 1:numel(descriptions), 'YTickLabel', descriptions);
xlabel('log_{10}( relative difference )');
title(sprintf('Step-by-step divergence: CQRRT\\_linop vs CQRRT\\_expl (%s)', ...
    strrep(header2(3:end), '_', '\_')));
grid on;
set(gca, 'YDir', 'reverse');  % Step 1 at top

% Add value labels on bars
for i = 1:numel(rel_diffs)
    text(bar_data(i) + 0.3, i, sprintf('%.1e', rel_diffs(i)), ...
        'VerticalAlignment', 'middle', 'FontSize', 8);
end

% Machine epsilon reference
hold on;
xline(log10(eps), '--r', '\epsilon_{mach}', 'LabelVerticalAlignment', 'bottom');
hold off;

fprintf('Done. Processed %d comparison files and %d diagnostic files.\n', ...
    numel(files), numel(diag_files));

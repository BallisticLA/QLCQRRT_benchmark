%% plot_orth_gap.m — CQRRT linop vs explicit orthogonality gap visualization
%
% Reads CSV files from benchmark-output/orth-gap/ and produces two figures:
%   Figure 1: Orthogonality error comparison (linop vs expl) across matrices
%   Figure 2: Step-by-step diagnostic — relative difference at each algorithm step
%
% Usage: run from /home/mymel/data/QLCQRRT_benchmark/

clear; close all;
data_dir = fullfile(pwd, 'benchmark-output', 'orth-gap');

% Colorblind-friendly palette (Wong 2011, Nature Methods)
color_linop = [0.90 0.60 0.00];  % orange
color_expl  = [0.00 0.45 0.70];  % blue
color_diag  = [0.00 0.62 0.45];  % teal

default_fontsize = 14;

%% ===== Figure 1: Orthogonality gap across matrices =====
% Collect all orth_gap_*.csv files (excluding diag files)
files = dir(fullfile(data_dir, 'orth_gap_*.csv'));
files = files(~contains({files.name}, 'diag'));

matrix_labels = {};
linop_orth    = [];
expl_orth     = [];
gap_ratios    = [];
precision_str = 'double';  % default; updated from CSV if available

for i = 1:numel(files)
    fpath = fullfile(data_dir, files(i).name);

    % Read header comments to extract matrix info
    fid = fopen(fpath, 'r');
    header1 = fgetl(fid);  % # comment line 1
    header2 = fgetl(fid);  % # comment line 2 with parameters
    fclose(fid);

    % Extract precision
    prec_match = regexp(header2, 'precision=(\w+)', 'tokens');
    if ~isempty(prec_match)
        precision_str = prec_match{1}{1};
    end

    % Extract display name (matrix name on line 1, kappa on line 2)
    fname = files(i).name;
    if contains(fname, '_gen_')
        % Synthetic: extract kappa and dimensions from header
        kappa_match = regexp(header2, 'kappa=([0-9.e+]+)', 'tokens');
        m_match = regexp(header2, 'm=(\d+)', 'tokens');
        n_match = regexp(header2, 'n=(\d+)', 'tokens');
        name_line = 'synthetic';
        if ~isempty(m_match) && ~isempty(n_match)
            name_line = sprintf('synthetic %sx%s', m_match{1}{1}, n_match{1}{1});
        end
        kappa_line = '';
        if ~isempty(kappa_match)
            kappa_line = sprintf('\\kappa = %.0e', str2double(kappa_match{1}{1}));
        end
        matrix_labels{end+1} = sprintf('%s\n%s', name_line, kappa_line);
    else
        % File mode: extract matrix name and dimensions
        parts = regexp(fname, 'orth_gap_(.+?)_\d{8}', 'tokens');
        m_match = regexp(header2, 'm=(\d+)', 'tokens');
        n_match = regexp(header2, 'n=(\d+)', 'tokens');
        if ~isempty(parts)
            name_str = parts{1}{1};
        else
            name_str = fname;
        end
        dim_str = '';
        if ~isempty(m_match) && ~isempty(n_match)
            dim_str = sprintf('%sx%s', m_match{1}{1}, n_match{1}{1});
        end
        matrix_labels{end+1} = sprintf('%s %s', name_str, dim_str);
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
matrix_labels = matrix_labels(idx);
linop_orth    = linop_orth(idx);
expl_orth     = expl_orth(idx);
gap_ratios    = gap_ratios(idx);

n_matrices = numel(matrix_labels);

% Machine epsilon depends on precision
if strcmp(precision_str, 'float')
    eps_val = single(eps('single'));
else
    eps_val = eps('double');
end

figure('Position', [100 100 1000 550]);

x = 1:n_matrices;
hold on;
b1 = bar(x - 0.17, log10(linop_orth), 0.3, 'FaceColor', color_linop, 'EdgeColor', 'none');
b2 = bar(x + 0.17, log10(expl_orth),  0.3, 'FaceColor', color_expl,  'EdgeColor', 'none');
yline(log10(eps_val), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'Label', sprintf('\\epsilon_{mach} (%s)', precision_str), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', default_fontsize - 2);
hold off;

set(gca, 'XTick', x, 'XTickLabel', matrix_labels, 'FontSize', default_fontsize);
ylabel('log_{10}( ||Q^TQ - I||_F )', 'FontSize', default_fontsize);
title('CQRRT: LinOp vs Explicit Orthogonality Error', 'FontSize', default_fontsize + 2);
legend([b1 b2], {'CQRRT\_linop', 'CQRRT\_expl'}, ...
    'Location', 'northwest', 'FontSize', default_fontsize);
grid on;
set(gca, 'TickLabelInterpreter', 'tex');

%% ===== Figure 2: Step-by-step diagnostic =====
diag_files = dir(fullfile(data_dir, 'orth_gap_diag_*.csv'));
if isempty(diag_files)
    fprintf('No diagnostic CSV files found in %s\n', data_dir);
    return;
end

% Read the most recent diagnostic file
diag_path = fullfile(data_dir, diag_files(end).name);
fid = fopen(diag_path, 'r');
header1 = fgetl(fid);  % # comment
header2 = fgetl(fid);  % # parameters
fclose(fid);

% Extract precision from diagnostic
prec_match = regexp(header2, 'precision=(\w+)', 'tokens');
if ~isempty(prec_match)
    precision_str = prec_match{1}{1};
end
if strcmp(precision_str, 'float')
    eps_val = single(eps('single'));
else
    eps_val = eps('double');
end

% Extract m, n for title
m_match = regexp(header2, 'm=(\d+)', 'tokens');
n_match = regexp(header2, 'n=(\d+)', 'tokens');
dim_str = '';
if ~isempty(m_match) && ~isempty(n_match)
    dim_str = sprintf('%s x %s', m_match{1}{1}, n_match{1}{1});
end

T = readtable(diag_path, 'CommentStyle', '#');

% Filter to steps 1-7 (exclude orth_linop, orth_expl, orth_gap rows)
mask = T.step < 8;
descriptions = T.description(mask);
rel_diffs = T.rel_diff(mask);

% Map CSV description keys to paper-consistent labels
label_map = containers.Map();
label_map('sketch')      = {'Sketch', 'M^{sk} = SM'};
label_map('qr_R_sk')     = {'QR', 'R^{sk}'};
label_map('R_pre')        = {'Preconditioner', 'R^{pre} = (R^{sk})^{-1}'};
label_map('A_pre')        = {'Preconditioning', 'MR^{pre}'};
label_map('AtApre_gemm')  = {'A^T \cdot A_{pre}', 'M^T(MR^{pre})'};
label_map('gram')         = {'Gram', '(R^{pre})^T M^T M R^{pre}'};
label_map('cholesky')     = {'Cholesky', 'R^{chol}'};
label_map('final_R')      = {'Final R', 'R^{chol} R^{sk}'};

step_labels = cell(numel(descriptions), 1);
for i = 1:numel(descriptions)
    key = descriptions{i};
    if isKey(label_map, key)
        pair = label_map(key);
        step_labels{i} = sprintf('%s\n%s', pair{1}, pair{2});
    else
        step_labels{i} = strrep(key, '_', ' ');
    end
end

figure('Position', [100 700 1000 500]);
bar_data = log10(rel_diffs);
b = bar(bar_data, 'FaceColor', color_diag, 'EdgeColor', 'none', 'BarWidth', 0.6);

set(gca, 'XTick', 1:numel(step_labels), 'XTickLabel', step_labels, ...
    'FontSize', default_fontsize, 'TickLabelInterpreter', 'tex');
ylabel('log_{10}( relative difference )', 'FontSize', default_fontsize);

% Build title with matrix name
% Try to infer matrix name from the file that produced the diagnostic
% (diag mode is always run on a specific file; name embedded in header or nearby comparison CSVs)
title_str = sprintf('CQRRT\\_linop vs CQRRT\\_expl: step-by-step divergence (%s, %s)', ...
    dim_str, precision_str);
% Check if photogrammetry2 comparison file exists to get matrix name
comp_files = dir(fullfile(data_dir, 'orth_gap_photogrammetry2_*.csv'));
if ~isempty(comp_files)
    title_str = sprintf('CQRRT\\_linop vs CQRRT\\_expl: step-by-step divergence\nphotogrammetry2 (%s, %s)', ...
        dim_str, precision_str);
end
title(title_str, 'FontSize', default_fontsize + 2);
grid on;

% Add value labels above bars
for i = 1:numel(rel_diffs)
    y_pos = bar_data(i);
    text(i, y_pos + 0.4, sprintf('%.1e', rel_diffs(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', default_fontsize - 2);
end

% Machine epsilon reference line
hold on;
yline(log10(eps_val), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
    'Label', sprintf('\\epsilon_{mach} (%s)', precision_str), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', default_fontsize - 2);
hold off;

fprintf('Done. Processed %d comparison files and %d diagnostic files.\n', ...
    numel(files), numel(diag_files));

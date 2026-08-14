function label = toeplitz_shape_label(csv_path, fallback)
%TOEPLITZ_SHAPE_LABEL  Human-readable "m x n" label for a Toeplitz-LS result CSV.
%
%   label = toeplitz_shape_label(csv_path, fallback)
%
% The per-size output folders are named after n only (n%05d), so a tab or caption
% built from the folder name reads as "n16000" while the underlying experiment is
% 32000 x 16000. That mismatch caused a real misreading of the figures on
% 2026-07-27. This helper reads the actual shape out of the CSV so labels and data
% can never disagree.
%
% Falls back to the supplied string if the CSV is missing, empty, or lacks the
% m/n columns, so callers never break on a partial campaign.

    label = fallback;
    if ~isfile(csv_path), return; end
    d = dir(csv_path);
    if isempty(d) || d.bytes == 0, return; end

    try
        n_skip = count_comment_lines(csv_path);
        opts = detectImportOptions(csv_path, 'NumHeaderLines', n_skip);
        T = readtable(csv_path, opts);
        if isempty(T) || ~all(ismember({'m','n'}, T.Properties.VariableNames))
            return;
        end
        label = sprintf('%d x %d', T.m(1), T.n(1));
    catch
        % Malformed CSV: keep the fallback rather than aborting the whole figure.
    end
end

%% resize_fem_pair.m — Generate a resized (K, V) FEM matrix pair
%
% Usage:
%   resize_fem_pair(K_file, V_file, output_dir, 'shrink', p)
%       Extract the p×p principal submatrix of K and first p rows of V.
%       Drops V columns that become entirely zero after the row cut.
%       Output: <base>_shrink<p>_K.mtx, <base>_shrink<p>_V.mtx
%
%   resize_fem_pair(K_file, V_file, output_dir, 'expand', k)
%       Stack k block-diagonal copies: K_big = kron(I_k, K), V_big = kron(I_k, V).
%       kappa(L^{-1}V) is exactly preserved by construction.
%       Output: <base>_expand<k>_K.mtx, <base>_expand<k>_V.mtx
%
%   <base> is derived from K_file by stripping the trailing _K.mtx suffix.
%
% Example:
%   dir = '/home/mymel/data/QLCQRRT_benchmark/input_matrices';
%   resize_fem_pair(fullfile(dir,'FEM_Problem_2_K.mtx'), ...
%                   fullfile(dir,'FEM_Problem_2_V.mtx'), dir, 'shrink', 30000);
%   resize_fem_pair(fullfile(dir,'FEM_Problem_2_K.mtx'), ...
%                   fullfile(dir,'FEM_Problem_2_V.mtx'), dir, 'expand', 3);

function resize_fem_pair(K_file, V_file, output_dir, mode, param)

if ~ismember(mode, {'shrink', 'expand'})
    error('resize_fem_pair: mode must be ''shrink'' or ''expand''');
end
if param ~= floor(param) || param <= 0
    error('resize_fem_pair: param must be a positive integer');
end

% --- read inputs ---
fprintf('Reading %s ... ', K_file);
K = read_mtx_sym(K_file);
fprintf('(%d x %d, nnz=%d)\n', size(K,1), size(K,2), nnz(K));

fprintf('Reading %s ... ', V_file);
V = read_mtx_gen(V_file);
fprintf('(%d x %d, nnz=%d)\n', size(V,1), size(V,2), nnz(V));

[m, n] = size(V);
assert(size(K,1) == m, 'K row count (%d) does not match V row count (%d)', size(K,1), m);

% --- apply operation ---
switch mode
    case 'shrink'
        p = param;
        if p >= m
            error('shrink size p=%d must be less than m=%d', p, m);
        end
        fprintf('\nShrinking to %d x %d ...\n', p, p);
        K_new = K(1:p, 1:p);
        V_new = V(1:p, :);
        keep  = any(V_new, 1);
        V_new = V_new(:, keep);
        n_new = sum(keep);
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), p, p, nnz(K_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d] (%d zero cols dropped)\n', ...
            m, n, nnz(V), p, n_new, nnz(V_new), n - n_new);
        tag = sprintf('shrink%d', p);

    case 'expand'
        k = param;
        fprintf('\nExpanding by factor %d ...\n', k);
        K_new = kron(speye(k), K);
        V_new = kron(speye(k), V);
        fprintf('  K: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, m, nnz(K), k*m, k*m, nnz(K_new));
        fprintf('  V: [%d x %d, nnz=%d] -> [%d x %d, nnz=%d]\n', ...
            m, n, nnz(V), k*m, k*n, nnz(V_new));
        fprintf('  Note: kappa(L^{-1}V) is exactly preserved.\n');
        tag = sprintf('expand%d', k);
end

% --- derive output filenames ---
[~, K_name, ~] = fileparts(K_file);
if length(K_name) > 2 && strcmp(K_name(end-1:end), '_K')
    base = K_name(1:end-2);
else
    base = K_name;
end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
K_out = fullfile(output_dir, sprintf('%s_%s_K.mtx', base, tag));
V_out = fullfile(output_dir, sprintf('%s_%s_V.mtx', base, tag));

% --- write outputs ---
fprintf('\nWriting %s ...\n', K_out);
write_mtx_sym(K_out, K_new);
fprintf('Writing %s ...\n', V_out);
write_mtx_gen(V_out, V_new);
fprintf('Done.\n');

end


% =========================================================================
% Local helpers
% =========================================================================

function A = read_mtx_sym(fname)
% Read a symmetric Matrix Market coordinate file.
% Reconstructs the full symmetric matrix from the stored lower-triangle entries.
    [rows, cols, vals, m, n] = read_mtx_coords(fname);
    A = sparse(rows, cols, vals, m, n);
    mask = rows ~= cols;
    A = A + sparse(cols(mask), rows(mask), vals(mask), m, n);
end

function A = read_mtx_gen(fname)
% Read a general (non-symmetric) Matrix Market coordinate file.
    [rows, cols, vals, m, n] = read_mtx_coords(fname);
    A = sparse(rows, cols, vals, m, n);
end

function [rows, cols, vals, m, n] = read_mtx_coords(fname)
% Core reader: skips % comment lines, parses size header, reads coordinate data.
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open: %s', fname); end
    line = fgetl(fid);
    while ~isempty(line) && line(1) == '%'
        line = fgetl(fid);
    end
    dims       = sscanf(line, '%d %d %d');
    m          = dims(1);
    n          = dims(2);
    nnz_count  = dims(3);
    data       = fscanf(fid, '%d %d %lf', [3, nnz_count]);
    fclose(fid);
    rows = data(1,:)';
    cols = data(2,:)';
    vals = data(3,:)';
end

function write_mtx_sym(fname, A)
% Write a sparse symmetric matrix in Matrix Market format (lower triangle only).
    [rows, cols, vals] = find(tril(A));
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real symmetric\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end

function write_mtx_gen(fname, A)
% Write a general sparse matrix in Matrix Market format.
    [rows, cols, vals] = find(A);
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end

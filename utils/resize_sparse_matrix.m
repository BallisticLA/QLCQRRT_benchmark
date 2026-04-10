%% resize_sparse_matrix.m — Generate upscaled versions of a single sparse matrix
%
% Applies kron(I_k, A) to produce k copies of A in a block-diagonal layout,
% scaling both m and n by k while exactly preserving the condition number.
%
% Usage:
%   resize_sparse_matrix(mtx_file, output_dir, k_values)
%
% Arguments:
%   mtx_file   — path to the base Matrix Market file (general, m x n, m >= n)
%   output_dir — directory to write output .mtx files
%   k_values   — vector of expansion factors, e.g. [2 4 8]
%
% Output files are named <base>_<m_out>x<n_out>.mtx
% (k=1 is the base matrix itself — not written, it already exists)
%
% Note: kron(I_k, A) has dimensions (k*m) x (k*n) and condition number
% exactly equal to that of A (singular values are k copies of each of A's).
%
% Example:
%   dir = '/home/mymel/matlab/QLCQRRT_benchmark/input_matrices/photogrammetry2';
%   resize_sparse_matrix(fullfile(dir, 'photogrammetry2.mtx'), dir, [2 4 8]);

function resize_sparse_matrix(mtx_file, output_dir, k_values)

fprintf('Reading %s ...\n', mtx_file);
A = read_mtx_gen(mtx_file);
[m, n] = size(A);
fprintf('  Base matrix: %d x %d, nnz=%d\n', m, n, nnz(A));

if m < n
    error('resize_sparse_matrix: matrix must be overdetermined (m >= n), got %dx%d', m, n);
end

[~, base_name, ~] = fileparts(mtx_file);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

for k = k_values(:)'
    A_big = kron(speye(k), A);
    m_out = k * m;
    n_out = k * n;
    out_file = fullfile(output_dir, sprintf('%s_%dx%d.mtx', base_name, m_out, n_out));
    fprintf('  k=%d: [%d x %d, nnz=%d] -> writing %s ...\n', ...
        k, m_out, n_out, nnz(A_big), out_file);
    write_mtx_gen(out_file, A_big);
end

fprintf('Done.\n');
end


% =========================================================================
% Local helpers
% =========================================================================

function A = read_mtx_gen(fname)
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open: %s', fname); end
    line = fgetl(fid);
    while ~isempty(line) && line(1) == '%'
        line = fgetl(fid);
    end
    dims     = sscanf(line, '%d %d %d');
    m        = dims(1); n = dims(2); nz = dims(3);
    data     = fscanf(fid, '%d %d %lf', [3, nz]);
    fclose(fid);
    A = sparse(data(1,:)', data(2,:)', data(3,:)', m, n);
end

function write_mtx_gen(fname, A)
    [rows, cols, vals] = find(A);
    fid = fopen(fname, 'w');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
    fprintf(fid, '%d %d %d\n', size(A,1), size(A,2), numel(vals));
    fprintf(fid, '%d %d %.15E\n', [rows(:)'; cols(:)'; vals(:)']);
    fclose(fid);
end

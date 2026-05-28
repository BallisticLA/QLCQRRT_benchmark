% trsm_direction_test.m
%
% Investigate the collaborator's claim about CQRRT TRSM_IDENTITY stability.
%
% Two ways to form R_sk^{-1} explicitly from R_sk:
%   M1 ("R_pre * R_sk = I"): R_inv = eye(n) / R_sk    (= (R_sk' \ eye)')
%                            <-- this is what current C++ does (Side::Right TRSM)
%   M2 ("R_sk * R_pre = I"): R_inv = R_sk \ eye(n)
%                            <-- proposed fix (Side::Left TRSM)
%
% Both are mathematically R_sk^{-1}, but the backward error of the
% subsequent product A * R_sk^{-1} differs.
%
% For each method we run a full sketch-preconditioned Cholesky-QR:
%   1. sketch S*A, R_sk = qr factor
%   2. form R_sk^{-1} via method M1 or M2
%   3. A_pre = A * R_sk^{-1}
%   4. G = A_pre' * A_pre
%   5. R_chol = chol(G) (may fail if G is not numerically SPD)
%   6. R = R_chol * R_sk
%   7. measure ||Q'Q - I||, ||A - Q R||, etc., where Q = A / R
%
% We also include a "direct" baseline that does NOT form R_sk^{-1} explicitly,
% mimicking the non-linop CQRRT path: A_pre = A / R_sk via one TRSM.
%
% Plus a "qrcp" stabilized baseline analogous to the GEQP3/BQRRP path in C++.

clear; clc;

%% Load matrix
mtx_path = 'input_matrices/photogrammetry2/photogrammetry2.mtx';
fprintf('Loading %s ...\n', mtx_path);
A = read_mtx_coord_real(mtx_path);
A = full(A);
[m, n] = size(A);
fprintf('  A is %d x %d\n', m, n);

%% Sketch
rng(42);
d = 2 * n;                  % embedding factor 2
S = randn(d, m) / sqrt(d);  % dense Gaussian sketch (simpler than SASO for matlab)
SA = S * A;
[~, R_sk] = qr(SA, 0);      % thin QR of sketch -> n x n upper triangular
% Force exact upper triangular (defensive)
R_sk = triu(R_sk);

cA   = cond(A);
cRsk = cond(R_sk);
fprintf('\ncond(A)    = %.3e\n', cA);
fprintf('cond(R_sk) = %.3e\n', cRsk);

In = eye(n);

%% Method 1: Side::Right TRSM (current C++)  ->  R_pre * R_sk = I
R_inv_M1 = In / R_sk;            % == (R_sk' \ In)'
err_inv_M1_right = norm(R_inv_M1 * R_sk - In, 'fro');
err_inv_M1_left  = norm(R_sk * R_inv_M1 - In, 'fro');

%% Method 2: Side::Left TRSM (proposed fix) ->  R_sk * R_pre = I
R_inv_M2 = R_sk \ In;
err_inv_M2_right = norm(R_inv_M2 * R_sk - In, 'fro');
err_inv_M2_left  = norm(R_sk * R_inv_M2 - In, 'fro');

fprintf('\n--- R_sk^{-1} forward errors ---\n');
fprintf('  M1 (Right TRSM, eye/R_sk):  ||R_inv*R_sk - I||=%.3e  ||R_sk*R_inv - I||=%.3e\n', ...
    err_inv_M1_right, err_inv_M1_left);
fprintf('  M2 (Left  TRSM, R_sk\\eye):  ||R_inv*R_sk - I||=%.3e  ||R_sk*R_inv - I||=%.3e\n', ...
    err_inv_M2_right, err_inv_M2_left);

%% Run full sketch-preconditioned Cholesky-QR for each variant
variants = {
    'M1_right_TRSM_eye',   @(A,Rsk) A * (eye(n) / Rsk);     % current C++
    'M2_left_TRSM_eye',    @(A,Rsk) A * (Rsk \ eye(n));     % proposed fix
    'direct_one_shot',     @(A,Rsk) A / Rsk;                % rl_cqrrt.hh style
};

fprintf('\n--- Sketch-preconditioned Cholesky QR results ---\n');
fprintf('%-22s | %-11s | %-11s | %-11s | %-11s\n', ...
    'variant', 'chol_ok', '||Q^TQ-I||', '||A-QR||/||A||', 'cond(R_final)');
fprintf('%s\n', repmat('-', 1, 80));

for k = 1:size(variants,1)
    name = variants{k,1};
    formA_pre = variants{k,2};

    A_pre = formA_pre(A, R_sk);

    G = A_pre' * A_pre;
    % symmetrize
    G = 0.5 * (G + G');

    [R_chol, p] = chol(G);
    chol_ok = (p == 0);
    if ~chol_ok
        fprintf('%-22s | %-11s | %-11s | %-11s | %-11s\n', ...
            name, sprintf('FAIL@%d',p), '-', '-', '-');
        continue;
    end

    R_final = R_chol * R_sk;
    cRf = cond(R_final);

    Q = A / R_final;                 % implicit Q
    orth = norm(Q' * Q - eye(n), 'fro');
    res  = norm(A - Q * R_final, 'fro') / norm(A, 'fro');

    fprintf('%-22s | %-11s | %-11.3e | %-11.3e | %-11.3e\n', ...
        name, 'yes', orth, res, cRf);
end

%% QRCP-based stabilized variant (analog of GEQP3 path in C++)
[Q_qrcp, R_qrcp, p] = qr(R_sk, 'vector');   % R_sk = Q*R*P^T (column pivot)
P = zeros(n);
for j = 1:n, P(p(j), j) = 1; end
% R_sk^{-1} = P * R_qrcp^{-1} * Q_qrcp^T
R_inv_qrcp = P * (R_qrcp \ Q_qrcp');
A_pre = A * R_inv_qrcp;
G = 0.5 * (A_pre' * A_pre + (A_pre' * A_pre)');
[R_chol, pp] = chol(G);
if pp == 0
    R_final = R_chol * R_sk;
    Q = A / R_final;
    orth = norm(Q' * Q - eye(n), 'fro');
    res  = norm(A - Q * R_final, 'fro') / norm(A, 'fro');
    fprintf('%-22s | %-11s | %-11.3e | %-11.3e | %-11.3e\n', ...
        'qrcp_stabilized', 'yes', orth, res, cond(R_final));
else
    fprintf('%-22s | %-11s | %-11s | %-11s | %-11s\n', ...
        'qrcp_stabilized', sprintf('FAIL@%d',pp), '-', '-', '-');
end

fprintf('\nDone.\n');

%% ------------------------------------------------------------------
function A = read_mtx_coord_real(fname)
    fid = fopen(fname, 'r');
    if fid < 0, error('cannot open %s', fname); end
    cleanup = onCleanup(@() fclose(fid));
    header = fgetl(fid);
    if isempty(strfind(header, 'coordinate')) || isempty(strfind(header, 'real'))
        error('expected "coordinate real" mtx, got: %s', header);
    end
    % skip comments and blank lines
    while true
        ln = fgetl(fid);
        if ischar(ln) && ~isempty(ln) && ln(1) ~= '%'
            break;
        end
    end
    sizes = sscanf(ln, '%d %d %d');
    M = sizes(1); N = sizes(2); nnz_ = sizes(3);
    data = fscanf(fid, '%d %d %f', [3, nnz_]);
    A = sparse(data(1,:), data(2,:), data(3,:), M, N);
end


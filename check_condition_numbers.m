%% Compare condition numbers: cond(K^{-1}V), cond(L^{-1}V), cond(L1^{-1}V)
%  L  = Cholesky factor from .mat file (permuted, non-triangular)
%  L1 = MATLAB's own chol(K) (proper triangular + permutation)
%
% Usage: run from /home/mymel/data/QLCQRRT_benchmark/

data_dir = 'input_matrices';

%% Load matrices from .mat file
fprintf('Loading FEM_Problem_2.mat...\n');
tic;
S = load(fullfile(data_dir, 'FEM_Problem_2.mat'));
P = S.Problem;
K = P.K;
V = P.V;
L = P.L;
clear S P;
fprintf('  Loaded in %.1fs\n', toc);
[m, n] = size(V);
fprintf('  K: %d x %d, nnz=%d\n', size(K,1), size(K,2), nnz(K));
fprintf('  V: %d x %d, nnz=%d\n', m, n, nnz(V));
fprintf('  L: %d x %d, nnz=%d\n', size(L,1), size(L,2), nnz(L));

%% Check L from .mat: should satisfy L*L' = K
fprintf('\nVerifying L from .mat file:\n');
fprintf('  ||L*L'' - K|| / ||K|| = %.4e\n', norm(L*L' - K, 'fro') / norm(K, 'fro'));
fprintf('  L is lower triangular: %d\n', istril(L));
fprintf('  # zero diagonal entries: %d / %d\n', sum(abs(diag(L)) == 0), m);

%% 1. MATLAB's own Cholesky of K (with fill-reducing ordering)
fprintf('\n--- MATLAB chol(K) ---\n');
tic;
[L1, flag, P_perm] = chol(K, 'lower', 'vector');
t_chol = toc;
if flag ~= 0
    fprintf('  WARNING: chol failed at column %d\n', flag);
else
    fprintf('  chol(K) succeeded in %.2fs\n', t_chol);
    fprintf('  L1: %d x %d, nnz=%d\n', size(L1,1), size(L1,2), nnz(L1));
    fprintf('  L1 is lower triangular: %d\n', istril(L1));
    % Verify: L1*L1' = K(P,P)
    K_perm = K(P_perm, P_perm);
    fprintf('  ||L1*L1'' - K(P,P)|| / ||K|| = %.4e\n', ...
            norm(L1*L1' - K_perm, 'fro') / norm(K, 'fro'));
end

V_dense = full(V);

%% 2. cond(K^{-1} V) via K\V
fprintf('\n--- K^{-1} V (direct solve) ---\n');
tic;
KiV = K \ V_dense;
t_solve = toc;
res_K = norm(K * KiV - V_dense, 'fro') / norm(V_dense, 'fro');
fprintf('  Solve time: %.1fs\n', t_solve);
fprintf('  Residual: ||K*K^{-1}V - V|| / ||V|| = %.4e\n', res_K);

tic;
s_KiV = svd(KiV, 'econ');
fprintf('  cond(K^{-1}V) = %.6e  (%.1fs)\n', s_KiV(1)/s_KiV(end), toc);

%% 3. cond(L^{-1} V) using L from .mat file
fprintf('\n--- L^{-1} V (L from .mat file, non-triangular) ---\n');
tic;
LiV = L \ V_dense;
t_solve = toc;
res_L = norm(L * LiV - V_dense, 'fro') / norm(V_dense, 'fro');
fprintf('  Solve time: %.1fs\n', t_solve);
fprintf('  Residual: ||L*L^{-1}V - V|| / ||V|| = %.4e\n', res_L);

tic;
s_LiV = svd(LiV, 'econ');
fprintf('  cond(L^{-1}V) = %.6e  (%.1fs)\n', s_LiV(1)/s_LiV(end), toc);
fprintf('  sigma_max = %.6e, sigma_min = %.6e\n', s_LiV(1), s_LiV(end));

%% 4. cond(L1^{-1} P V) using MATLAB's own Cholesky factor
fprintf('\n--- L1^{-1} P V (MATLAB chol, proper triangular) ---\n');
PV = V_dense(P_perm, :);   % Apply permutation
tic;
L1iV = L1 \ PV;            % Triangular solve (forward substitution)
t_solve = toc;
res_L1 = norm(L1 * L1iV - PV, 'fro') / norm(PV, 'fro');
fprintf('  Solve time: %.1fs\n', t_solve);
fprintf('  Residual: ||L1*L1^{-1}PV - PV|| / ||PV|| = %.4e\n', res_L1);

tic;
s_L1iV = svd(L1iV, 'econ');
fprintf('  cond(L1^{-1}PV) = %.6e  (%.1fs)\n', s_L1iV(1)/s_L1iV(end), toc);
fprintf('  sigma_max = %.6e, sigma_min = %.6e\n', s_L1iV(1), s_L1iV(end));

%% 5. Condition number of K itself
fprintf('\n--- cond(K) ---\n');
eig_max = eigs(K, 1, 'largestabs');
eig_min = eigs(K, 1, 'smallestabs');
fprintf('  lambda_max(K) = %.6e\n', eig_max);
fprintf('  lambda_min(K) = %.6e\n', eig_min);
fprintf('  cond(K) = %.6e\n', eig_max / eig_min);
fprintf('  cond(L) = sqrt(cond(K)) = %.6e\n', sqrt(eig_max / eig_min));

%% Summary
fprintf('\n========== SUMMARY ==========\n');
fprintf('  cond(K^{-1} V)                            = %.6e  (res=%.1e)\n', ...
        s_KiV(1)/s_KiV(end), res_K);
fprintf('  cond(L^{-1} V)   [.mat L, non-triangular] = %.6e  (res=%.1e)\n', ...
        s_LiV(1)/s_LiV(end), res_L);
fprintf('  cond(L1^{-1} PV) [MATLAB chol, triangular] = %.6e  (res=%.1e)\n', ...
        s_L1iV(1)/s_L1iV(end), res_L1);
fprintf('  C++ benchmark (Eigen TRSM)                = 1.436138e+03\n');
fprintf('  Collaborator                              ~ 8.3e+07\n');

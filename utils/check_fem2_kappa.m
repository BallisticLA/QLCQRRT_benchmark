function check_fem2_kappa(K_file, M_file, V_file, col_spreads)
% check_fem2_kappa - verify the ACTUAL condition number of the FEM2 IR-LSQ
% operator J = L^{-1} K V (L = chol(M), lower), optionally with geometric
% column-scaling of V applied (spread = injected multiplier on kappa).
%
% Uses dense svd VALUES-ONLY (s = svd(full(J))), i.e. the LAPACK gesdd jobz='N'
% path. This resolves sigma_min down to ~u*sigma_max, so it reads kappa up to
% ~1e16. Do NOT use the Gram/eig route (eig of J'J): it squares the conditioning
% and goes numerically singular for kappa(J) > ~1e8, so it cannot see 1e9.
%
% SMALL problem only: J is dense ~5 GB for the 75824x8304 case (needs ~16 GB RAM).
% Medium/large are kron(I_k, .) expansions -> identical kappa, so small suffices.
%
% Inputs:
%   K_file,M_file,V_file : Matrix Market paths (needs mmread on the MATLAB path)
%   col_spreads          : vector of column-scaling spreads to test.
%                          1 = native FEM2 (no scaling). e.g. [1 100] for
%                          target kappa ~ {native, native*100}.
%
% Example:
%   check_fem2_kappa('FEM_Problem_2_corrected_75824x8304_K.mtx', ...
%                    'FEM_Problem_2_corrected_75824x8304_M.mtx', ...
%                    'FEM_Problem_2_corrected_75824x8304_V.mtx', [1 100]);

if nargin < 4 || isempty(col_spreads), col_spreads = [1 100]; end

fprintf('Loading K, M, V ...\n');
K = mmread(K_file);
M = mmread(M_file);
V = mmread(V_file);
[mK, nV] = size(V);
fprintf('  K %dx%d (nnz %d), M %dx%d, V %dx%d\n', ...
        size(K,1), size(K,2), nnz(K), size(M,1), size(M,2), mK, nV);

fprintf('Cholesky M = L L^T ...\n');
L  = chol(M, 'lower');     % sparse lower factor
KV = K * V;                % sparse m_K x n_V

fprintf('\n%-12s %-12s %-12s %-12s  spectrum (head | tail)\n', ...
        'col_spread', 'kappa(J)', 'sigma_max', 'sigma_min');
for sp = col_spreads(:)'
    if sp > 1
        d  = sp .^ ((0:nV-1)' / (nV-1));   % geometric, columns span [1, sp]
        Js = L \ (KV * spdiags(d, 0, nV, nV));
    else
        Js = L \ KV;                        % native FEM2
    end
    s   = svd(full(Js));                    % values only (no U,V) -> cheap path
    kap = s(1) / s(end);
    head = sprintf('%.3g %.3g %.3g', s(1), s(2), s(3));
    tail = sprintf('%.3g %.3g %.3g', s(end-2), s(end-1), s(end));
    fprintf('%-12g %-12.3e %-12.3e %-12.3e  [%s | %s]\n', ...
            sp, kap, s(1), s(end), head, tail);
    clear Js                                % free the ~5 GB before the next spread
end
fprintf('\nNote: kappa(J) here = the true condition number of the operator the\n');
fprintf('benchmark builds. Confirm native (~1e7) and that spread 100 lands ~1e9.\n');
end

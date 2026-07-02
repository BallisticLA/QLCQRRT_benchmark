function [kap, kap_colnorm, smin, smax, colspread] = compute_kappa_variants(K, M, V)
%COMPUTE_KAPPA_VARIANTS  Condition numbers of the FEM operator J = L^{-1} K V.
%
% Single source of truth for the two condition-number metrics of the IR-LSQ
% operator. Dense values-only SVD -> SMALL problems only (J is m x n dense).
%
%   kap         = kappa(J) = sigma_max/sigma_min            (nominal condition number)
%   kap_colnorm = kappa of J with unit-normalized columns   (the COLUMN-EQUILIBRATED
%                 condition number -- what CholeskyQR actually feels; the honest
%                 difficulty, since CholeskyQR factors diagonal column scaling out)
%   smin, smax  = extreme singular values of J
%   colspread   = max/min column norm of J (the scale-away-able part;
%                 roughly  kap ~ kap_colnorm * colspread)
%
% Background: a smooth-contrast FEM field gives kap huge but kap_colnorm ~ tens
% (pure column scaling, easy). A genuinely hard test needs kap_colnorm large.
%
% Usage:
%   [k, kc, smin, smax, spread] = compute_kappa_variants(K, M, V);

    p  = symamd(M);                         % fill-reducing order for the mass Cholesky
    Lp = chol(M(p, p), 'lower');
    J  = full(Lp \ (K(p, :) * V));          % J = L^{-1} K V, row-permuted (kappa invariant)

    s    = svd(J);                          % values only (gesdd 'N')
    smax = s(1);  smin = s(end);  kap = smax / smin;

    cn = vecnorm(J);                        % column norms (1 x n)
    cn(cn == 0) = 1;                        % guard (full-rank J has no zero columns)
    colspread = max(cn) / min(cn);

    sn = svd(J ./ cn);                      % columns normalized to unit norm
    kap_colnorm = sn(1) / sn(end);
end

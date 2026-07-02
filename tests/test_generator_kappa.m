% test_generator_kappa.m
% Verify the hardened FEM generator produces a GENUINELY hard operator:
% kappa^colnorm (the column-equilibrated condition number that CholeskyQR actually
% feels) tracks the target {1e7, 1e9, 1e11}, the operator stays full rank, the
% hardness survives column normalization (kappa^colnorm ~ kappa), and -- crucially --
% the target is hit MESH-INDEPENDENTLY, so a value measured at a small proxy transfers
% to the (SVD-infeasible) near-original size.
%
% Construction (see gen_fem2_hard_variants): uniform-coefficient FEM base (well-cond
% K, M), coarse basis V_hard = (V .* g.^linspace(0,1,n)) * W_butterfly. The geometric
% column scaling g injects the conditioning; the sparse butterfly W spreads it so it
% survives equilibration; g sets kappa^colnorm directly with no h-dependent calibration.
%
% Run:  matlab -batch "run('tests/test_generator_kappa.m')"

here   = fileparts(mfilename('fullpath'));
parent = fileparts(here);
addpath(parent);  addpath(fullfile(parent, 'utils'));

targets = [1e7 1e9 1e11];    % target kappa^colnorm (the honest difficulty axis)
r = 3;  SEED = 1;  NST = 3;
TOL_LO = 0.1;  TOL_HI = 10;  % accept kappa^colnorm within [g/10, 10g] (mesh scatter ~1-7x)

nfail = 0;
fprintf('\n');
for kt = targets
    g = kt;
    % --- two meshes (48x24, 96x48) for the mesh-independence check ---
    [kapc1, kap1, smin1, n1] = build_and_measure(48, 24, r, g, SEED, NST);
    [kapc2, kap2, smin2, n2] = build_and_measure(96, 48, r, g, SEED, NST);

    % Kron x2 on the 96x48 operator must preserve kappa^colnorm exactly
    [K, M, Vh] = build_op(96, 48, r, g, SEED, NST);
    [~, kapc2k] = compute_kappa_variants(kron(speye(2),K), kron(speye(2),M), kron(speye(2),Vh));

    ok_target = kapc2 >= g*TOL_LO && kapc2 <= g*TOL_HI;   % near g
    ok_rank   = smin2 > 0;                                 % full rank
    ok_hard   = kapc2 / kap2 > 0.05;                       % survives equilibration
    ok_mesh   = abs(log10(kapc1/kapc2)) < 1;               % <10x scatter across meshes (no drift)
    ok_kron   = abs(log10(kapc2k/kapc2)) < 0.02;           % preserved under Kron

    pass = ok_target && ok_rank && ok_hard && ok_mesh && ok_kron;
    fprintf(['g=%.0e -> kcol(96)=%.2e kcol(48)=%.2e kappa=%.2e smin=%.1e kron(x2)=%.2e | ' ...
             'tgt=%d rank=%d hard=%d mesh=%d kron=%d  %s\n'], ...
            g, kapc2, kapc1, kap2, smin2, kapc2k, ...
            ok_target, ok_rank, ok_hard, ok_mesh, ok_kron, string(pass));
    if ~pass, nfail = nfail + 1; end
end

fprintf('\n');
if nfail > 0
    error('test_generator_kappa: %d of %d targets FAILED', nfail, numel(targets));
end
disp('test_generator_kappa: ALL PASS');


function [K, M, Vh] = build_op(nxf, nyf, r, g, seed, nst)
    [K, M, V] = gen_fem2_small_local([], 1, nxf, nyf, r, 'smooth', false);  % uniform base
    n  = size(V, 2);
    Vh = apply_butterfly(V .* (g .^ linspace(0,1,n)), nst, seed);
    thr = 1e-14 * max(abs(nonzeros(Vh)));  Vh(abs(Vh) < thr) = 0;
end

function [kapc, kap, smin, n] = build_and_measure(nxf, nyf, r, g, seed, nst)
    [K, M, Vh] = build_op(nxf, nyf, r, g, seed, nst);
    n = size(Vh, 2);
    [kap, kapc, smin] = compute_kappa_variants(K, M, Vh);
end

function V = apply_butterfly(V, nstages, seed)
%APPLY_BUTTERFLY  Right-multiply V by a sparse orthogonal butterfly W (V -> V*W).
%
% W is a product of `nstages` Givens layers. Each layer randomly pairs the n
% columns into n/2 disjoint pairs and rotates each pair by a random angle, so a
% single layer is orthogonal with exactly 2 nnz per column. Over a few layers the
% random pairings give LONG-RANGE column mixing, which is what converts the cheap
% column-scaling ill-conditioning of the FEM coarse basis into genuine,
% equilibration-surviving hardness (kappa^colnorm ~ kappa) while keeping V*W SPARSE
% -- unlike a full dense random orthogonal W, which hardens equally but densifies V.
%
% kappa(J*W) = kappa(J) exactly (W orthogonal); only kappa^colnorm rises. 2-3 stages
% suffice; fill per column converges to a small n-independent constant, so V*W stays
% writable at the original operator size. Reproducible via `seed`.
%
%   V = apply_butterfly(V, nstages, seed)
%
% See also: gen_fem2_hard_variants, compute_kappa_variants.

    if nargin < 2 || isempty(nstages), nstages = 3; end
    if nargin < 3 || isempty(seed),    seed    = 1; end
    n  = size(V, 2);
    rs = RandStream('threefry', 'Seed', seed);

    for L = 1:nstages
        p  = randperm(rs, n);
        np = floor(n / 2);
        a  = p(1:2:2*np);            % first  element of each pair
        b  = p(2:2:2*np);            % second element of each pair
        th = 2*pi * rand(rs, 1, np);
        c  = cos(th);  s = sin(th);
        Va = V(:, a);  Vb = V(:, b);
        % (V*G) columns for a Givens G with G(a,a)=c,G(a,b)=-s,G(b,a)=s,G(b,b)=c:
        V(:, a) =  Va .* c + Vb .* s;
        V(:, b) = -Va .* s + Vb .* c;
        % odd column out (n odd) keeps identity -> untouched
    end
end

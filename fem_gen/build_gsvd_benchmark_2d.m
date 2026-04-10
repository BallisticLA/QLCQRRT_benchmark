function data = build_gsvd_benchmark_2d(varargin)
%BUILD_GSVD_BENCHMARK_2D
% Build a realistic 2D FEM benchmark for generalized/weighted SVD tests.
%
% Outputs a struct with:
%   K, M              : fine-grid reduced (free-DOF) stiffness/mass matrices
%   P                 : tall prolongation (fine_free x coarse_free)
%   Kc_asm, Mc_asm    : directly assembled coarse reduced matrices
%   Kc_gal, Mc_gal    : Galerkin coarse matrices P' K P, P' M P
%   Kmetric           : regularized metric matrix for inverse-inner-product tests
%   freeF, freeC      : free DOF indices on fine/coarse grids (w.r.t. full nodal indexing)
%   nodesF, elemsF    : fine mesh nodes/elements (full)
%   nodesC, elemsC    : coarse mesh nodes/elements (full)
%   Pfull             : full nodal interpolation (fine_full x coarse_full)
%   Bload, Csensor    : optional source/sensor matrices (on free DOFs)
%
% Example:
%   data = build_gsvd_benchmark_2d('HmaxFine',0.05,'HmaxCoarse',0.18,'Plot',true);
%
%   % Largest singular values (model reduction) from snapshots
%   A = data.Kmetric \ data.Bload;              % solution snapshots
%   Rm = chol(data.M,'lower');
%   [Um,Sm,Vm] = svd(Rm' * A, 'econ');
%   Phi = Rm' \ Um;                             % M-orthonormal modes
%   fprintf('||Phi'' M Phi - I|| = %.2e\n', norm(Phi' * data.M * Phi - eye(size(Phi,2))));
%
%   % Smallest singular values (sensitivity / observability proxy)
%   G = data.Csensor * (data.Kmetric \ data.Bload);
%   s = svd(full(G));
%   fprintf('Smallest singular value of transfer = %.3e\n', s(end));
%
% Requires MATLAB PDE Toolbox.

% ----------------------------
% Parse inputs
% ----------------------------
p = inputParser;
p.addParameter('HmaxFine',   0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('HmaxCoarse', 0.18, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Hgrad',      1.35, @(x)isnumeric(x)&&isscalar(x)&&x>1);
p.addParameter('Plot',       true, @(x)islogical(x)||ismember(x,[0 1]));
p.addParameter('AlphaReg',   1e-8, @(x)isnumeric(x)&&isscalar(x)&&x>=0); % Kmetric = K + alpha*scale*M
p.parse(varargin{:});
opts = p.Results;

% ----------------------------
% Geometry: complex plate with holes + slit-like cutouts + re-entrant corner
% Domain roughly [0,4] x [0,2]
% ----------------------------
[dl, geomInfo] = makeComplexGeometryCSG();

% Fine model
modelF = createpde(1);
geometryFromEdges(modelF, dl);

% Coefficient field (high contrast / nontrivial)
specifyCoefficients(modelF, ...
    "m", 0, ...
    "d", 1, ...   % ensures a nonzero mass matrix M in assembly
    "c", @cCoeffIllConditioned, ...
    "a", 0, ...
    "f", 0);
%"c", @cCoeffHeterogeneous, ...
% Zero Dirichlet on all boundaries => SPD after elimination
applyBoundaryCondition(modelF, "dirichlet", ...
    "Edge", 1:modelF.Geometry.NumEdges, "u", 0);

% Linear triangles make nodal interpolation P straightforward
generateMesh(modelF, ...
    "Hmax", opts.HmaxFine, ...
    "Hgrad", opts.Hgrad, ...
    "GeometricOrder", "linear");

meshF = modelF.Mesh;

% Assemble full matrices (we will eliminate boundary DOFs ourselves)
FEMF = assembleFEMatrices(modelF, "none");
KfullF = FEMF.K;
MfullF = FEMF.M;

% Coarse model (same geometry / coefficients)
modelC = createpde(1);
geometryFromEdges(modelC, dl);
specifyCoefficients(modelC, ...
    "m", 0, ...
    "d", 1, ...
    "c", @cCoeffIllConditioned, ...
    "a", 0, ...
    "f", 0);
applyBoundaryCondition(modelC, "dirichlet", ...
    "Edge", 1:modelC.Geometry.NumEdges, "u", 0);

generateMesh(modelC, ...
    "Hmax", opts.HmaxCoarse, ...
    "Hgrad", opts.Hgrad, ...
    "GeometricOrder", "linear");

meshC = modelC.Mesh;

FEMC = assembleFEMatrices(modelC, "none");
KfullC = FEMC.K;
MfullC = FEMC.M;

% ----------------------------
% Identify boundary/free DOFs from triangulation topology
% (robust and independent of edge-number queries)
% ----------------------------
TF = double(meshF.Elements(1:3,:));
TC = double(meshC.Elements(1:3,:));

nF = size(meshF.Nodes,2);
nC = size(meshC.Nodes,2);

bndF = boundaryNodesFromTriangles(TF);
bndC = boundaryNodesFromTriangles(TC);

freeF = setdiff((1:nF)', bndF(:));
freeC = setdiff((1:nC)', bndC(:));

% Reduced fine/coarse matrices
K = KfullF(freeF, freeF);
M = MfullF(freeF, freeF);

Kc_asm = KfullC(freeC, freeC);
Mc_asm = MfullC(freeC, freeC);

% ----------------------------
% Build prolongation / coarse basis matrix P via nodal interpolation
% Pfull: maps coarse nodal values -> fine nodal values
% ----------------------------
Pfull = nodalProlongationFromMesh(meshC, meshF);   % (nF x nC)

% Restrict to free DOFs on both grids
P = Pfull(freeF, freeC);

% Galerkin coarse matrices induced by P (useful comparison)
Kc_gal = P' * K * P;
Mc_gal = P' * M * P;

% Regularized metric for inverse-inner-product tests (safety)
% (if user wants K^{-1}-type inner product but with robust SPD numerics)
scaleKM = trace(K) / max(trace(M), eps);
Kmetric = K + opts.AlphaReg * scaleKM * M;

% ----------------------------
% Optional source/sensor matrices for SVD tests
% Bload  : localized loads (columns) on free DOFs
% Csensor: point-sampling sensors (rows) on free DOFs
% ----------------------------
Xf = meshF.Nodes(:, freeF)';   % (#freeF x 2)

srcCenters = [ ...
    0.45 0.35;
    0.75 1.65;
    1.95 0.28;
    2.55 1.70;
    3.35 0.45;
    3.55 1.35;
    2.05 1.05];   % near crack-ish / reentrant regions too

snsCenters = [ ...
    1.70 1.00;    % near horizontal slit tip
    2.55 1.25;    % near vertical slit neighborhood
    3.05 0.95;    % near reentrant corner ligament
    1.20 0.55;
    3.55 0.55;
    2.20 1.75];

srcIdx = nearestFreeNodes(Xf, srcCenters);
snsIdx = nearestFreeNodes(Xf, snsCenters);

% Consistent localized loads using mass columns (smeared nodal basis loads)
Bload = M(:, srcIdx);

% Point sensors (can replace with weighted averaging if preferred)
Csensor = sparse(1:numel(snsIdx), snsIdx, 1, numel(snsIdx), size(M,1));

% ----------------------------
% Optional plots
% ----------------------------
if opts.Plot
    figure('Name','Geometry (edge labels)'); 
    pdegplot(modelF, "EdgeLabels", "on"); axis equal tight
    title('Complex 2D geometry (holes + slit-like cutouts + reentrant corner)');

    figure('Name','Fine mesh');
    pdemesh(modelF); axis equal tight
    title(sprintf('Fine mesh (Hmax = %.3g), n = %d nodes', opts.HmaxFine, nF));

    figure('Name','Coarse mesh');
    pdemesh(modelC); axis equal tight
    title(sprintf('Coarse mesh (Hmax = %.3g), n = %d nodes', opts.HmaxCoarse, nC));

    % Plot coefficient field sampled on fine nodes (for inspection only)
    xy = meshF.Nodes;
    cvals = log(cCoeffIllConditioned(struct('x',xy(1,:),'y',xy(2,:)), []));
    figure('Name','Illconditioned coefficient c(x,y)');
    pdeplot(modelF, "XYData", cvals, "Mesh", "off"); axis equal tight
    colorbar; title('Conductivity / stiffness coefficient c(x,y)');
end

% ----------------------------
% Package outputs
% ----------------------------
data = struct();
data.geometry = geomInfo;

data.K = K;
data.M = M;
data.Kmetric = Kmetric;
data.P = P;

data.Kc_asm = Kc_asm;
data.Mc_asm = Mc_asm;
data.Kc_gal = Kc_gal;
data.Mc_gal = Mc_gal;

data.KfullF = KfullF;
data.MfullF = MfullF;
data.KfullC = KfullC;
data.MfullC = MfullC;

data.Pfull = Pfull;
data.freeF = freeF;
data.freeC = freeC;
data.bndF = bndF;
data.bndC = bndC;

data.nodesF = meshF.Nodes;         % 2 x nF
data.elemsF = meshF.Elements(1:3,:); % 3 x nElemF (linear triangles)
data.nodesC = meshC.Nodes;
data.elemsC = meshC.Elements(1:3,:);

data.Bload = Bload;
data.Csensor = Csensor;
data.srcIdx = srcIdx;
data.snsIdx = snsIdx;

% quick sanity checks
data.checks = struct();
data.checks.P_size = size(P);
data.checks.nnzP = nnz(P);
data.checks.symK = norm(K-K','fro')/max(norm(K,'fro'),eps);
data.checks.symM = norm(M-M','fro')/max(norm(M,'fro'),eps);
data.checks.rowSumsPfull = [min(sum(Pfull,2)), max(sum(Pfull,2))]; % ~1 except numerical misses fixed by fallback

fprintf('\n=== Benchmark summary ===\n');
fprintf('Fine nodes (full/free):   %d / %d\n', nF, numel(freeF));
fprintf('Coarse nodes (full/free): %d / %d\n', nC, numel(freeC));
fprintf('P size (free):            %d x %d\n', size(P,1), size(P,2));
fprintf('nnz(P):                   %d\n', nnz(P));
fprintf('rel symmetry error K:     %.2e\n', data.checks.symK);
fprintf('rel symmetry error M:     %.2e\n', data.checks.symM);

end

% ========================================================================
% Local functions
% ========================================================================

function [dl, info] = makeComplexGeometryCSG()
% Construct a complex 2D domain using CSG primitives.
%
% Base plate:      R1
% Removed slits:   R2, R3 (thin rectangles -> "crack-like" cutouts)
% Reentrant notch: R4
% Holes:           C1, C2, C3

% Base rectangle R1 = [0,4] x [0,2]
R1 = [3 4 0 4 4 0 0 0 2 2]';

% Thin horizontal slit-like void (internal notch)
R2 = [3 4 0.30 1.70 1.70 0.30 0.97 0.97 1.003 1.003]';

% Thin vertical slit-like void
R3 = [3 4 2.59 2.60 2.60 2.59 1.015 1.015 1.85 1.85]';

% Reentrant corner cutout (upper-right bite)
R4 = [3 4 2.95 4.00 4.00 2.95 1.05 1.05 2.00 2.00]';

% Circular holes
C1 = [1 1.10 0.55 0.016 0 0 0 0 0 0]';
C2 = [1 2.15 0.55 0.12 0 0 0 0 0 0]';
C3 = [1 3.35 0.55 0.18 0 0 0 0 0 0]';

gd = [R1 R2 R3 R4 C1 C2 C3];

% Names: each column corresponds to one geometry object
ns = char('R1','R2','R3','R4','C1','C2','C3')';
sf = 'R1 - R2 - R3 - R4 - C1 - C2 - C3';

[dl, bt] = decsg(gd, sf, ns); %#ok<ASGLU>

info = struct();
info.setFormula = sf;
info.base = 'Plate [0,4]x[0,2] with holes, slit-like voids, and reentrant corner';
end

function c = cCoeffHeterogeneous(location, ~)
% High-contrast scalar coefficient c(x,y) for -div(c grad u).
% Designed to create nontrivial spectra and localized modes.

x = location.x;
y = location.y;
n = numel(x);
c = ones(1,n);

% Soft spots near selected regions (can create small singular values / sensitivity)
soft1 = ((x-1.65).^2 + (y-1.00).^2) < 0.10^2;   % near slit tip
soft2 = ((x-2.56).^2 + (y-1.18).^2) < 0.08^2;   % near vertical slit base
soft3 = ((x-3.00).^2 + (y-1.05).^2) < 0.12^2;   % reentrant corner neighborhood
c(soft1) = 1.e-8;
c(soft2) = 1.e-7;
c(soft3) = 1.e-6;

% Very stiff channels / inclusions (contrast broadens spectrum)
chan1 = abs(y - (0.35 + 0.10*sin(2*pi*x/4))) < 0.035 & x>0.15 & x<3.8;
chan2 = abs(x - 2.35) < 0.045 & y>0.15 & y<1.8;
chan3 = abs(y - 1.78) < 0.03  & x>1.4  & x<2.8;
c(chan1) = 1.e8;
c(chan2) = 1.e7;
c(chan3) = 1.e6;

% Moderately soft annuli around holes (more local structure)
r1 = hypot(x-1.10, y-0.55);
r2 = hypot(x-2.15, y-0.55);
r3 = hypot(x-3.35, y-0.55);
ring = (abs(r1-0.22)<0.03) | (abs(r2-0.18)<0.025) | (abs(r3-0.24)<0.03);
c(ring) = min(c(ring), 1.e-8);

% Clamp (safety)
c = max(c, 1.e-10);
%c = c*0+1;
end

function bnd = boundaryNodesFromTriangles(T)
% T: 3 x nElem integer node indices (triangles)
T = double(T);
E = [T([1 2],:)'; T([2 3],:)'; T([3 1],:)'];
E = sort(E, 2);

[Eu, ~, ic] = unique(E, 'rows');
counts = accumarray(ic, 1);
bndEdges = Eu(counts == 1, :);
bnd = unique(bndEdges(:));
end

function P = nodalProlongationFromMesh(meshC, meshF)
% Build full nodal interpolation matrix P (fine <- coarse) using barycentric coords.
%
% For each fine node x_f:
%   - find containing coarse triangle
%   - write x_f as barycentric combo of the triangle's coarse vertices
% Then P(i,:) applied to coarse nodal values yields interpolated fine value at node i.

XC = meshC.Nodes';                  % nC x 2
TC = double(meshC.Elements(1:3,:))';% nElemC x 3
XF = meshF.Nodes';                  % nF x 2

TR = triangulation(TC, XC);

[tid, bc] = pointLocation(TR, XF);  % tid: nF x 1, bc: nF x 3
inside = ~isnan(tid);

nF = size(XF,1);
nC = size(XC,1);

% Main barycentric interpolation entries
ii = find(inside);
triVerts = TC(tid(inside), :);      % nInside x 3
bcInside = bc(inside, :);           % nInside x 3

rows = repelem(ii, 3, 1);
cols = reshape(triVerts.', [], 1);
vals = reshape(bcInside.', [], 1);

P = sparse(rows, cols, vals, nF, nC);

% Fallback for missed points (usually boundary roundoff)
miss = find(~inside);
if ~isempty(miss)
    warning('nodalProlongationFromMesh:pointMiss', ...
        '%d fine nodes were not located in coarse mesh triangles; using nearest coarse node fallback.', ...
        numel(miss));
    for k = 1:numel(miss)
        j = nearestPointIndex(XC, XF(miss(k),:));
        P(miss(k), j) = 1;
    end
end
end

function idx = nearestFreeNodes(X, centers)
% X       : n x 2 free-node coordinates
% centers : q x 2 target points
q = size(centers,1);
idx = zeros(q,1);
for j = 1:q
    idx(j) = nearestPointIndex(X, centers(j,:));
end
end

function j = nearestPointIndex(X, x0)
d2 = (X(:,1)-x0(1)).^2 + (X(:,2)-x0(2)).^2;
[~, j] = min(d2);
end

function c = cCoeffIllConditioned(location, state) %#ok<INUSD>
% Scalar conductivity c(x,y) for PDE Toolbox:
% -div(c grad u) = ...
% Returns a 1 x N row vector, strictly positive.
%
% Designed to make K ill-conditioned (with Dirichlet BCs still SPD).

x = location.x;
y = location.y;
N = numel(x);

% Baseline
c = ones(1, N)*10^(-2);

% ------------------------------------------------------------
% A) Ultra-soft thin ligament (dominant source of ill-conditioning)
%    Almost disconnects left/right regions while keeping SPD.
% ------------------------------------------------------------
lig = (x > 0.995 & x < 1.005 & y > 0.02 & y < 1.98);
c(lig) = min(c(lig), 1e-8);

% Even softer "neck" in the middle of the ligament
neck = (x > 0.992 & x < 1.008 & abs(y - 1.0) < 0.025);
c(neck) = min(c(neck), 1e-8);

% ------------------------------------------------------------
% B) Soft pockets (localized weak regions, e.g. near slit/crack tips)
% ------------------------------------------------------------
soft1 = ((((x - 1.70).^2 + (y - 1.00).^2) < (0.12)^2) & y>1);
soft2 = ((x - 2.58).^2 + (y - 1.20).^2) < (0.09)^2;
soft3 = ((x - 3.02).^2 + (y - 1.08).^2) < (0.13)^2;
soft4 = ((((x - 1.70).^2 + (y - 1.00).^2) < (0.12)^2) & y<=1);


c(soft1) = min(c(soft1), 1e-4);
c(soft2) = min(c(soft2), 1e-5);
c(soft3) = min(c(soft3), 1e-6);
c(soft4) = max(c(soft4), 1e4);
% ------------------------------------------------------------
% C) Very stiff channels (increases contrast and spectral spread)
% ------------------------------------------------------------
chan1 = abs(y - (0.35 + 0.10*sin(2*pi*x/4))) < 0.03 & x > 0.15 & x < 3.85;
chan2 = abs(x - 2.35) < 0.04 & y > 0.12 & y < 1.88;
chan3 = abs(y - 1.78) < 0.025 & x > 1.35 & x < 2.90;

c(chan1) = max(c(chan1), 1e5);
c(chan2) = max(c(chan2), 1e6);
c(chan3) = max(c(chan3), 1e7);

% ------------------------------------------------------------
% D) Moderately soft annuli around holes (if your geometry has these holes)
% ------------------------------------------------------------
r1 = hypot(x - 1.10, y - 0.55);
r2 = hypot(x - 2.15, y - 0.55);
r3 = hypot(x - 3.35, y - 0.55);

ring = (abs(r1 - 0.22) < 0.25) | (abs(r2 - 0.18) < 0.02) | (abs(r3 - 0.24) < 0.025);
c(ring) = min(c(ring), 1e-6);

ring = abs(r1 - 0.22) < 0.1;
c(ring) = 1e6;

% Strict positivity safeguard
c = max(c, 1e-8);
end
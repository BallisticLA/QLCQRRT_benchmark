%% FEM matrices generator

opts = struct('HmaxFine',   0.01, ...
     'HmaxCoarse', 0.03, ...
     'AlphaReg',   1.e-14, ...
     'Plot',       false);

data = build_gsvd_benchmark_2d(opts);

Problem.K = data.Kmetric;
% L = chol(Problem.K);
% nnz(L)/size(L,1)

[R,flag,P] = chol(Problem.K,'vector');
nnz(R)/size(R,1)
norm(R'*R-Problem.K(P,P),'fro')/norm(Problem.K,'fro')
Pc = [];
Pc(P) = 1:length(P);
L = R(:,Pc)';
norm(L*L'-Problem.K(:,:),'fro')/norm(Problem.K,'fro')

Problem.V = data.P;
Problem.L = L;

nnz(Problem.K\Problem.V)/size(Problem.V,1)

d = diag(full(vecnorm(Problem.L\Problem.V)));
W = Problem.L\full(Problem.V);
G = W'*W;
condest(G)
W = W./diag(d)';
d_W = vecnorm(W);
G = W'*W;
condest(G)
Problem.options = opts;
save('FEM_Problem.mat', 'Problem');
% 




%% Cholesky QR test
load('FEM_Problem.mat')
d = diag(full(vecnorm(Problem.L\Problem.V)));

W = Problem.L\full(Problem.V);
G = W'*W;
condest(G)
R = chol(G);
Q = W/R;
norm(Q'*Q-eye(size(Q,2),size(Q,2)))


W = W./diag(d)';
d_W = vecnorm(W);
G = W'*W;
condest(G)






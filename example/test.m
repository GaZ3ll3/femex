function  fem = test(prec, min_area)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('~/Documents/github/femex'));

fem = FEM(prec, min_area);
num_elems = size(fem.Promoted.elems, 2);
num_fqnodes = size(fem.Facet.Qnodes, 2);
num_edges = size(fem.Promoted.edges, 2);
num_eqnodes = size(fem.Edge.Qnodes, 2);
N = size(fem.Promoted.nodes, 2);
numofnodes = fem.Num_nodes;
k = 1;
kk = 1;
for i = 1: size(fem.Promoted.edges, 2)
    if abs(fem.Promoted.nodes(1, fem.Promoted.edges(1, i)) - 1) < 1e-8 && abs(fem.Promoted.nodes(1, fem.Promoted.edges(2, i)) - 1) < 1e-8
        bc_1.bc(:, k) = fem.Promoted.edges(:, i);
        k = k + 1;
    else
        bc_2.bc(:, kk) = fem.Promoted.edges(:, i);
        kk = kk + 1;
    end
end
bc_1.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc_1.bc);
bc_2.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bc_2.bc);
M = fem.assema(ones(num_elems, num_fqnodes));
S = fem.assems(ones(num_elems, num_fqnodes));
Q = fem.assemlbc(ones(num_edges, num_eqnodes), fem.Promoted.edges);
% boundary condition
f = @(x,y)(6*x);
bc_1.robin = @(x,y)(x.*x.*x + 3);
bc_2.robin = @(x,y)(x.*x.*x);
Load = f(fem.Qnodes(1,:), fem.Qnodes(2,:));
bc_1.Load = bc_1.robin(bc_1.qnodes1D(1,:), bc_1.qnodes1D(2, :));
bc_2.Load = bc_2.robin(bc_2.qnodes1D(1,:), bc_2.qnodes1D(2, :));
LoadVector = fem.asseml(Load);
bc_1.LoadVector = fem.assemrbc(bc_1.Load, bc_1.bc);
bc_2.LoadVector = fem.assemrbc(bc_2.Load, bc_2.bc);
LoadVector = LoadVector - bc_1.LoadVector - bc_2.LoadVector;
Dirichlet = [];
boundary = Boundary(Dirichlet);
dofs = boundary.dofs(N, Dirichlet);
fem.Solution = zeros(N, 1);
%true
v = fem.Promoted.nodes(1,:).^3;
R = S+Q;
% apply Dirichlet
LoadVector = LoadVector + R*fem.Solution;


solver = Solver('ilu');
tic;
fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();
disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

% solver = Solver('umfpack');
solver = Solver('gmres');

tic;
fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));
solver = Solver('umfpack');
% solver = Solver('gmres');

tic;
fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

% trimesh(FEM.TriMesh', FEM.Promoted.nodes(1,1:numofnodes), ...
%     FEM.Promoted.nodes(2, 1:numofnodes), FEM.Solution(1:numofnodes));
disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));






end


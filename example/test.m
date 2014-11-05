function  FEM = test(prec, min_area)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

FEM = fem(prec, min_area);
num_elems = size(FEM.Promoted.elems, 2);
num_fqnodes = size(FEM.Facet.Qnodes, 2);
num_edges = size(FEM.Promoted.edges, 2);
num_eqnodes = size(FEM.Edge.Qnodes, 2);
N = size(FEM.Promoted.nodes, 2);
numofnodes = FEM.Num_nodes;
k = 1;
kk = 1;
for i = 1: size(FEM.Promoted.edges, 2)
    if abs(FEM.Promoted.nodes(1, FEM.Promoted.edges(1, i)) - 1) < 1e-8 && abs(FEM.Promoted.nodes(1, FEM.Promoted.edges(2, i)) - 1) < 1e-8
        bc_1.bc(:, k) = FEM.Promoted.edges(:, i);
        k = k + 1;
    else
        bc_2.bc(:, kk) = FEM.Promoted.edges(:, i);
        kk = kk + 1;
    end
end
bc_1.qnodes1D = FEM.Assembler.qnodes1D(FEM.Promoted.nodes, FEM.Edge.Qnodes, bc_1.bc);
bc_2.qnodes1D = FEM.Assembler.qnodes1D(FEM.Promoted.nodes, FEM.Edge.Qnodes ,bc_2.bc);
M = FEM.assema(ones(num_elems, num_fqnodes));
S = FEM.assems(ones(num_elems, num_fqnodes));
Q = FEM.assemlbc(ones(num_edges, num_eqnodes), FEM.Promoted.edges);
% boundary condition
f = @(x,y)(6*x);
bc_1.robin = @(x,y)(x.*x.*x + 3);
bc_2.robin = @(x,y)(x.*x.*x);
Load = f(FEM.Qnodes(1,:), FEM.Qnodes(2,:));
bc_1.Load = bc_1.robin(bc_1.qnodes1D(1,:), bc_1.qnodes1D(2, :));
bc_2.Load = bc_2.robin(bc_2.qnodes1D(1,:), bc_2.qnodes1D(2, :));
LoadVector = FEM.asseml(Load);
bc_1.LoadVector = FEM.assemrbc(bc_1.Load, bc_1.bc);
bc_2.LoadVector = FEM.assemrbc(bc_2.Load, bc_2.bc);
LoadVector = LoadVector - bc_1.LoadVector - bc_2.LoadVector;
Dirichlet = [];
boundary = Boundary(Dirichlet);
dofs = boundary.dofs(N, Dirichlet);
FEM.Solution = zeros(N, 1);
%true
v = FEM.Promoted.nodes(1,:).^3;
R = S+Q;
% apply Dirichlet
LoadVector = LoadVector + R*FEM.Solution;
FEM.Solution(dofs) = -R(dofs, dofs)\LoadVector(dofs);
% trimesh(FEM.TriMesh', FEM.Promoted.nodes(1,1:numofnodes), ...
%     FEM.Promoted.nodes(2, 1:numofnodes), FEM.Solution(1:numofnodes));
disp(norm(FEM.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));






end


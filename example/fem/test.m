function  fem = test(prec, min_area)
% TEST solves Neumann/Dirichlet/Robin boundary condition PDE with h-p
% finite element method.

addpath(genpath('~/Documents/github/femex'));

fem = FEM([0 0 1 0 1 1 0 1]', prec, min_area);
N = size(fem.Promoted.nodes, 2);
numofnodes = fem.Num_nodes;

Dirichlet = [];
% boundary = Boundary(Dirichlet);
boundary = Boundary();
boundary.setDirichlet(Dirichlet);
dofs = boundary.dofs(N);

boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');

tic;
bcs = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);
toc;
bc1.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bcs{1});
bc2.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bcs{2});
bc3.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bcs{3});
bc4.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bcs{4});

M = fem.assema(1);
S = fem.assems(1);
Q = fem.assemlbc(1, fem.Promoted.edges);

f = @(x)(6*x(1,:));
bc1.robin = @(x)(x(1,:).^3 + 3);
bc2.robin = @(x)(x(1,:).^3);
bc3.robin = @(x)(x(1,:).^3);
bc4.robin = @(x)(x(1,:).^3);

Load = f(fem.Qnodes);
bc1.Load = bc1.robin(bc1.qnodes1D);
bc2.Load = bc2.robin(bc2.qnodes1D);
bc3.Load = bc3.robin(bc3.qnodes1D);
bc4.Load = bc4.robin(bc4.qnodes1D);

LoadVector = fem.asseml(Load);

bc1.LoadVector = fem.assemrbc(bc1.Load, bcs{1});
bc2.LoadVector = fem.assemrbc(bc2.Load, bcs{2});
bc3.LoadVector = fem.assemrbc(bc3.Load, bcs{3});
bc4.LoadVector = fem.assemrbc(bc4.Load, bcs{4});

LoadVector = LoadVector - bc1.LoadVector - bc2.LoadVector - bc3.LoadVector - bc4.LoadVector;

% preallocate memory
fem.Solution = zeros(N, 1);
%true
v = fem.Promoted.nodes(1,:).^3;
R = S+Q;
% apply Dirichlet
LoadVector = LoadVector + R*fem.Solution;

% 
% solver = Solver('ilu');
% tic;
% fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
% toc;
% solver.delete();
% disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

solver = Solver('umfpack');
tic;
fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));


solver = Solver('agmg');
tic;
fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

trimesh(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), fem.Solution(1:numofnodes));
end


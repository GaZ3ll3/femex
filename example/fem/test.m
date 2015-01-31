function  [fem, neumann] = test(prec, min_area)
% TEST solves Neumann/Dirichlet/Robin boundary condition PDE with h-p
% finite element method.

addpath(genpath('~/Documents/github/femex'));

fem = FEM([0 0 1 0 1 1 0 1]', prec, min_area);

N = size(fem.Promoted.nodes, 2);
numofnodes = fem.Num_nodes;

% boundary = Boundary(Dirichlet);
boundary = Boundary();


boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');

[bc1.bc, bc2.bc, bc3.bc, bc4.bc] = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);

boundary.setDirichlet(bc1.bc);
boundary.setDirichlet(bc2.bc);

[dofs, ndofs] = boundary.dofs(N);

% bc1.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc1.bc);
% bc2.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bc2.bc);
bc3.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc3.bc);
bc4.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bc4.bc);

% M = fem.assema(1);
S = fem.assems(1);
Q = fem.assemlbc(1, bc3.bc) + fem.assemlbc(1, bc4.bc);

f = @(x)(6*x(1,:));
% bc1.robin = @(x)(x(1,:).^3 + 3);
% bc2.robin = @(x)(x(1,:).^3);
bc3.robin = @(x)(x(1,:).^3);
bc4.robin = @(x)(x(1,:).^3);

Load = f(fem.Qnodes);
% bc1.Load = bc1.robin(bc1.qnodes1D);
% bc2.Load = bc2.robin(bc2.qnodes1D);
bc3.Load = bc3.robin(bc3.qnodes1D);
bc4.Load = bc4.robin(bc4.qnodes1D);

LoadVector = fem.asseml(Load);

% bc1.LoadVector = fem.assemrbc(bc1.Load, bc1.bc);
% bc2.LoadVector = fem.assemrbc(bc2.Load, bc2.bc);
bc3.LoadVector = fem.assemrbc(bc3.Load, bc3.bc);
bc4.LoadVector = fem.assemrbc(bc4.Load, bc4.bc);

LoadVector = LoadVector - ...
 bc3.LoadVector - bc4.LoadVector; % - bc1.LoadVector - bc2.LoadVector 

% preallocate memory
fem.Solution = zeros(N, 1);
%true
v = fem.Promoted.nodes(1,:).^3;
fem.Solution(ndofs) = v(ndofs);
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


A = R(dofs, dofs);
b = LoadVector(dofs);

solver = Solver('umfpack');
tic;
fem.Solution(dofs) = - solver.solve(A, b);
toc;
[DX, DY] = solver.reference(fem.Ref_points);
[GX, GY] = solver.grad(fem.Solution, fem.Promoted.nodes, fem.Promoted.elems, DX, DY);
solver.delete();

disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));


%%% use MA57
control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;

tic;
[str,~,~] = ma57_factor(A,control);
[fem.Solution(dofs),~,~] = ma57_solve(A,-b,str);
toc;
disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

%%%


%%% use Mi20

% tic;
% hsl_mi20_startup;
% hsl_control = hsl_mi20_control;
% inform = hsl_mi20_setup(A,hsl_control);
% fem.Solution(dofs) = pcg(A, -b, 1e-12, 1000,'hsl_mi20_precondition');
% hsl_mi20_finalize;
% toc;
% disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));


%%%


%%% pardiso

% info = pardisoinit(-2,0);
% 
% verbose = false;
% 
% tic;
% info = pardisoreorder(tril(A),info,verbose);
% info = pardisofactor(tril(A),info,verbose);
% % Compute the solutions X using the symbolic factorization.
% 
% [fem.Solution(dofs), info] = pardisosolve(tril(A),-b,info,verbose);
% toc;
% fprintf('PARDISO performed %d iterative refinement steps.\n',info.iparm(7));
% 
% 
% disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));
% 
% pardisofree(info);
% clear info


%%%





solver = Solver('agmg');
tic;
fem.Solution(dofs) = - solver.solve(A, b);
toc;


disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

% needs to allocate memory for u and v.

u = zeros(1, N);
v = zeros(1, N);


for i = 1:size(fem.Promoted.elems, 2)
u(fem.Promoted.elems(:, i)) = GX(:, i);
v(fem.Promoted.elems(:, i)) = GY(:, i);
end

neumann = solver.Neumann(GX, GY, fem.Promoted.nodes, fem.Promoted.elems, fem.Promoted.edges, fem.Promoted.indices);

solver.delete();
quiver3(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), fem.Solution', u, v, -ones(size(fem.Solution, 2),size(fem.Solution, 1)));

hold on
trimesh(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), fem.Solution(1:numofnodes));
end


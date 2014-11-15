function [fem] = forward_solver(fem)

num_elems = size(fem.Promoted.elems, 2);
num_fqnodes = size(fem.Facet.Qnodes, 2);
num_edges = size(fem.Promoted.edges, 2);
num_eqnodes = size(fem.Edge.Qnodes, 2);
N = size(fem.Promoted.nodes, 2);
numofnodes = fem.Num_nodes;


% Neumann boundary condition
k1 = 0; k2 = 0;
k3 = 0; k4 = 0;

% dynamically increase
for i = 1: size(fem.Promoted.edges, 2)
    % right side
    if abs(fem.Promoted.nodes(1, fem.Promoted.edges(1, i)) - 1) < 1e-8 && abs(fem.Promoted.nodes(1, fem.Promoted.edges(2, i)) - 1) < 1e-8
        k1 = k1 + 1;
        bc_1.bc(:, k1) = fem.Promoted.edges(:, i);
    % top side
    elseif abs(fem.Promoted.nodes(2, fem.Promoted.edges(1, i)) - 1) < 1e-8 && abs(fem.Promoted.nodes(2, fem.Promoted.edges(2, i)) - 1) < 1e-8
        k2 = k2 + 1;
        bc_2.bc(:, k2) = fem.Promoted.edges(:, i);  
    % left side
    elseif abs(fem.Promoted.nodes(1, fem.Promoted.edges(1, i)) ) < 1e-8 && abs(fem.Promoted.nodes(1, fem.Promoted.edges(2, i)) ) < 1e-8
        k3 = k3 + 1;
        bc_3.bc(:, k3) = fem.Promoted.edges(:, i);
    % bottom side
    elseif  abs(fem.Promoted.nodes(2, fem.Promoted.edges(1, i)) ) < 1e-8 && abs(fem.Promoted.nodes(2, fem.Promoted.edges(2, i)) ) < 1e-8
        k4 = k4 + 1;
        bc_4.bc(:, k4) = fem.Promoted.edges(:, i);        
    end
    
end

bc_1.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc_1.bc);
bc_2.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes ,bc_2.bc);
bc_3.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc_3.bc);
bc_4.qnodes1D = fem.Assembler.qnodes1D(fem.Promoted.nodes, fem.Edge.Qnodes, bc_4.bc);

% use function to project onto qnodes
M = fem.assema(ones(num_elems, num_fqnodes));
S = fem.assems(ones(num_elems, num_fqnodes));
% Q = fem.assemlbc(ones(num_edges, num_eqnodes), fem.Promoted.edges);

% boundary condition
f = @(x,y)(6*x + 0.1*x.*x.*x);

% bc_1.robin = @(x,y)(x.*x.*x + 3);
% bc_2.robin = @(x,y)(x.*x.*x);


bc_1.robin = @(x, y)(3*ones(size(x, 1), size(x,2)));
bc_2.robin = @(x, y)(zeros(size(x, 1), size(x,2)));
bc_3.robin = @(x, y)(zeros(size(x, 1), size(x,2)));
bc_4.robin = @(x, y)(zeros(size(x, 1), size(x,2)));


Load = f(fem.Qnodes(1,:), fem.Qnodes(2,:));
bc_1.Load = bc_1.robin(bc_1.qnodes1D(1,:), bc_1.qnodes1D(2, :));
bc_2.Load = bc_2.robin(bc_2.qnodes1D(1,:), bc_2.qnodes1D(2, :));
bc_3.Load = bc_3.robin(bc_3.qnodes1D(1,:), bc_3.qnodes1D(2, :));
bc_4.Load = bc_4.robin(bc_4.qnodes1D(1,:), bc_4.qnodes1D(2, :));


LoadVector = fem.asseml(Load);

bc_1.LoadVector = fem.assemrbc(bc_1.Load, bc_1.bc);
bc_2.LoadVector = fem.assemrbc(bc_2.Load, bc_2.bc);
bc_3.LoadVector = fem.assemrbc(bc_3.Load, bc_3.bc);
bc_4.LoadVector = fem.assemrbc(bc_4.Load, bc_4.bc);

LoadVector = LoadVector - bc_1.LoadVector - bc_2.LoadVector - bc_3.LoadVector - bc_4.LoadVector;

Dirichlet = [];
boundary = Boundary(Dirichlet);
dofs = boundary.dofs(N, Dirichlet);


fem.Solution = zeros(N, 1);
%true
v = fem.Promoted.nodes(1,:).^3;
R = 0.1*M - S ;
% apply Dirichlet
LoadVector = LoadVector + R*fem.Solution;

% 
% solver = Solver('ilu');
% tic;
% fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
% toc;
% solver.delete();
% disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));


% solver = Solver('umfpack');
% solver = Solver('gmres');
% 
% tic;
% fem.Solution(dofs) = - solver.solve(R(dofs, dofs), LoadVector(dofs));
% toc;
% solver.delete();
% 
% disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));
solver = Solver('umfpack');
% solver = Solver('gmres');

tic;
fem.Solution(dofs) = solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

% trimesh(FEM.TriMesh', FEM.Promoted.nodes(1,1:numofnodes), ...
%     FEM.Promoted.nodes(2, 1:numofnodes), FEM.Solution(1:numofnodes));
disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

solver = Solver('agmg');
tic;
fem.Solution(dofs) = solver.solve(R(dofs, dofs), LoadVector(dofs));
toc;
solver.delete();

disp(norm(fem.Solution(1:numofnodes) - v(1:numofnodes)')/sqrt(double(numofnodes)));

trimesh(fem.TriMesh', fem.Promoted.nodes(1,1:numofnodes), ...
    fem.Promoted.nodes(2, 1:numofnodes), fem.Solution(1:numofnodes));


end


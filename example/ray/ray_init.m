function [fem, dom] = ray_init()

%addpath(genpath('../../'));

global bc1 bc2 bc3 bc4
fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 80 * 80), []');

boundary = Boundary();
boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');

[bc1, bc2, bc3, bc4] = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);

boundary.setDirichlet(bc1);
boundary.setDirichlet(bc2);
boundary.setDirichlet(bc3);
boundary.setDirichlet(bc4);

dom = DOM(16);

% sweeping
tic;
dom.rayint(fem.Promoted.nodes, fem.Promoted.elems, fem.Promoted.neighbors);
toc;


end


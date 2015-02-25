fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 64 * 64), []');
boundary = Boundary();
boundary.set_boundary('(x - 1) * (y - 1) * x * y');
bc = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 1);
dom = DOM(16);
tic;
dom.rayint(fem.Promoted.nodes, fem.Promoted.elems, fem.Promoted.neighbors);
toc;

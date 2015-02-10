function hobj = focus_init(prec, min_area, sigma, freq)
% initialize all setups, without any PML.


addpath(genpath('~/Documents/github/femex'));
hobj.fem = FEM([-0.5, -.5, .5, -.5, .5, .5, -.5, .5]', ...
    prec, min_area);

hobj.ase = AssemblerExtension();


hobj.freq = freq;
hobj.sigma = sigma;

hobj.D = 1;
hobj.mu_a = 0.3;


% piecewise constant
hobj.rho = ones(size(hobj.fem.Promoted.elems, 2), 1);


hobj.M = hobj.fem.assema(1);
hobj.S = hobj.fem.assems(1);

hobj.Kernel = (hobj.freq * hobj.freq + 1i * hobj.freq * sigma) * hobj.M - hobj.S;

% boundary = Boundary(Dirichlet);
hobj.boundary = Boundary();
hobj.boundary.set_boundary('x - 0.5');
hobj.boundary.set_boundary('y - 0.5');
hobj.boundary.set_boundary('x + 0.5');
hobj.boundary.set_boundary('y + 0.5');


[bc1, bc2, bc3, bc4] = hobj.boundary.get_boundary(hobj.fem.Promoted.edges, hobj.fem.Promoted.nodes, 4);

hobj.boundary.setDirichlet(bc1);
hobj.boundary.setDirichlet(bc2);
hobj.boundary.setDirichlet(bc3);
hobj.boundary.setDirichlet(bc4);


hobj.Q = hobj.fem.assemlbc(1, bc1) + hobj.fem.assemlbc(1, bc2) + hobj.fem.assemlbc(1, bc3) + hobj.fem.assemlbc(1, bc4);


N = size(hobj.fem.Promoted.nodes, 2);
[hobj.dofs, hobj.ndofs] = hobj.boundary.dofs(N);
hobj.sol = zeros(N, 1);



ndofs = hobj.ndofs;
dofs = hobj.dofs;

hobj.DKernel = hobj.S + hobj.mu_a * hobj.M -0.5 * hobj.Q;

hobj.F = (hobj.DKernel(ndofs, ndofs) - hobj.DKernel(ndofs, dofs) * (hobj.DKernel(dofs, dofs)\hobj.DKernel(dofs, ndofs)));
hobj.G = hobj.DKernel(dofs, dofs)\hobj.DKernel(dofs, ndofs);




end

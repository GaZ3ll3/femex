function Energy = focus_solve(hobj, center)

Kernel = (hobj.freq * hobj.freq + 1i * hobj.freq * hobj.sigma) * hobj.M - hobj.S;


% using Bessel function as a experimental trial function for boundary condition 

ndofs = hobj.ndofs;
dofs  = hobj.dofs;

hobj.sol(ndofs) = besselj(0.,...
    hobj.freq * sqrt((hobj.fem.Promoted.nodes(1, ndofs) - center(1)).^2 + (hobj.fem.Promoted.nodes(2, ndofs) - center(2)).^2));

hobj.sol(dofs) = -Kernel(dofs, dofs)\(Kernel(dofs, ndofs) * hobj.sol(ndofs));

% figure;
% trimesh(hobj.fem.TriMesh', hobj.fem.Promoted.nodes(1,1:hobj.fem.Num_nodes),...
%     hobj.fem.Promoted.nodes(2, 1:hobj.fem.Num_nodes), conj(hobj.sol(1:hobj.fem.Num_nodes)).*hobj.sol(1:hobj.fem.Num_nodes));

Energy = conj(hobj.sol).* hobj.sol;

% mapping the solution onto all qnodes.
% hobj.qsol = focus_mapping(hobj.sol, hobj.fem.Promoted.elems, hobj.fem.Facet.Ref');

end


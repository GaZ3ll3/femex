function [LHS, RHS] = focus_diffusive(hobj, Energy)

dofs = hobj.dofs;
ndofs = hobj.ndofs;

interp_energy = hobj.sigma * focus_mapping(Energy, hobj.fem.Promoted.elems, hobj.fem.Facet.Ref');


[I, J ,V] = hobj.ase.assemble_ex_load_matrix(hobj.fem.Promoted.nodes, hobj.fem.Promoted.elems, hobj.fem.Facet.Ref, hobj.fem.Facet.Weights, interp_energy);

Ref = sparse(I, J ,V);

solution = hobj.DKernel\(Ref * hobj.rho);
RHS = solution(ndofs);


LHS = hobj.F\full(Ref(ndofs, :) - hobj.G'*Ref(dofs, :));

end


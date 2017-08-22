function forward_solve(o, p)
%FORWARD_SOLVE forward solve optical tomography.

    qD =...
        o.mapping(p.D, o.model.Promoted.elems, o.model.Facet.Ref');
    qsigma = ...
        o.mapping(p.sigma, o.model.Promoted.elems, o.model.Facet.Ref');

    % assembled matrix for finite element. 
    % mass matrix + stiff matrix + edge part.
    assemb_ot = o.model.assema(qsigma) + o.model.assems(qD) + (1.0/o.kappa) * o.matrices.Edge +...
        sqrt(-1) * o.omega * o.matrices.Mass;
    
    tic;
    u = assemb_ot \ o.loads;
    for l_id = 1:length(o.ndofs)
        o.J{l_id} =  u(o.ndofs, l_id) .*(1 + 0.00*( rand(size(o.ndofs, 1), 1) - 0.5));
    end
    toc;

end


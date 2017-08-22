function  forward_solve( o, p )
%FORWARD_SOLVE solves PDE for measurements.

    qD =...
        o.mapping(p.D, o.model.Promoted.elems, o.model.Facet.Ref');
    qsigma = ...
        o.mapping(p.sigma, o.model.Promoted.elems, o.model.Facet.Ref');

    % assembled matrix for finite element. 
    % mass matrix + stiff matrix + edge part.
    assemb_pat = o.model.assema(qsigma) + o.model.assems(qD) + (1.0/o.kappa) * o.matrices.Edge;
    
    tic;
    for l_id = 1:length(o.loads)
        u = assemb_pat\o.loads{l_id};
        % multipication noises
        o.H{l_id} = p.Gamma .* p.sigma .* u .*(1 + 0.00*( rand(o.tds, 1) - 0.5));
    end
    toc;

end


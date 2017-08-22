function [ f, g ] = objective_gradient( o, pv )
%OBJECTIVE_GRADIENT calculate objective and Frechet derivatives.

    N = o.tds;
    % load local parameters.
    o.local_param.sigma     = pv(1 : N);
    o.local_param.D         = pv(N + 1 : 2 * N);

    % construct adjoint solver.
    p = o.local_param;

    % interpolate D and sigma to quadrature points.
    qD =...
        o.mapping(p.D, o.model.Promoted.elems, o.model.Facet.Ref');

    qsigma = ...
        o.mapping(p.sigma, o.model.Promoted.elems, o.model.Facet.Ref');
    
    % assembly solver matrix.
    assemb_ot  = o.model.assema(qsigma) + o.model.assems(qD) + (1.0/o.kappa) * o.matrices.Edge +...
        sqrt(-1) * o.omega * o.matrices.Mass;
    
    
    current =  assemb_ot \ o.loads;
    
    
    for l_id = 1:length(o.ndofs)
        o.local_J{l_id} =  current(o.ndofs, l_id);
    end
    
    % initialize f and g.
    f = 0; g_s = zeros(N, 1); g_D = zeros(N, 1); 

        
    % regularization coefficients.
    alpha = 1e-3;
    beta = 1e-3;
    
    
    for l_id = 1:size(o.ndofs,1)
        mismatch = o.local_J{l_id} - o.J{l_id}; 
        
        % objective functionals      
        f= f   + 0.5 * (mismatch' * mismatch);
       
        
        load = zeros(N, 1);
        load(o.ndofs) = (mismatch);
        phi = conj(assemb_ot)\load;
        
        rp = real(phi); cp = imag(phi); 
        rc = real(current(:, l_id)); cc = imag(current(:, l_id)); 

        g_D = g_D - ...
            o.model.assemnode(rp, rc,  ones(N,1), zeros(N,1) ) - ...
            o.model.assemnode(cp, cc,  ones(N,1), zeros(N,1) );
        g_s = g_s - ...
            o.model.assemnode(rp, rc, zeros(N,1), ones(N,1) ) - ...
            o.model.assemnode(cp, cc, zeros(N,1), ones(N,1) );
        
    end
    
    
    
    
   reg = 0.5 * beta  * ( o.local_param.D' * o.matrices.Stiff * o.local_param.D ) + ...
          0.5 * alpha * ( o.local_param.sigma' * o.matrices.Stiff * o.local_param.sigma );
    
    f = f + reg;
    
    g_s =  g_s  + alpha * o.matrices.Stiff *o.local_param.sigma;
    g_D =  g_D  + beta * o.matrices.Stiff *o.local_param.D;
    
    
    g_D(o.ndofs) = 0;
    g_s(o.ndofs) = 0;
    
    g = [g_s; g_D];
end


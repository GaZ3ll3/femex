function [f, g] = objective_gradient(o, pv)
%OBJECTIVE_GRADIENT compute objective function and Frechet derivatives, 

    N = o.tds;
    % load local parameters.
    o.local_param.sigma     = pv(1 : N);
    o.local_param.D         = pv(N + 1 : 2 * N);
    
    % Gamma is not explicitly included in reconstruction.
    o.local_param.Gamma     = ones(N, 1);

    % construct adjoint solver.
    p = o.local_param;

    % interpolate D and sigma to quadrature points.
    qD =...
        o.mapping(p.D, o.model.Promoted.elems, o.model.Facet.Ref');

    qsigma = ...
        o.mapping(p.sigma, o.model.Promoted.elems, o.model.Facet.Ref');
    
    % assembly solver matrix.
    assemb_pat = o.model.assema(qsigma) + o.model.assems(qD) + (1.0/o.kappa) * o.matrices.Edge;

    % preallocate memory for solutions.
    fluence = cell(length(o.loads), 1);
    
    for l_id = 1:length(o.loads)
        fluence{l_id} = assemb_pat\o.loads{l_id};
        o.local_H{l_id} = p.Gamma .* p.sigma .* fluence{l_id};
    end

    % initialize f and g.
    f = 0; g_s = zeros(N, 1); g_D = zeros(N, 1); 

        
    % regularization coefficients.
    alpha = 1e-5;
    beta = 1e-5;
    
    for l_id = 1:length(o.loads)
        feedback = o.local_H{l_id}.* o.H{1} - o.H{l_id}.* o.local_H{1};
        
        % objective functionals incremetal.
        f = f + 0.5 * (feedback'* feedback);        
        
        % gradient for qPAT.                
        g_s = g_s + p.Gamma .* (fluence{l_id} .* o.H{1} - fluence{1} .* o.H{l_id}) .* feedback;
        
        % adjoint method for gradient
        phi = assemb_pat \ (p.sigma .* p.Gamma .* o.H{1} .* feedback);
        psi = assemb_pat \ (p.sigma .* p.Gamma .* o.H{l_id} .* feedback);
        
        % assemble gradient.
        g_D = g_D - o.model.assemnode(phi, fluence{l_id}, ones(N,1), zeros(N,1)) +...
            o.model.assemnode(psi, fluence{1}, ones(N,1), zeros(N,1));
        g_s = g_s - o.model.assemnode(phi, fluence{l_id}, zeros(N,1),  ones(N,1)) + ...
            o.model.assemnode(psi, fluence{1}, zeros(N,1),  ones(N,1));
    end
    
    
    normalization = 1/mean(o.H{1}.^4);
        
    reg = 0.5 * beta  * ( o.local_param.D' * o.matrices.Stiff * o.local_param.D ) + ...
          0.5 * alpha * ( o.local_param.sigma' * o.matrices.Stiff * o.local_param.sigma );
    
    f = f * normalization + reg;
    
    g_s =  g_s * normalization + alpha * o.matrices.Stiff *o.local_param.sigma;
    g_D =  g_D * normalization + beta * o.matrices.Stiff *o.local_param.D;
    
    
    g_D(o.ndofs) = 0;
    g_s(o.ndofs) = 0;
    g = [g_s; g_D];

end


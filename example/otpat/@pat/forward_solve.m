function [m, info] = forward_solve(this, p)
    % solves the forward problem
    % input: 
    %   p is parameter struct, contains at least 'D', 'sigma',
    %   'Gamma'
    % output:
    %   m is struct for measurement
    %   info shows utility information.

    % linear interpolation from node points to quadrature points.

    qD =...
        this.mapping(p.D, this.model.Promoted.elems, this.model.Facet.Ref');
    qsigma = ...
        this.mapping(p.sigma, this.model.Promoted.elems, this.model.Facet.Ref');

    % assembled matrix for finite element. 
    % mass matrix + stiff matrix + edge part.
    assemb_pat = this.model.assema(qsigma) + this.model.assems(qD) + (1.0/this.kappa) * this.matrices.Edge;
    assemb_ot = assemb_pat + sqrt(-1) * this.omega * this.matrices.Mass;
    
    tic;
    for l_id = 1:length(this.loads)
        u = assemb_pat\this.loads{l_id};
        
        v = assemb_ot \this.loads{l_id}; % complex
        % multipication noises
        m.H{l_id} = p.Gamma .* p.sigma .* u .*(1 + 0.00*( rand(this.geometry.N, 1) - 0.5));
        m.J{l_id} = v(this.geometry.ndofs)  .*(1 + 0.00*( rand(size(this.geometry.ndofs, 1),1) - 0.5));
    end
    t = toc;

    % timing information.
    info.t = t;
end
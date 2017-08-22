function [f,g] = objective_gradient(this, pv)
    % pv is parameter input in vector form. 
    % f is objective functional to minimize.
    % g is gradient w.r.t each parameters at current parameter
    % setting.


    % load local parameters.
    this.local_parameters = this.parameters;

    this.local_parameters.sigma     = pv(1:this.geometry.N);
    this.local_parameters.D = pv((this.geometry.N + 1): 2 * this.geometry.N);

    % construct adjoint solver.

    p = this.local_parameters;

    qD =...
        this.mapping(p.D, this.model.Promoted.elems, this.model.Facet.Ref');


    qsigma = ...
        this.mapping(p.sigma, this.model.Promoted.elems, this.model.Facet.Ref');

    % assembled matrix for finite element. 
    % mass matrix + stiff matrix + edge part.
    assemb_pat = this.model.assema(qsigma) + this.model.assems(qD) + (1.0/this.kappa) * this.matrices.Edge;
    assemb_ot  = assemb_pat +  sqrt(-1) * this.omega * this.matrices.Mass;
    
    % preprocessing for OT assemble matrix by Schur complement, L is
    % symmetric, so we do not have to transpose it in later use.

    ndofs = this.geometry.ndofs;
    N   = this.geometry.N;
    
    fluence = cell(length(this.loads), 1);
    current = cell(length(this.loads), 1);
    
    for l_id = 1:length(this.loads)
        fluence{l_id} = assemb_pat\this.loads{l_id};
        current{l_id} = assemb_ot \this.loads{l_id};
        this.local_measurements.H{l_id} = p.Gamma .* p.sigma .* fluence{l_id};
        this.local_measurements.J{l_id} = current{l_id}(ndofs);

    end

    % calculate objective and gradient     
    % gradient from PAT.
    obj_pat = 0;
    g_s = zeros(N, 1);
    g_D = zeros(N, 1); 

    % gradient from OT.
    obj_ot = 0;
    h_s = zeros(N, 1);
    h_D = zeros(N ,1);
    
    % regularization coefficients.
    alpha = 1e-2;
    beta = 0;
    
    % skip first measurement for quotient use.
    for l_id = 1:length(this.loads)
        feedback = this.local_measurements.H{l_id}.* this.measurements.H{1} - ...
            this.measurements.H{l_id}.* this.local_measurements.H{1};
        
        mismatch = this.local_measurements.J{l_id} - this.measurements.J{l_id}; 
        
        % objective functionals
        obj_pat = obj_pat + 0.5 * (feedback'* feedback);        
        obj_ot = obj_ot   + 0.5 * (mismatch' * mismatch);
        
        % gradient for PAT.                
        g_s = g_s + p.Gamma .* (fluence{l_id} .* this.measurements.H{1} - fluence{1} .* this.measurements.H{l_id}) .* feedback;
        
        phi = assemb_pat \ (p.sigma .* p.Gamma .* this.measurements.H{1} .* feedback);
        psi = assemb_pat \ (p.sigma .* p.Gamma .* this.measurements.H{l_id} .* feedback);
        
        g_D = g_D - this.model.assemnode(phi, fluence{l_id}, ones(N,1), zeros(N,1)) +...
            this.model.assemnode(psi, fluence{1}, ones(N,1), zeros(N,1));
        g_s = g_s - this.model.assemnode(phi, fluence{l_id}, zeros(N,1),  ones(N,1)) + ...
            this.model.assemnode(psi, fluence{1}, zeros(N,1),  ones(N,1));
        
        % gradient for OT.
        load = zeros(N, 1);
        load(ndofs) = conj(mismatch);
        phi = assemb_ot\load;
        
        h_D = h_D - this.model.assemnode(real(phi), real(current{l_id}),  ones(N,1), zeros(N,1) ) +...
            this.model.assemnode(imag(phi), imag(current{l_id}),  ones(N,1), zeros(N,1) );
        h_s = h_s - this.model.assemnode(real(phi), real(current{l_id}), zeros(N,1), ones(N,1) ) + ...
            this.model.assemnode(imag(phi), imag(current{l_id}), zeros(N,1), ones(N,1) );
        
    end

    reg = 0.5*alpha *( this.local_parameters.D' * this.matrices.Stiff * this.local_parameters.D ) + ...
        0.5*beta *( this.local_parameters.sigma' * this.matrices.Stiff * this.local_parameters.sigma );
    

    pp = 1/norm(this.measurements.H{1})^4;
    
    f = ...% pp * obj_pat +
    obj_ot + reg;
    
    grad_s =  h_s + beta * this.matrices.Stiff *this.local_parameters.sigma;
    grad_D =  h_D + alpha * this.matrices.Stiff *this.local_parameters.D;

    
    grad_D(this.geometry.ndofs) = 0;
    grad_s(this.geometry.ndofs) = 0;
    
%     disp(obj_pat);
%     disp( obj_ot);
%     disp(reg)

    g = [grad_s;grad_D];
end

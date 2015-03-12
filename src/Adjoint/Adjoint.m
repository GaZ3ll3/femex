classdef Adjoint < handle
    
    properties (Access = private)
        
        thres
        
        is_u_updated
        is_m_updated
        is_lambda_updated
    end
    
    properties (Access = public)
        J_func
        dJdu_func
        dJdm_func
        forward_func
        adjoint_func
        binary_op_func
        
        
        J
        dJdu
        dJdm
        gradient
        u
        m    
        lambda
    end

    methods
        function this = Adjoint(metadata)
            this.J_func = metadata{1};
            this.dJdu_func = metadata{2};
            this.dJdm_func = metadata{3};
            this.forward_func = metadata{4};
            this.adjoint_func = metadata{5};
            this.binary_op_func = metadata{6};
            this.m  = metadata{7};
            this.thres = metadata{8};
            this.is_m_updated = 0;
            this.is_u_updated = 0;
            this.is_lambda_updated = 0;
        end
        
        function forward(this)
            if (this.is_m_updated)
                this.u = this.forward_func(this.m);
                this.dJdu = this.dJdu_func(this.u, this.m);
                this.dJdm = this.dJdm_func(this.u, this.m);
                this.J = this.J_func(this.u, this.m);
                this.is_u_updated = 1;
                this.is_m_updated = 0;
            end
        end
        
        function adjoint(this)
            if (this.is_u_updated)
                this.lambda = this.adjoint_func(this.dJdu);
                this.is_lambda_updated = 1;
                this.is_u_updated = 0;
            end
        end
        
        function update(this, m0)     
            if (norm(this.m - m0, inf) > this.thres)
                this.m = m0;
                this.is_m_updated = 1;
            end
            
        end
        
        function grad(this)
            if (this.is_lambda_updated)
                this.gradient = this.binary_op_func(this.u, this.lambda);
                this.is_lambda_updated = 0;
            end
        end
    end
end
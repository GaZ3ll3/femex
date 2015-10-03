classdef radiative < handle
    %OPT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        fem 
        dom
        
        sol
        sol_
        source
        
        Mass
        Stiff
        alpha
        
        m_
        sigma_t_
        Q_
        Qt_
        
        
        
        sigma_a
        sigma_s
    
        sigma_a_0
        sigma_s_0
        sigma_s_
        sigma_a_
    end
    
    
    
    methods 
        function this = radiative()
            this.fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 40 * 40), []');
            this.dom = DOM(16);
            
            this.dom.rayint(this.fem.Promoted.nodes, this.fem.Promoted.elems, this.fem.Promoted.neighbors);
            this.source = this.source_fcn(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            
            this.sigma_a = this.Sigma_a_Fcn(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            this.sigma_s = this.Sigma_s_Fcn(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:))';
            
            this.alpha = 1e-8;
            this.Mass = this.fem.assema(1);
            this.Stiff =  this.fem.assems(1);
            
%             this.sigma_s_0 = this.sigma_s .* (1.0 + 0.02* (rand(size(this.sigma_s)) - 0.5));
            this.sigma_s_0 = 5.0 * ones(size(this.sigma_s));
            this.sigma_a_0 = 0.1 * ones(size(this.sigma_a));
            
            sigma_t = this.sigma_a + this.sigma_s;
            
            m = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, sigma_t);
            
            [this.sol,~,~,~,~]= gmres(speye(size(m, 1)) - m' * sparse(1:size(m,1), 1:size(m,1), this.sigma_s), m' * (this.source),10 , 1e-14, 400);
        end
        
        function f = objective(this, ret) 
            
            this.sigma_t_ = this.sigma_s + ret;
            this.m_ = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, this.sigma_t_);
            this.Q_ = speye(size(this.m_, 1)) - this.m_' * sparse(1:size(this.m_, 1), 1:size(this.m_, 1), this.sigma_s);
            this.Qt_ = speye(size(this.m_, 1)) - this.m_ * sparse(1:size(this.m_, 1), 1:size(this.m_, 1), this.sigma_s);
            [this.sol_, ~, ~, ~, ~] = gmres(this.Q_, this.m_' * (this.source), 10, 1e-14, 400);
            
            
            r = 0.5 * ret' * this.Stiff * ret;
            f = 0.5 * (this.sol_ - this.sol)' * this.Mass * (this.sol_ - this.sol) + this.alpha * r;
            
        end
        
        function g = gradient(this, ret)
            tmp = this.Mass * (this.sol_ - this.sol);
            h = this.Stiff * ret;
            [L, ~, ~,~,~] = gmres(this.Qt_, tmp, 10, 1e-14, 400); 
            
            g = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
                ...
                (this.source + ret .* this.sol_), L, this.sigma_t_) + this.alpha *h;
        end
        
        function [f, g] = objective_gradient(this, ret)
            
            sigma_t = this.sigma_s + ret;
            m = this.dom.si_build_omp(this.fem.Promoted.nodes, this.fem.Promoted.elems, sigma_t);
            
            Q = speye(size(m, 1)) - m' * sparse(1:size(m, 1), 1:size(m, 1), this.sigma_s);
            Qt = speye(size(m, 1)) - m * sparse(1:size(m, 1), 1:size(m, 1), this.sigma_s);
            
            
            this.sol_ = (Q)\(m' * (this.source));
            

            r = 0.5 * ret' * this.Stiff * ret;
            h = this.Stiff * ret;
            
            f = 0.5 * (this.sol_ - this.sol)' * this.Mass * (this.sol_ - this.sol) + this.alpha * r;
            tmp = this.Mass * (this.sol_ - this.sol);
             
            L = Qt\tmp; 
%             R = m * L;
            
%             g = this.sol_ .* R + this.alpha * h;
%             t = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
%                 ...
%                 (this.source + ret .* this.sol_), L, sigma_t);
            

%             g = g + t;
            
            g = this.dom.ray_scatter_grad(this.fem.Promoted.nodes, this.fem.Promoted.elems,...
                ...
                (this.source + ret .* this.sol_), L, sigma_t) + this.alpha *h;
        end
        
        

        function LBFGS(this, objstr, gradstr, cbstr)
            start = this.sigma_a_0;
            lb = 0. * ones(size(start));
            ub = Inf * ones(size(start));
            
            ret = lbfgsb(start,lb,ub, objstr, gradstr,...
                    [],cbstr,'maxiter',1e4,'m',8,'factr',1e-12,...
                    'pgtol',1e-12);
       
            this.sigma_a_ = ret;            
        end    

        
        function optimize(this)
            
            if (nargin == 1)
                start = this.sigma_a_0;
            end
            
            options = optimoptions('fminunc','Display','iter','Algorithm',...
            'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 800, 'TolFun',...
            1e-20, 'TolX',1e-20,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');

            problem.options = options;
            problem.x0 = start;
            problem.objective = @this.objective_gradient;
            problem.solver = 'fminunc';

            ret = fminunc(problem);
            
            this.sigma_a_ = ret;
            
            
        end
        
        function plot(this)
            figure(1);
            numofnodes = size(this.fem.Promoted.nodes, 2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(2);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sigma_a(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(3);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
            figure(4);
            trisurf(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
            this.fem.Promoted.nodes(2, 1:numofnodes), this.sol_(1:numofnodes),'EdgeColor', 'None');shading interp;
            colorbar;colormap jet;
        end
        
    end 
    
    methods (Static)
        
        function val = source_fcn(x, y)
            val = ring(x, y);
        end
        
        function val = Sigma_a_Fcn(x, y)
            val =(0.1 * (1.0 + 0.45 .* (x > 0.5)...
                .* (y > 0.5) .* (x < 0.75)...
                .* (y < 0.75) +...
                0.45 .* (x < 0.4) .*(y < 0.4)...
                .* (x > 0.2) .* (y > 0.2)));
        end
        
        function val = Sigma_s_Fcn(x, y)
            val = (5.0 * (1.0 + 0.45 .* (x > 0.5)...
                .* (y > 0.5) .* (x < 0.75)...
                .* (y < 0.75) +...
                0.45 .* (x < 0.4) .*(y < 0.4)...
                .* (x > 0.2) .* (y > 0.2)));
        end
        
    end
end


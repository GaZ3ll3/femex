classdef optics < handle
    % optics class solves inverse problem for photo-acoustic imaging.
    % the inverse problem is under diffusion approximation. 
    %
    %   -\nabla\cdot (D\nabla u) + \sigma u = 0
    %   u + k n\cdot D\nabla u = g on boundary.
    %
    %   data collected from
    %       H = \Gamma \sigma u
    %       J = u on boundary
    %
    %   this class holds several properties (public). 
    %
    %       1. parameters: D, sigma, Gamma are column vectors on quadrature
    %       points.
    %        
    %       2. kappa: equation related constant coefficient.
    %
    %       3. geoemtry: stores indices for boundary measurement.
    %
    %       4. measurements: H, J are cell array of column vectors for each
    %       sources, resp.
    %
    %       5. model: holds all necesary utility functions for finite
    %       element method.
    %
    %       6. sources: holds all source load.
    %
    %       7. result: holds information on reconstruction.
    %
    %
    %  References: 
    %       1 . H.Gao, K.Ren, H.Zhao, A hybrid reconstruction method for quantitative PAT, SIAM J. Imag. Sci., 6, 32-55, 2013
    
    
    
    properties (Access = private)
        matrices
        local_parameters
        local_measurements
    end
    properties (Access = public)
        measurements
        parameters
        kappa % kappa for Neumann condition.
        sources
        model
        loads % loads stay unchanged.
        geometry
        
        result
    end
    
    methods
        function this = optics()
            this.measurements = struct('H', {},  'J', {});
            this.parameters = struct('D', [], 'Gamma', [], 'sigma', []);
            this.matrices = struct('Edge', [], 'Mass',[], 'Stiff', []);
            this.geometry = struct('dofs', 0, 'ndofs', 0, 'N', 0);
            this.sources = {};
            this.loads = {};
            this.kappa = 1.0;
            
            this.result = {};
            
            % model utilizes finite element method on unit square mesh.
            this.model= FEM([0 0 1 0 1 1 0 1]', 1, 1.0/2/128/128, []', 3);
            bnd_eq = {'x-1', 'y-1', 'x', 'y'};
            
            boundary = Boundary();

            for b_id = 1:length(bnd_eq)
                boundary.set_boundary(bnd_eq{b_id});
            end

            [bcs{1}.bc, bcs{2}.bc, bcs{3}.bc, bcs{4}.bc] =...
                boundary.get_boundary(this.model.Promoted.edges,...
                this.model.Promoted.nodes, length(bnd_eq));
            

            % Edge matrix stay unchanged.
            this.matrices.Edge =...
                this.model.assemlbc(1, bcs{1}.bc) + ...
                this.model.assemlbc(1, bcs{2}.bc) + ...
                this.model.assemlbc(1, bcs{3}.bc) + ...
                this.model.assemlbc(1, bcs{4}.bc);
            
            this.matrices.Mass = this.model.assema(1);
            this.matrices.Stiff = this.model.assems(1);
            
            for b_id = 1:length(bnd_eq)
                bcs{b_id}.qnodes1D = ...
                    this.model.Assembler.qnodes1D(this.model.Promoted.nodes,...
                    this.model.Edge.Qnodes, bcs{b_id}.bc);
            end
            
            this.load_source_functions(@external_sources);
            
            for s_id = 1:length(this.sources)
                this.loads{s_id} = 0;
                for b_id = 1: length(bnd_eq)
                    cur_load = this.sources{s_id}(bcs{b_id}.qnodes1D');
                    this.loads{s_id} = this.loads{s_id}+ this.model.assemrbc(cur_load, bcs{b_id}.bc);   
                    
                end
            end
            
            this.load_parameter_functions(struct('D', @diffusionF, 'sigma', @absorptionF, 'Gamma', @GruneisenF));
            
            % extract boundary node indices.
            for b_id = 1:length(bnd_eq)
                boundary.setDirichlet(bcs{b_id}.bc);
            end
            this.geometry.N = size(this.model.Promoted.nodes, 2);
            [this.geometry.dofs, this.geometry.ndofs] = boundary.dofs(this.geometry.N);
            
        end
        
        function delete(this)
            this.model.delete()
        end
        
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
            assemb = this.model.assema(qsigma) + this.model.assems(qD) + (1.0/this.kappa) * this.matrices.Edge;
            
            
            
            tic;
            for l_id = 1:length(this.loads)
                u = assemb\this.loads{l_id};
                m.H{l_id} = p.Gamma .* p.sigma .* u .*(1 + 0.00*( rand(size(u)) - 0.5));
                m.J{l_id} = u(this.geometry.ndofs); % on boundary value.
                
            end
            t = toc;
            
            % timing information.
            info.t = t;
        end
        
        function [f,g] = objective_gradient(this, pv)
            % pv is parameter input in vector form. 
            % f is objective functional to minimize.
            % g is gradient w.r.t each parameters at current parameter
            % setting.
            
            
            % load local parameters.
            this.local_parameters = this.parameters;
            
            this.local_parameters.D = pv;
            
            % construct adjoint solver.
       
            p = this.local_parameters;
            
            qD =...
                this.mapping(p.D, this.model.Promoted.elems, this.model.Facet.Ref');
            
            
            qsigma = ...
                this.mapping(p.sigma, this.model.Promoted.elems, this.model.Facet.Ref');
          
            % assembled matrix for finite element. 
            % mass matrix + stiff matrix + edge part.
            assemb = this.model.assema(qsigma) + this.model.assems(qD) + (1.0/this.kappa) * this.matrices.Edge;
            
            fluence = cell(length(this.loads), 1);
            
            for l_id = 1:length(this.loads)
                fluence{l_id} = assemb\this.loads{l_id};
                this.local_measurements.H{l_id} = p.Gamma .* p.sigma .* fluence{l_id};
                this.local_measurements.J{l_id} = fluence{l_id}(this.geometry.ndofs); % on boundary value.
            end
            
            % calculate objective and gradient 
            
            f = 0;      
            g = zeros(this.geometry.N, 1);
            
            alpha = 1e-6;
           
            for l_id = 1:length(this.loads)
                feedback = this.local_measurements.H{l_id} - this.measurements.H{l_id};
                f = f+ 0.5 * (feedback'* feedback);
%                 g = g + p.Gamma .* fluence{l_id} .* feedback;            
                v = assemb \ (p.Gamma.*p.sigma.*feedback);
                g = g - this.model.assemnode(v, fluence{l_id}, ones(size(feedback)), zeros(size(feedback))  );
            end
            
            f = f + 0.5*alpha *( pv'* this.matrices.Stiff * pv);
            g = g + alpha * this.matrices.Stiff *pv;
            
            g(this.geometry.ndofs) = 0;
        end
        
        function backward_solve(this, initial_guess)
            if (nargin == 1)
                initial_guess = zeros(this.geometry.N, 1);
            end
            
            initial_guess(this.geometry.ndofs) = this.parameters.D(this.geometry.ndofs);
            
%             options = optimoptions('fminunc','Display','iter','Algorithm',...
%             'quasi-newton', 'HessUpdate', 'bfgs', 'GradObj','On', 'MaxIter', 1000, 'TolFun',...
%             1e-8, 'TolX',1e-5,'MaxFunEvals', 1e5, 'DerivativeCheck', 'Off');
% 
%             problem.options = options;
%             problem.x0 = initial_guess;
%             problem.objective = @this.objective_gradient;
%             problem.solver = 'fminunc';
% 
%             this.result = fminunc(problem);

            opts    = struct( 'factr', 1e7, 'pgtol', 1e-8, 'm', 50, 'x0', initial_guess, 'maxIts', 1e3, 'maxTotalIts', 1e5);
            opts.printEvery     = 1;

            [this.result, ~, ~] =...
                lbfgsb_c(@this.objective_gradient, zeros(this.geometry.N, 1), inf * ones(this.geometry.N, 1), opts);
            
        end
            
        
        function load_parameter_functions(this, func)
            this.parameters.D = func.D(this.model.Promoted.nodes');
            this.parameters.sigma = func.sigma(this.model.Promoted.nodes');
            this.parameters.Gamma = func.Gamma(this.model.Promoted.nodes');
        end
        
        function load_source_functions(this, func_holder) 
            funcs = func_holder();
            for f_id = 1:length(funcs)
                this.sources{f_id} = funcs{f_id};
            end
        end
        
        function visualize(this, u)
            figure;
            trisurf(this.model.TriMesh', ...
                this.model.Promoted.nodes(1,:)', ...
                this.model.Promoted.nodes(2,:)', ...
                u, 'EdgeColor', 'none' ); shading interp; colormap jet; colorbar;view(2);
        end
    end
    
    methods (Static)   
        function [interpolate] = mapping(func, elems, trans_ref)
            % interpolation mapping
            
            % allocate memory
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
    end
    
end


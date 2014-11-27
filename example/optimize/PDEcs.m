classdef PDEcs < handle
    % PDE constrained problem solver
    properties (Access = public)
        % settings (static), usually will not change during the iterations,
        % not including x0, thought x0 is also static
        static_var
        % settings (volatile), will change all the way until the end,
        % x is not included 
        volatile_var   
        % Ipopt options
        options         
        % Ipopt funcs setup
        funcs        
        % regularization
        beta
        % start point
        x0        
        % final solution
        x        
        % output information
        info
    end
    
    properties (Access = private)
        % current iteration
        iter
    end
    
    
    % methods depend on the functions will or won't change the object
    methods
        % static_var is always not relevant to the construction of object
        function this = PDEcs(auxdata)
            this.set_static_var(auxdata)

            
            this.set_options();
            this.set_target(1.0);
            
            this.set_funcs();
            this.init_volatile_var();
        end
        
        % one pass
        function set_static_var(this, auxdata)
            % static_var
            % *    numofnodes
            % *    S
            % *    Q
            % *    LoadVector
            % *    fem
            % *    data
            
            
            this.static_var.fem = FEM([0 0 1 0 1 1 0 1]', auxdata{1}, auxdata{2});
            
            toy_fem = FEM([0 0 1 0 1 1 0 1]', auxdata{1}, 0.5);
            
            this.static_var.numofnodes = this.static_var.fem.Num_nodes;
            
            % all Neumann boundary does not need to specify dofs
            
            this.static_var.S = this.static_var.fem.assems(1);
            this.static_var.M = this.static_var.fem.assema(1);
            this.static_var.Q = this.static_var.fem.assemlbc(1,this.static_var.fem.Promoted.edges);
            this.static_var.R = toy_fem.assema(1);
            
            f = auxdata{3};
            Load = f(this.static_var.fem.Qnodes);
  
            % suppose Neumann data as vanishing data.
            
            this.static_var.LoadVector = this.static_var.fem.asseml(Load);
            this.x0         = auxdata{4};
            this.beta       = auxdata{5};
            
                        
        end
        
        function set_target(this, var) 
            % var as true solution
            this.static_var.data  = -(this.static_var.S - var*this.static_var.M)\this.static_var.LoadVector;
        end
        
        function res = get_var_length(this)
            res = size(this.static_var.fem.Promoted.elems,2);
        end
        
        function set_options(this)
            
           

            % there are no constraints, thus no jac needed.
            this.options.ipopt.max_iter              = 400;
            this.options.ipopt.print_level           = 5;
            
            this.options.ipopt.hessian_approximation = 'limited-memory';
            this.options.ipopt.limited_memory_update_type = 'bfgs'; % sr1
            this.options.ipopt.limited_memory_max_history = 6;
            
            this.options.ipopt.mu_strategy           = 'adaptive';
            this.options.mu_oracle                   = 'probing';
            this.options.ipopt.tol                   = 1e-8;
            
            
%             this.options.ipopt.derivative_test       = 'first-order';
%             this.options.ipopt.derivative_test_tol   = 1e-7;
%             this.options.ipopt.derivative_test_perturbation = 1e-8;
%             this.options.ipopt.derivative_test_print_all = 'no';
        end
        
        function set_funcs(this)
            this.funcs.objective = @this.objective;
            this.funcs.gradient  = @this.gradient;
%             this.funcs.iterfunc  = @this.callback;
        end

        
        % preallocate memory for volatile_var
        function init_volatile_var(this)
            % declaration as placehold
            % assuming x0 as constant 1.
            this.volatile_var.M = ...
                this.static_var.fem.assema(this.x0);
            
        end
        
        
        function set_volatile_var(this, x)
            this.volatile_var.M = ...
                this.static_var.fem.assema(x);
        end
        
        
        
        
        
        function f = objective(this, x)
                        
            this.set_volatile_var(x);
            
            this.volatile_var.u = -(this.static_var.S - this.volatile_var.M)\this.static_var.LoadVector;
            
            diff = this.volatile_var.u - this.static_var.data;
            
            % with regularization term
            f = 0.5*diff'* this.static_var.Q*diff + 0.5*this.beta*x*x;
            
        end
        
        function g = gradient(this, x)
            
            this.set_volatile_var(x);
            
            this.volatile_var.u = -(this.static_var.S - this.volatile_var.M)\this.static_var.LoadVector;
            diff = this.volatile_var.u - this.static_var.data;
            this.volatile_var.v = (this.static_var.S - this.volatile_var.M)\(this.static_var.Q*diff);
            

            g = this.volatile_var.u'* this.static_var.M*this.volatile_var.v...
                + this.beta*x;
            
        end
        
%         function b = callback(this, nIter, f)
%             % only returns true or false
%         end
        
        function solve(this)
            [this.x, this.info] = ipopt(this.x0, this.funcs, this.options);
        end
    end
    
    methods (Static)
        
    end
    
end


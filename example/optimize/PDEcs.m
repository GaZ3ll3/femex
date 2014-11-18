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
        function obj = PDEcs()
            obj.set_static_var()
        end
        
        % one pass
        function set_static_var(obj)
            obj.static_var = 0;
            
            
        end
        
        function set_options(this, x0_var,  cl_var, cu_var)
            
            this.x0         = x0_var;
            this.options.cl = cl_var;
            this.options.cu = cu_var;
            
            this.options.ipopt.jac_c_constant        = 'yes';
            this.options.ipopt.hessian_approximation = 'limited-memory';
            this.options.ipopt.mu_strategy           = 'adaptive';
            this.options.ipopt.tol                   = 1e-7;
            this.options.ipopt.derivative_test       = 'first-order';

            this.options.ipopt.derivative_test_tol   = 1e-7;
            this.options.ipopt.derivative_test_perturbation = 1e-8;

            this.options.ipopt.derivative_test_print_all = 'yes';
        end

        
        % preallocate memory for volatile_var
        function init_volatile_var(obj)
            obj.volatile_var = 0;
        end
        
    end
    
    methods (Static)

    end
    
end


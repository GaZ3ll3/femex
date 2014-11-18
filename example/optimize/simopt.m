classdef simopt < handle
    %SIMOPT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        options
        funcs
        
        x
        info
    end
    
    methods
        function this = simopt(x0_var, cl_var, cu_var)
            this.x0 = x0_var;
            this.options.cl = cl_var;
            this.options.cu = cu_var;
        end
        
        function set_options(this)
              this.options.ipopt.jac_c_constant        = 'yes';
              this.options.ipopt.hessian_approximation = 'limited-memory';
              this.options.ipopt.mu_strategy           = 'adaptive';
              this.options.ipopt.tol                   = 1e-7;
              this.options.ipopt.derivative_test       = 'first-order';

              this.options.ipopt.derivative_test_tol   = 1e-7;
              this.options.ipopt.derivative_test_perturbation = 1e-8;

              this.options.ipopt.derivative_test_print_all = 'yes';
        end
        
        
        function set_func(this)
              this.funcs.objective         = @this.objective;
              this.funcs.constraints       = @this.constraints;
              this.funcs.gradient          = @this.gradient;
              this.funcs.jacobian          = @this.jacobian;
              this.funcs.jacobianstructure = @this.jacobian;
        end
        
        function solve(this)
            [this.x, this.info] = ipopt(this.x0,this.funcs,this.options);
        end
        
    end
    
    
    methods (Static)
        function f = objective (x)
              f = (x(1) - x(2))^2 + ...
                  (x(2) + x(3) - 2)^2 + ...
                  (x(4) - 1)^2 + (x(5) - 1)^2;
              
        end

% ----------------------------------------------------------------------
        function g = gradient (x)
            % if after objective, g is called multiple times with different
            % x as input param, then better with a assignment here, or do
            % not use assignment at all, simply remove this from member
            % functions.
          g = 2*[ x(1) - x(2);
              x(2) + x(3) - 2 - x(1) + x(2);
              x(2) + x(3) - 2;
              x(4) - 1;
              x(5) - 1 ];

        end

% ----------------------------------------------------------------------
        function c = constraints (x)
          c = [ x(1) + 3*x(2);
                x(3) + x(4) - 2*x(5);
                x(2) - x(5) ];

        end
% ----------------------------------------------------------------------
        function J = jacobian (~)  
          J = sparse([ 1  3  0  0  0;
                   0  0  1  1 -2;
                   0  1  0  0 -1 ]);

        end
    end
    
end


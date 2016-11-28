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
            this.set_target();
            this.set_funcs();
            this.init_volatile_var();
        end
        
        
        [] = set_static_var(this, auxdata);
        [] = set_target(this, var);
        [] = set_options(this);
        [] = set_funcs(this);
        
        
        

        

    end
    
    
end


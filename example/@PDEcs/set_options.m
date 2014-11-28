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


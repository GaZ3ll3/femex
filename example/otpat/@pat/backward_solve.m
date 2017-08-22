function backward_solve(this, initial_guess)
    if (nargin == 1)
        initial_guess_s = 0.1 *  ones(this.geometry.N, 1);
        initial_guess_D = 0.01 * ones(this.geometry.N, 1);
        
        initial_guess = [initial_guess_s; initial_guess_D];
    end
    
    % enforce boundary
    initial_guess(                  this.geometry.ndofs) = this.parameters.sigma(this.geometry.ndofs);
    initial_guess(this.geometry.N + this.geometry.ndofs) = this.parameters.D(this.geometry.ndofs);
    
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

    opts    = struct( 'factr', 1e-8, 'pgtol', 1e-12, 'm', 400, 'x0', initial_guess, 'maxIts', 2e2, 'maxTotalIts', 1e5);
    opts.printEvery     = 10;

    [this.result, ~, hist] =...
        lbfgsb_c(@this.objective_gradient, zeros(size(initial_guess)), inf * ones(size(initial_guess)), opts);
    disp(hist);
end

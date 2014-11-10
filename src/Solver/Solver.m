classdef Solver < handle

properties (Access = private)
id_ 
end

properties (Access = public)
ilu
gmres
end



methods
function this = Solver(type)
    if strcmp(type, 'ilu')
        this.ilu = 1;
        this.gmres= 0;
    elseif (strcmp(type, 'gmres'))
        this.ilu = 0;
        this.gmres = 1;
    else 
    	this.ilu = 0;
        this.gmres = 0;
    end
this.id_ = Solver_('new');
end


function delete(this)
Solver_('delete', this.id_);
end



function [x] = solve(this, A, b)

if (this.ilu == 1)
	% not parallelized, faster than MATLAB's gmres
    disp('solver uses ilu'); 
    options = AMGinit(A);
    [Prec, options] = AMGfactor(A, options);
    [x, ~] = AMGsolver(A, Prec, options,  b);

elseif (this.gmres == 1)
    % use gmres
    % parallelized with internal MATLAB
	disp('solver uses gmres');
    x = gmres(A, b,30, 1e-12, 4000);
else
	% parallelized with internal MATLAB
    disp('solver uses umfpack');
    x = A\b;
end


end
end


end
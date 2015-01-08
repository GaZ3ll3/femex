classdef Solver < handle

properties (Access = private)
	id_ 
end

properties (Access = public)
	ilu
	agmg
end



methods
function this = Solver(type)
	if strcmp(type, 'ilu')
		this.ilu = 1;
			this.agmg= 0;
	elseif (strcmp(type, 'agmg'))
		this.ilu = 0;
		this.agmg = 1;
	else 
		this.ilu = 0;
		this.agmg = 0;
	end
	this.id_ = Solver_('new');
end


function delete(this)
	Solver_('delete', this.id_);
end


function [DX, DY] = reference(this, points)
	[DX, DY] = Solver_('reference', this.id_, points);
end

function [gradX, gradY] = grad(this, soln, nodes, elems, DX, DY)
	[gradX, gradY] = Solver_('grad', this.id_, soln, nodes, elems, DX, DY);
end

function [x] = solve(this, A, b)

if (this.ilu == 1)
	% not parallelized, faster than MATLAB's gmres
    disp('solver uses ilu'); 
    options = AMGinit(A);
    [Prec, options] = AMGfactor(A, options);
    [x, ~] = AMGsolver(A, Prec, options,  b);
	Prec = AMGdelete(Prec);
elseif (this.agmg == 1)
    % use agmg
    % parallelized with internal MATLAB
	disp('solver uses agmg');
    [x] = agmg(A, b,30, 1e-12, 4000, 0);
else
	% parallelized with internal MATLAB
    disp('solver uses umfpack');
    x = A\b;
end


end
end


end

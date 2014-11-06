classdef Solver < handle

properties (Access = private)
id_ 
end

properties (Access = public)
ilu
umfpack
end



methods
function this = Solver(type)
    if (strcmp(type, 'ilu'))
        
        this.ilu = 1;
        this.umfpack = 0;
    end
    
    if (strcmp(type, 'umfpack'))
        this.ilu = 0;
        this.umfpack = 1;
    end
this.id_ = Solver_('new');
end


function delete(this)
Solver_('delete', this.id_);
end



function [x] = solve(this, A, b)

if (this.ilu == 1)
    disp('solver uses ilu');
    
    options = AMGinit(A);
    [Prec, options] = AMGfactor(A, options);
    [x, ~] = AMGsolver(A, Prec, options,  b);

elseif (this.umfpack == 1)
    disp('solver uses umfpack');
    x = A\b;
else
    x = A\b;
end


end
end


end
classdef Boundary < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Boundary(varargin)
  	if nargin == 1
    	this.id_ = Boundary_('new', varargin{1});
    elseif nargin == 0
    	this.id_ = Boundary_('placeholder');
    end
  end

  function setDirichlet(this, edges)
      Boundary_('set_dirichlet', this.id_, edges);
  end
  
  function delete(this)
      Boundary_('delete', this.id_);
  end
    
  
  function report(this)
      Boundary_('report', this.id_);
  end
  
  function [dof, ndof] = dofs(this, N)
      [dof, ndof]  = Boundary_('dofs', this.id_, N);
  end
  
  function [ndof] = getDirichlet(this)
      ndof = Boundary_('get_dirichlet', this.id_);
  end
  function [] = set_boundary(this, expr)
  	Boundary_('set_boundary', this.id_, expr);
  end
  
  function [varargout] = get_boundary(this, edges, nodes, NargOut)
     varargout = cell(NargOut, 1);
 
    [varargout{:}] = Boundary_('get_boundary', this.id_, edges, nodes);
  end
end
end

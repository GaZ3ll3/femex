classdef Boundary < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Boundary(edge)
    this.id_ = Boundary_('new', edge);
  end

  function delete(this)
  %DELETE Destructor.
    Boundary_('delete', this.id_);
  end
    
  
  function report(this)
      Boundary_('report', this.id_);
  end
  
  function dof = dofs(this, N)
      dof  = Boundary_('dofs', this.id_, N);
  end
  
  function [] = set_boundary(this, expr)
  	Boundary_('set_boundary', this.id_, expr);
  end
  
  function [bcs] = get_boundary(this, edges, nodes, NargOut)
     bcs = cell(NargOut, 1);
 
    [bcs{:}] = Boundary_('get_boundary', this.id_, edges, nodes);
  end
end
end

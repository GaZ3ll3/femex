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
  
  function dof = dofs(this, N, pedge)
      dof  = Boundary_('dofs', this.id_, N, pedge);
  end
  % Other methods...
end
end

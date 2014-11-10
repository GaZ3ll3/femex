classdef Integrator < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Integrator(dim, degree)
    assert(isnumeric(degree));
    this.id_ = Integrator_('new', dim, degree);
  end

  function delete(this)
  %DELETE Destructor.
    Integrator_('delete', this.id_);
  end
  
  function [qpts, qwts] = export(this)
  	[qpts, qwts] = Integrator_('export', this.id_);
  end
  
  
  % Other methods...
end
end

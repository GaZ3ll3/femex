classdef Tree < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Tree(position, size)
    this.id_ = Cell_('new', position, size);
  end

  function delete(this)
  %DELETE Destructor.
    Cell_('delete', this.id_);
  end
  
  function import(this, particles)
    Cell_('import', this.id_, particles);
  end
  
  function split(this)
    Cell_('split', this.id_); 
  end
  
  function m = buildmatrix(this, sigma_t, theta)
    m = Cell_('buildmatrix', this.id_, sigma_t, theta);
  end
  
end
end


classdef QuadTree < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = QuadTree(position, size)
    this.id_ = QuadTree_('new', position, size);
  end

  function delete(this)
  %DELETE Destructor.
    QuadTree_('delete', this.id_);
  end
  
  function import(this, particles)
    QuadTree_('import', this.id_, particles);
  end
  
  function split(this)
    QuadTree_('split', this.id_); 
  end
 
  function m = buildmatrix(this, sigma_t, theta)
    m = QuadTree_('buildmatrix', this.id_, sigma_t, theta);
  end
    
end
end


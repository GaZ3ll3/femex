classdef Treecode < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Treecode(x, y, l, v)
    this.id_ = Treecode_('new', x, y, l, v);
  end

  function delete(this)
  %DELETE Destructor.
    Treecode_('delete', this.id_);
  end
  
  function set(this, values)
    Treecode_('setAttribute', this.id_, values);
  end
  
  function split(this)
    Treecode_('split', this.id_); 
  end
 
  function m = buildmatrix(this, theta)
    m = Treecode_('buildmatrix', this.id_,  theta);
  end
    
end
end
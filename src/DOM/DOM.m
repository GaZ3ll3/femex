classdef DOM < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = DOM(nA)
    this.id_ = DiscreteOrinates_('new', nA);
  end

  function delete(this)
  %DELETE Destructor.
    DiscreteOrinates_('delete', this.id_);
  end
  
  function  rayint(this, nodes, elems, neighbors)
  	DiscreteOrinates_('rayint', this.id_, nodes, elems, neighbors);
  end
  
  function rayshow(this)
      DiscreteOrinates_('rayshow', this.id_);
  end
  
  function si_init(this, source, sigma_t, sigma_s)
      DiscreteOrinates_('si_init', this.id_, source, sigma_t, sigma_s);
  end
  
  % Other methods...
end
end

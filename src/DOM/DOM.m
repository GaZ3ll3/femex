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
  
  function [ret] = rayint(this, nodes, elems, neighbors, edges, weights, Fcn)
  	ret = DiscreteOrinates_('rayint', this.id_, nodes, elems, neighbors, edges, weights, Fcn);
  end
  
  
  % Other methods...
end
end


classdef DOM < handle

properties (Access = private)
  id_ % ID of the session instance.
  nAngle
end

methods
  function this = DOM(nA)
    this.id_ = DiscreteOrinates_('new', nA);
    this.nAngle = nA;
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
  
  function si_init(this)
      DiscreteOrinates_('si_init', this.id_);
  end
  
  function si_import(this, source, sigma_t, sigma_s)
      DiscreteOrinates_('si_import', this.id_, source, sigma_t, sigma_s);
  end
  
  function si_iter(this, nodes, elems)
      DiscreteOrinates_('si_iter', this.id_, nodes, elems);
  end
  
  function [ret] = si_output(this)
      [ret] = DiscreteOrinates_('si_output', this.id_);
  end
  
  function si_dsa(this, delta)
      DiscreteOrinates_('si_dsa', this.id_, delta);
  end
  
  function si_set(this, delta)
      DiscreteOrinates_('si_set', this.id_, delta);
  end
  
  function [T] = si_build(this, nodes, elems, sigma_t)
      T = DiscreteOrinates_('si_build', this.id_, nodes, elems, sigma_t);
      T = T/this.nAngle;
  end
  
  function [T] = si_build_omp(this, nodes, elems, sigma_t)
      T = DiscreteOrinates_('si_build_omp', this.id_, nodes, elems, sigma_t);
      T = T/this.nAngle;
  end
  
  function [T, M] = ray_build_omp(this, nodes, elems, sigma_t, sigma_s)
      [T, M] = DiscreteOrinates_('ray_build_omp', this.id_, nodes, elems, sigma_t, sigma_s);
      T = T/this.nAngle;
      M = M/this.nAngle;
  end
  
  function g = ray_scatter_grad(this, nodes, elems, u, v, sigma_t)
      g = DiscreteOrinates_('ray_scatter_grad', this.id_, nodes, elems, u, v, sigma_t);
  end
  % Other methods...
end
end


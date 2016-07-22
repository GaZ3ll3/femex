classdef MCRT < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
function this = MCRT(photons, n_nodes, mu_s, mu_t)
  	this.id_ = MCRT_('new', photons, n_nodes, mu_s, mu_t);
end

function delete(this)
  MCRT_('delete',this.id_);
end

function ret = simulate(this)
  ret = MCRT_('simulate', this.id_);
end

end
end

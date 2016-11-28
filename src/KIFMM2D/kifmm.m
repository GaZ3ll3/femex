classdef kifmm < handle

properties (Access = private)
    id_
end

methods
  function this = kifmm()
    [this.id_] = kifmm_('new');
  end
  
  function delete(this)
  %DELETE Destructor.
    kifmm_('delete', this.id_);
  end
  
  function ret = calcc(this, ncheb, charge, z, N, rank, mu_t) 
    ret =  kifmm_('calc_cache', this.id_,ncheb, charge, z, N, rank, mu_t);
  end
  
  function ret = calcf(this, ncheb, charge, z, N, rank, mu_t) 
    ret =  kifmm_('calc_fast', this.id_,ncheb, charge, z, N, rank, mu_t);
  end
  
  function debug(this)
      kifmm_('debug', this.id_);
  end
      
  
end

end


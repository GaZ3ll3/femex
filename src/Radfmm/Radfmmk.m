classdef Radfmmk < handle

properties (Access = private)
    id_
end

methods
  function this = Radfmmk()
    [this.id_] = Radfmmk_('new');
  end
  
  function delete(this)
  %DELETE Destructor.
    Radfmmk_('delete', this.id_);
  end
  
  function ret = calcc(this, ncheb, charge, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_cache', this.id_,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcf(this, ncheb,x, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_fast', this.id_,ncheb,x, z, N, m, mu_t);
  end
  
  function ret = calccs(this, ncheb, charge, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_cache_svd', this.id_,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfs(this, ncheb,x, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_fast_svd', this.id_,ncheb,x, z, N, m, mu_t);
  end
  
  function disp(this)
    Radfmmk_('disp', this.id_);
  end
  
end

end


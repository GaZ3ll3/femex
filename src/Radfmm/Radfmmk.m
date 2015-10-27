classdef Radfmmk < handle

properties (Access = private)
    id_
    id__
end

methods
  function this = Radfmmk()
    [this.id_] = Radfmmk_('new');
  end
  
  function delete(this)
  %DELETE Destructor.
    Radfmmk_('delete', this.id_);
  end
  
  function ret = calc(this, N, m, mu_t) 
    ret =  Radfmmk_('calc', this.id_, id_, N, m, mu_t);
  end
  
  function ret = calcc(this, ncheb, charge, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_cache', this.id_,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcf(this, ncheb,x, z, N, m, mu_t) 
    ret =  Radfmmk_('calc_fast', this.id_,ncheb,x, z, N, m, mu_t);
  end
end

end


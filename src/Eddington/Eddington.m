classdef Eddington < handle

properties (Access = private)
    id_uu
    id_uc
    id_us
    id_cc
    id_cs
    id_ss
end

methods
  function this = Eddington()
    [this.id_uu] = RadfmmUU_('new');
    [this.id_uc] = RadfmmUC_('new');
    [this.id_us] = RadfmmUS_('new');
    [this.id_cc] = RadfmmCC_('new');
    [this.id_cs] = RadfmmCS_('new');
    [this.id_ss] = RadfmmSS_('new');
  end
  
  function delete(this)
  %DELETE Destructor.
    RadfmmUU_('delete', this.id_uu);
    RadfmmUC_('delete', this.id_uc);
    RadfmmUS_('delete', this.id_us);
    RadfmmCC_('delete', this.id_cc);
    RadfmmCS_('delete', this.id_cs);
    RadfmmSS_('delete', this.id_ss);
  end
  
% UU
  function ret = calccUU(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmUU_('calc_cache', this.id_uu,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfUU(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmUU_('calc_fast', this.id_uu,ncheb,x, z, N, m, mu_t);
  end
% UC
  function ret = calccUC(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmUC_('calc_cache', this.id_uc,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfUC(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmUC_('calc_fast', this.id_uc,ncheb,x, z, N, m, mu_t);
  end
% US
  function ret = calccUS(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmUS_('calc_cache', this.id_us,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfUS(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmUS_('calc_fast', this.id_us,ncheb,x, z, N, m, mu_t);
  end
% CC
  function ret = calccCC(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmCC_('calc_cache', this.id_cc,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfCC(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmCC_('calc_fast', this.id_cc,ncheb,x, z, N, m, mu_t);
  end
% CS
  function ret = calccCS(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmCS_('calc_cache', this.id_cs,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfCS(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmCS_('calc_fast', this.id_cs,ncheb,x, z, N, m, mu_t);
  end 

% SS
  function ret = calccSS(this, ncheb, charge, z, N, m, mu_t) 
    ret =  RadfmmSS_('calc_cache', this.id_ss,ncheb, charge, z, N, m, mu_t);
  end
  
  function ret = calcfSS(this, ncheb,x, z, N, m, mu_t) 
    ret =  RadfmmSS_('calc_fast', this.id_ss,ncheb,x, z, N, m, mu_t);
  end 

  function disp(this)
    RadfmmUU_('disp', this.id_uu);
    RadfmmUC_('disp', this.id_uc);
    RadfmmUS_('disp', this.id_us);
    RadfmmCC_('disp', this.id_cc);
    RadfmmCS_('disp', this.id_cs);
    RadfmmSS_('disp', this.id_ss);
  end
  
end

end


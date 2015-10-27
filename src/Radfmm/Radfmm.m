classdef Radfmm < handle

properties (Access = public)
    id_
end

methods
  function this = Radfmm(ncheb,charge,loc, N, m)
    [this.id_] = Radfmm_('new', ncheb,charge,loc, N, m);
  end
  
  function delete(this)
  %DELETE Destructor.
    Radfmm_('delete', this.id_);
  end
end

end


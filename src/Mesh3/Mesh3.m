classdef Mesh3 < handle

  properties (Access = private)
    id_ % ID of the session instance.
  end

  methods
    function this = Mesh3(points, min_vol, PML)
      assert(isnumeric(points));
      if nargin == 2
      	this.id_ = Mesh3_('new', points, [],  min_vol);
      else
        this.id_ = Mesh3_('new', points, PML,  min_vol);
      end
    end

    function delete(this)
      Mesh3_('delete',this.id_);
    end
    
    function report(this)
      Mesh3_('report',this.id_);
    end
    
    function [points,connectivity] = meshdata(this)
        [points, connectivity] = Mesh3_('meshdata', this.id_);
    end
  end
end

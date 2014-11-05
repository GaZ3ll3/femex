classdef Mesh < handle
%DATABASE Hypothetical Matlab database API.
properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Mesh(boundary, min_area)
  %DATABASE Create a new database.
    assert(isnumeric(boundary));
    this.id_ = Mesh_('new', boundary, min_area);
  end

  function delete(this)
  %DELETE Destructor.
    Mesh_('delete', this.id_);
  end
  
  function clear(this)
      assert(isscalar(this));
      Mesh_('clear', this.id_);
  end

  function report(this, verbose)
      assert(isscalar(this));

      Mesh_('report', this.id_);
      
      
      
  end
  
  function [nodes, elems, segs] = promote(this, degree)
  %PUT Save something to the database.
    assert(isscalar(this));
    [nodes, elems, segs] = Mesh_('promote', this.id_, degree);
  end

  function [np, nt, ne] = export(this)
      [np, nt, ne] = Mesh_('export', this.id_);
  end
  
  
  % Other methods...
end
end

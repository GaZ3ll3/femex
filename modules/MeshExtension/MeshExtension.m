classdef MeshExtension < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = MeshExtension()
    this.id_ = MeshExtension_('new');
  end

  function delete(this)
    MeshExtension_('delete', this.id_);
  end
  
  
  function ret = layer(this, nodes, elems, boundary)
  	ret = MeshExtension_('layer', this.id_, nodes, elems, boundary);
  end
  
end

end

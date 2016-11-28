classdef Mesh < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
function this = Mesh(boundary, min_area, PML)
  assert(isnumeric(boundary));
  if nargin == 2
  	this.id_ = Mesh_('new', boundary, [],  min_area);
  else
    this.id_ = Mesh_('new', boundary, PML,  min_area);
  end
end

function delete(this)
  Mesh_('delete',this.id_);
end
  
function clear(this)
  assert(isscalar(this));
  Mesh_('clear', this.id_);
end

function report(this, verbose)
  assert(isscalar(this));
  Mesh_('report', this.id_); 
end
  
function [nodes, elems, segs, inds, neighbors] = promote(this, degree)
  assert(isscalar(this));
  [nodes, elems, segs, inds, neighbors] = Mesh_('promote', this.id_, degree);
end

function [np, nt, ne] = export(this)
  [np, nt, ne] = Mesh_('export', this.id_);
end
  
  
  % Other methods...
end
end

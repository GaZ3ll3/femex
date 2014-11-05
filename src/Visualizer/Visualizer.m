classdef Visualizer < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Visualizer()
    this.id_ = Visualizer__('new');
  end

  function delete(this)
  %DELETE Destructor.
    Visualizer_('delete', this.id_);
  end
  
end
end
classdef AssemblerExtension < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = AssemblerExtension()
    this.id_ = AssemblerExtension_('new');
  end

  function delete(this)
    AssemblerExtension_('delete', this.id_);
  end
  
  
  function [I, J ,V] = assemble_ex_grad_x_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  AssemblerExtension_('assemex_gradfunc_x',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V] = assemble_ex_grad_y_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  AssemblerExtension_('assemex_gradfunc_y',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V, W] = assemble_ex_grad_xy_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y)
     [I, J, V, W] =  AssemblerExtension_('assemex_gradfunc_xy',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y);
  end
  
  function [I, J ,V] = assemble_ex_load_matrix(this, pnodes, pelems, ref, weights, extern)
     [I, J, V] =  AssemblerExtension_('assemex_lm',  this.id_, pnodes, pelems, ref, weights, extern);
  end
  
  
end

end

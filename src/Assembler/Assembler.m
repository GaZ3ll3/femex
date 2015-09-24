classdef Assembler < handle

properties (Access = private)
  id_ % ID of the session instance.
end

methods
  function this = Assembler()
    this.id_ = Assembler_('new');
  end

  function delete(this)
  %DELETE Destructor.
    Assembler_('delete', this.id_);
  end
  
  function [F, DX, DY] = reference2D(this, nodes, qnodes)
  	[F, DX, DY] = Assembler_('reference2D', this.id_, nodes, qnodes);
  end
  
  function [F, DX] = reference1D(this, deg, qnodes)
  	[F, DX] = Assembler_('reference1D', this.id_, deg, qnodes);
  end
  
  function [I, J, V] = assema(this, pnodes, pelems, ref_fnk, weights, extern)
  	[I, J, V] = Assembler_('assema', this.id_, pnodes, pelems, ref_fnk, weights, extern);
  end
  
  function [I, J, V] = assems(this, pnodes, pelems, ref_gradx, ref_grady, weights, extern)
  	[I, J, V] = Assembler_('assems', this.id_, pnodes, pelems, ref_gradx, ref_grady, weights, extern);
  end
  
  function [L] = asseml(this, pnodes, qnodes, pelems, ref, weights, extern) 
    [L] = Assembler_('asseml', this.id_, pnodes, qnodes, pelems, ref, weights, extern);
  end
  
  function [L] = assemrbc(this, pnodes, qnodes, pedges, ref, weights, extern) 
    [L] = Assembler_('assemrbc', this.id_, pnodes, qnodes, pedges, ref, weights, extern);
  end
  
  function [I, J, V] = assemlbc(this, pnodes, pedges, ref, weights, extern) 
    [I, J, V] = Assembler_('assemlbc', this.id_, pnodes, pedges, ref, weights, extern);
  end
  
  function [C] = qnodes2D(this, pnodes, qnodes, pelems) 
  	[C] = Assembler_('qnodes2D', this.id_, pnodes, qnodes, pelems);
  end
  
  function [C] = qnodes1D(this, pnodes, qnodes, pedges) 
  	[C] = Assembler_('qnodes1D', this.id_, pnodes, qnodes, pedges);
  end
  
  function [I, J ,V] = assemble_grad_x_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  Assembler_('assemex_gradfunc_x',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V] = assemble_grad_y_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern)
     [I, J, V] =  Assembler_('assemex_gradfunc_y',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern);
  end
  
  function [I, J ,V, W] = assemble_grad_xy_func(this, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y)
     [I, J, V, W] =  Assembler_('assemex_gradfunc_xy',  this.id_, pnodes, pelems, ref, ref_gradx, ref_grady, weights, extern_x, extern_y);
  end
  
  function [I, J ,V] = assemble_load_matrix(this, pnodes, pelems, ref, weights, extern)
     [I, J, V] =  Assembler_('assemex_lm',  this.id_, pnodes, pelems, ref, weights, extern);
  end
  
  function [w] = assemble_elem(this, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v)
    w = Assembler_('assemloe', this.id_, pnodes, pelems, ref, refx, refy, weights, extern_s, extern_a, u, v);
  end
  
  % Other methods...
end
end

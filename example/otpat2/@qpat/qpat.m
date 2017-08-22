classdef qpat < handle
    %QPAT qpat solves inverse problem for quantitative photoacoustic
    %imaging, the inverse problem is under diffusion approximation for 2D.
    %
    %    -nabla\cdot (D\nabla u) + \sigma u = 0
    %      u + k n \cdot D\nabla u = g on boundary
    
    properties (Access = public)
        local_H
        local_param
        
        H     % H for measurements
        param % param: D, Gamma, sigma
        kappa % kappa for Robin boundary
        
        sources % Robin sources
        
        model % primary driver for assembly.
        
        loads % load vectors precomputed.
        
        ndofs % nodes not degrees of freedom.
        dofs  % nodes are degrees of freedom.
        tds   % total degrees in system.
        matrices; % matrices for reuse.
        
        result
    end
    
    methods
        function o = qpat(mdl, bnd_eq)
            % initialization.
            o.H = {};
            o.param = struct('D', [], 'Gamma', [], 'sigma', []);
            o.matrices = struct('Edge', [], 'Mass', [], 'Stiff', []);
            
            o.ndofs = 0; o.dofs = 0; o.tds = 0;
            
            o.kappa = 1.0;
            
            o.model = mdl;
            
            boundary = Boundary();

            for b_id = 1:length(bnd_eq)
                boundary.set_boundary(bnd_eq{b_id});
            end

            [bcs{1}.bc, bcs{2}.bc, bcs{3}.bc, bcs{4}.bc] =...
                boundary.get_boundary(o.model.Promoted.edges,...
                o.model.Promoted.nodes, length(bnd_eq));
            

            % precompute all matrices.
            o.matrices.Edge =...
                o.model.assemlbc(1, bcs{1}.bc) + ...
                o.model.assemlbc(1, bcs{2}.bc) + ...
                o.model.assemlbc(1, bcs{3}.bc) + ...
                o.model.assemlbc(1, bcs{4}.bc);
            
            o.matrices.Mass = o.model.assema(1);
            o.matrices.Stiff = o.model.assems(1);
            
            for b_id = 1:length(bnd_eq)
                bcs{b_id}.qnodes1D = ...
                    o.model.Assembler.qnodes1D(o.model.Promoted.nodes,...
                    o.model.Edge.Qnodes, bcs{b_id}.bc);
            end
            
            o.load_source_functions(@external_sources);
            
            for s_id = 1:length(o.sources)
                o.loads{s_id} = 0;
                for b_id = 1: length(bnd_eq)
                    cur_load = o.sources{s_id}(bcs{b_id}.qnodes1D');
                    o.loads{s_id} = o.loads{s_id}+ o.model.assemrbc(cur_load, bcs{b_id}.bc);   
                end
            end
            
            o.load_parameter_functions(struct('D', @diffusionF, 'sigma', @absorptionF, 'Gamma', @GruneisenF));
            
            % extract boundary node indices.
            for b_id = 1:length(bnd_eq)
                boundary.setDirichlet(bcs{b_id}.bc);
            end
            o.tds = size(o.model.Promoted.nodes, 2);
            [o.dofs, o.ndofs] = boundary.dofs(o.tds);
  
            
        end
        
        function delete(o)
            % nothing.
        end
        
        % methods in separated files
        [m, info] = forward_solve(o, p)
        [f,g] = objective_gradient(o, pv)
        load_parameter_functions(o, func)
        load_source_functions(o, func_holder) 
        visualize(o, p)
    end
    
    methods (Static)
        [interpolate] = mapping(func, elems, trans_ref)
    end
    
end


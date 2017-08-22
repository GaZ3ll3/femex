classdef qdot < handle
    %QDOT qdot solves inverse optical tomography in 2D.
    
    properties
        local_J
        local_param
        
        J     % H for measurements
        param % param: D, sigma
        kappa % kappa for Robin boundary
        omega % modularized frequency
        
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
        function o = qdot(qp)
            % load shared information from qpat.
            o.model = qp.model;
            o.ndofs = qp.ndofs; o.dofs = qp.ndofs; o.tds = qp.tds;
            o.matrices = qp.matrices;
            o.kappa = qp.kappa;
            
            % unshared information initialize.
            o.J = {};
            o.param = struct('D', [], 'sigma', []);
            o.omega = 100 * (2 * pi * 10^6) / (3 * 10^8);
            
            % copy synthetic coefficients into qdot.
            o.param.D = qp.param.D;
            o.param.sigma = qp.param.sigma;
            
            % create synthetic loads for each boundary node.
            % use sparse vector to save memory.
            tmp = eye(o.tds);
            o.loads = tmp(:, o.ndofs);
        end
        
        function delete(o)
            % nothing
        end
        
        [m, info] = forward_solve(o, p)
        [f,g] = objective_gradient(o, pv)
    end
    
    methods (Static)
        [interpolate] = mapping(func, elems, trans_ref)
    end
    
end


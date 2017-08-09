classdef pat < handle
    % pat class solves inverse problem for photo-acoustic imaging.
    % the inverse problem is under diffusion approximation. 
    %
    %      -\nabla\cdot (D\nabla u) + \sigma u = 0
    %      u + k n\cdot D\nabla u = g on boundary.
    %     
    % coupling with frequency domain OT.
    %   
    %      iw v/c -\nabla\cdot (D\nabla v) + \sigma v = 0 
    %      v + k n\cdot D\nabla v = h on boundary.      
    %
    %   data collected from
    %       H = \Gamma \sigma u
    %       J = v on boundary.
    %
    %   this class holds several properties (public). 
    %
    %       1. parameters: D, sigma, Gamma are column vectors on nodal
    %       points.
    %        
    %       2. kappa: equation related constant coefficient.
    %
    %       3. omega: frequency for OT.
    %
    %       4. geometry: stores indices for boundary measurement.
    %
    %       5. measurements: 'H' is cell array of column vectors for each
    %       source, resp.
    %
    %       6. model: holds all necesary utility functions for finite
    %       element method.
    %
    %       7. sources: holds all source load.
    %
    %       8. result: holds information on reconstruction.
    %
    %
    %  References: 
    %       1 . H.Gao, K.Ren, H.Zhao, A hybrid reconstruction method for quantitative PAT, SIAM J. Imag. Sci., 6, 32-55, 2013
    
    
    properties (Access = private)
        % these information should not be modified during computation, they
        % are used for updating.
        matrices
        local_parameters
        local_measurements
    end
    properties (Access = public)
        measurements
        parameters 
        kappa % kappa for Neumann condition.
        omega % frequency domain parameter.

        sources % provided from external sources for PAT.
                % also provided from external boundary conditions for OT.
        model % primary component for computing in finite element.
        loads % loads are computed from sources, which stay unchanged.
        geometry % geoemtry information for finite element, staty unchanged.
        
        result
    end
    
    methods
        function this = pat()
            this.measurements = struct('H', {}, 'J', {});
            this.parameters = struct('D', [], 'Gamma', [], 'sigma', []);
            this.matrices = struct('Edge', [], 'Mass',[], 'Stiff', []);
            this.geometry = struct('dofs', 0, 'ndofs', 0, 'N', 0);
            this.sources = {};
            this.loads = {};
            this.kappa = 1.0;
            this.omega = 200 * (2*pi * 10^6) / (3.0 * 10^8); % 200Hz
            
            this.result = {};
            
            % model utilizes finite element method on unit square mesh.
            % todo: use new version of meshing function.
            %
            this.model= FEM([0 0 1 0 1 1 0 1]', 2, 1.0/2/40/40, []', 6);
            bnd_eq = {'x-1', 'y-1', 'x', 'y'};
            
            boundary = Boundary();

            for b_id = 1:length(bnd_eq)
                boundary.set_boundary(bnd_eq{b_id});
            end

            [bcs{1}.bc, bcs{2}.bc, bcs{3}.bc, bcs{4}.bc] =...
                boundary.get_boundary(this.model.Promoted.edges,...
                this.model.Promoted.nodes, length(bnd_eq));
            

            % Edge matrix stay unchanged.
            this.matrices.Edge =...
                this.model.assemlbc(1, bcs{1}.bc) + ...
                this.model.assemlbc(1, bcs{2}.bc) + ...
                this.model.assemlbc(1, bcs{3}.bc) + ...
                this.model.assemlbc(1, bcs{4}.bc);
            
            this.matrices.Mass = this.model.assema(1);
            this.matrices.Stiff = this.model.assems(1);
            
            for b_id = 1:length(bnd_eq)
                bcs{b_id}.qnodes1D = ...
                    this.model.Assembler.qnodes1D(this.model.Promoted.nodes,...
                    this.model.Edge.Qnodes, bcs{b_id}.bc);
            end
            
            this.load_source_functions(@external_sources);
            
            for s_id = 1:length(this.sources)
                this.loads{s_id} = 0;
                for b_id = 1: length(bnd_eq)
                    cur_load = this.sources{s_id}(bcs{b_id}.qnodes1D');
                    this.loads{s_id} = this.loads{s_id}+ this.model.assemrbc(cur_load, bcs{b_id}.bc);   
                    
                end
            end
            
            this.load_parameter_functions(struct('D', @diffusionF, 'sigma', @absorptionF, 'Gamma', @GruneisenF));
            
            % extract boundary node indices.
            for b_id = 1:length(bnd_eq)
                boundary.setDirichlet(bcs{b_id}.bc);
            end
            this.geometry.N = size(this.model.Promoted.nodes, 2);
            [this.geometry.dofs, this.geometry.ndofs] = boundary.dofs(this.geometry.N);
            
        end        
        function delete(this)
            this.model.delete()
        end
        
        % methods in separated files
        [m, info] = forward_solve(this, p)
        backward_solve(this, initial_guess)
        [f,g] = objective_gradient(this, pv)
        load_parameter_functions(this, func)
        load_source_functions(this, func_holder) 
    end
    
     methods (Static)   
        [interpolate] = mapping(func, elems, trans_ref)
    end
    
end


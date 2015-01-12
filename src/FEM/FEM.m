classdef FEM < handle
    %FEM Calculated Finite Element Method
    %   Provided constructors and many methods
    
    properties (Access = public)
       
        Assembler
        Facet
        Edge
        Promoted
        
        TriMesh
        % all qnodes
        Qnodes
        
        Num_nodes
        Num_elems
        Num_edges
        
        Solution
        
        Ref_points
    end
    
    properties (Access = private)
        Ref_mesh
        Mesh
    end
    
    methods
        function this = FEM(edge_points, prec, min_area)
            this.Ref_mesh = Mesh([0 0 1 0 0 1]', 0.5);
            [this.Ref_points, ~, ~, ~] = this.Ref_mesh.promote(prec);
            
            this.Assembler = Assembler();
            
            this.Facet.Integrator = Integrator(2, 2*prec);
            [this.Facet.Qnodes, this.Facet.Weights] = this.Facet.Integrator.export();
            [this.Facet.Ref, this.Facet.RefX, this.Facet.RefY] = ...
                this.Assembler.reference2D(this.Ref_points, this.Facet.Qnodes);
        
            this.Edge.Integrator = Integrator(1, 2*prec);
            [this.Edge.Qnodes, this.Edge.Weights] = this.Edge.Integrator.export();
            [this.Edge.Ref, this.Edge.RefX] = ...
                this.Assembler.reference1D(prec, this.Edge.Qnodes);
              
            if size(edge_points,2) ~= 1
            	this.Mesh =  Mesh(edge_points', min_area);
            else
            	this.Mesh =  Mesh(edge_points, min_area);
            end
            [this.Num_nodes, this.Num_elems, this.Num_edges] = ...
                this.Mesh.export();
            
            
            [this.Promoted.nodes, this.Promoted.elems, this.Promoted.edges, this.Promoted.indices] = ...
                this.Mesh.promote(prec);
            
            this.TriMesh = this.Promoted.elems(1:3, :);
            this.Qnodes  = this.Assembler.qnodes2D(this.Promoted.nodes,...
                this.Facet.Qnodes, this.Promoted.elems);
            
        end
        
        function [Mass] = assema(this, Fcn)
            [I, J , U] = ...
                this.Assembler.assema(this.Promoted.nodes, ...
                this.Promoted.elems, this.Facet.Ref, this.Facet.Weights, Fcn);
            Mass = sparse(I, J ,U);
        end
        
        function [Stiff] = assems(this ,Fcn)
            [I, J ,V] = ...
                this.Assembler.assems(this.Promoted.nodes, ...
                this.Promoted.elems, this.Facet.RefX, this.Facet.RefY,...
                this.Facet.Weights, Fcn);
            
            Stiff = sparse(I,J , V);
        end
        
        function [Robin] = assemlbc(this, Fcn, BC)
            [I, J , W] = ...
                this.Assembler.assemlbc(this.Promoted.nodes, ...
                BC, this.Edge.Ref, this.Edge.Weights, Fcn);
            
            Robin = sparse(I, J ,W, size(this.Promoted.nodes, 2), ...
                size(this.Promoted.nodes, 2) );
        end
        
        function [LoadVector] = asseml(this, Fcn)
            LoadVector = this.Assembler.asseml(this.Promoted.nodes,...
                this.Facet.Qnodes, this.Promoted.elems, this.Facet.Ref, ...
                this.Facet.Weights, Fcn);
        end
        
        function [BCLoadVector] = assemrbc(this, Fcn, BC)
            BCLoadVector = this.Assembler.assemrbc(this.Promoted.nodes,...
                this.Edge.Qnodes, BC, this.Edge.Ref, this.Edge.Weights, Fcn);
        end
        
    end
    
end

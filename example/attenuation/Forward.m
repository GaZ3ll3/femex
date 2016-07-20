classdef Forward < handle
    properties
        fem
        g
        k
        sigma
        D
        u
        S
        M
        dofs
        ndofs
        phi
        f
        
    end
    
    properties(Access = private)
        ds
        db
        du

        sigma_nodes
    end
    
    methods
        function this = Forward()
            this.fem = FEM([-1, -1, 1, -1, 1, 1, -1, 1]', 1, 1/6400,[]');
            
        end
        
        function setDs(this, s, d)
            this.D = d;
            this.sigma =s;           
            
            this.sigma_nodes = this.sigma(this.fem.Promoted.nodes(1, :), this.fem.Promoted.nodes(2,:))';
            
            
            sigma_qnodes = this.sigma(this.fem.Qnodes(1,:), this.fem.Qnodes(2,:));
            D_qnodes = this.D(this.fem.Qnodes(1,:), this.fem.Qnodes(2,:));
            this.S = this.fem.assems(D_qnodes);
            this.M = this.fem.assema(sigma_qnodes);
            
        end
        
        function setg(this, func)
            this.g = func;
            N = size(this.fem.Promoted.nodes, 2);
            numofnodes = this.fem.Num_nodes;

            boundary = Boundary();
            this.u = zeros(N, 1);
            

            boundary.set_boundary('x - 1');
            boundary.set_boundary('y - 1');
            boundary.set_boundary('x + 1');
            boundary.set_boundary('y + 1');

            [bc1.bc, bc2.bc, bc3.bc, bc4.bc] = boundary.get_boundary(this.fem.Promoted.edges, this.fem.Promoted.nodes, 4);

            boundary.setDirichlet(bc1.bc);
            boundary.setDirichlet(bc2.bc);
            boundary.setDirichlet(bc3.bc);
            boundary.setDirichlet(bc4.bc);

            [this.dofs, this.ndofs] = boundary.dofs(N);
            tmp = this.g(this.fem.Promoted.nodes(1, :), this.fem.Promoted.nodes(2,:));
            this.u(this.ndofs) = tmp(this.ndofs);
            R = this.S + this.M;
            load = -R * this.u;
            this.u(this.dofs) = R(this.dofs, this.dofs)\load(this.dofs);
            
            figure(1);
            trimesh(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
                this.fem.Promoted.nodes(2, 1:numofnodes), this.u(1:numofnodes));
            
        end
        
        
        function setdb(this, db)
            this.db = db;
        end
        
        function setds(this, func)
            this.ds = func;
            N = size(this.fem.Promoted.nodes, 2);
            numofnodes = this.fem.Num_nodes;
            ds_nodes = this.ds(this.fem.Promoted.nodes(1,:), this.fem.Promoted.nodes(2,:));
            ds_qnodes = focus_mapping(ds_nodes, this.fem.Promoted.elems, this.fem.Facet.Ref');
            u_qnodes = focus_mapping(this.u, this.fem.Promoted.elems, this.fem.Facet.Ref');
            
            rhs = -ds_qnodes.* u_qnodes;
            
            load = this.fem.asseml(rhs);
          
            this.du = zeros(N, 1);
            this.du(this.ndofs) = 0;
            
            R = this.S + this.M;
            
            this.du(this.dofs) = R(this.dofs, this.dofs)\load(this.dofs);
            
            figure(2);
            trimesh(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
                this.fem.Promoted.nodes(2, 1:numofnodes), this.du(1:numofnodes));
            
            this.phi = this.u .* ds_nodes + this.du .* this.sigma_nodes;
            
            figure(3);
            trimesh(this.fem.TriMesh', this.fem.Promoted.nodes(1,1:numofnodes), ...
                this.fem.Promoted.nodes(2, 1:numofnodes), this.phi(1:numofnodes));
            
            
        end
        
        

        
        
        
    end
    
    
    methods(Static)
        function ret = helmholtz_solver2d(fem, load)
            
            ret = 0;
        end
    end
    
    
end


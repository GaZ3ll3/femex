% one pass
function set_static_var(this, auxdata)
    % static_var
    % *    numofnodes
    % *    S
    % *    Q
    % *    LoadVector
    % *    fem
    % *    data


    this.static_var.fem = FEM([0 0 1 0 1 1 0 1]', auxdata{1}, auxdata{2});

    toy_fem = FEM([0 0 1 0 1 1 0 1]', auxdata{1}, 0.5);

    this.static_var.numofnodes = this.static_var.fem.Num_nodes;

    % all Neumann boundary does not need to specify dofs

    this.static_var.S = this.static_var.fem.assems(1);
    this.static_var.M = this.static_var.fem.assema(1);
    this.static_var.Q = this.static_var.fem.assemlbc(1,this.static_var.fem.Promoted.edges);
    this.static_var.R = toy_fem.assema(1);

    f = auxdata{3};
    Load = f(this.static_var.fem.Qnodes);

    % suppose Neumann data as vanishing data.

    this.static_var.LoadVector = this.static_var.fem.asseml(Load);
    this.x0                    = auxdata{4};
    this.static_var.true       = auxdata{5};
    this.beta                  = auxdata{6};
    


end


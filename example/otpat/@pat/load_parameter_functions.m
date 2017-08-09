function load_parameter_functions(this, func)
    this.parameters.D = func.D(this.model.Promoted.nodes');
    this.parameters.sigma = func.sigma(this.model.Promoted.nodes');
    this.parameters.Gamma = func.Gamma(this.model.Promoted.nodes');
end
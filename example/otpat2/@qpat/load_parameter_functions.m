function load_parameter_functions(o, func)
%LOAD_PARAMETER_FUNCTIONS load all parameter functions and evaluate at all
%nodes.
    o.param.D = func.D(o.model.Promoted.nodes');
    o.param.sigma = func.sigma(o.model.Promoted.nodes');
    o.param.Gamma = func.Gamma(o.model.Promoted.nodes');
end


op = optics();
[op.measurements , info] = op.forward_solve(op.parameters);
op.backward_solve(ones(size(op.parameters.D)));
op.visualize(op.result);
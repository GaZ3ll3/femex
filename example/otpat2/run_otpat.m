function  [res, qp] = run_otpat(initial)
%RUN_OTPAT solves inverse problem of hybrid problem
global qPAT qDOT rp rd


    pat = 1; dot = 0;
    model= FEM([0 0 1 0 1 1 0 1]', 1, 1.0/2/30/30, []', 3);
    qPAT = qpat(model,  {'x-1', 'y-1', 'x', 'y'});
    qDOT = qdot(qPAT);
    
    qPAT.forward_solve(qPAT.param);
    qDOT.forward_solve(qDOT.param);
    
    if (nargin < 1)
        initial = [0.1 * ones(qPAT.tds, 1); 0.01 * ones(qPAT.tds, 1)];
    end
    
    if pat
        rp = 1;
    else
        rp = 0;
    end
    
    if dot 
        rd = 1;
    else
        rd = 0;
    end
    
    
    opts    = struct( 'factr', 1e2, 'pgtol', 1e-12, 'm', 400, 'x0', initial, 'maxIts', 8e2, 'maxTotalIts', 1e5);
    
    opts.printEvery     = 1;

    res = lbfgsb_c(@obj_grad, zeros(2*qPAT.tds, 1), inf * ones(2*qPAT.tds, 1), opts);
    
    qPAT.visualize(res(1:qPAT.tds));
    qPAT.visualize(res(qPAT.tds+1:end));
    
    qp = qPAT;
    
end

function [f,g]= obj_grad(pv)
global qPAT qDOT rp rd

    f = 0; g= zeros(size(pv));
    if rp 
        [qPAT_f, qPAT_g] = qPAT.objective_gradient(pv);
    end
    
    if rd 
        [qDOT_f, qDOT_g] = qDOT.objective_gradient(pv);
    end
    
    if rp 
        f = f+ qPAT_f;  
        g = g + qPAT_g;
    end
    
    if rd 
        f = f+ qDOT_f;
        g = g + qDOT_g;
    end
end

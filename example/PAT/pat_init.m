function hobj = pat_init(prec, min_area)

%addpath(genpath('../../'));
L = 0.1;
sigma0 = 80;

hobj.fem = FEM([-0.5, -.5, .5, -.5, .5, .5, -.5, .5]', ...
    prec, min_area, [-0.5-L, -0.5-L, .5+L, -.5-L, .5+L, .5+L, -.5-L, .5+L]');

hobj.sigma_x = @sigma_x;
hobj.sigma_y = @sigma_y;
hobj.mu_xy   = @mu_xy;


% all qnodes
sigma_x_of_qnodes     = hobj.sigma_x(hobj.fem.Qnodes(1,:), hobj.fem.Qnodes(2,:));
sigma_y_of_qnodes     = hobj.sigma_y(hobj.fem.Qnodes(1,:), hobj.fem.Qnodes(2,:));
sigma_xy_of_qnodes    = sigma_x_of_qnodes.* sigma_y_of_qnodes;
mu_xy_of_qnodes       = hobj.mu_xy(hobj.fem.Qnodes(1,:), hobj.fem.Qnodes(2,:));
mu_sigma_y_of_qnodes  = mu_xy_of_qnodes.* sigma_y_of_qnodes;
mu_sigma_x_of_qnodes  = mu_xy_of_qnodes.* sigma_x_of_qnodes;


hobj.M  = hobj.fem.assema(1);


% lump mass matrix
hobj.M = sparse(diag(sum(hobj.M)));


hobj.Sm = hobj.fem.assems(mu_xy_of_qnodes);
hobj.Mx = hobj.fem.assema(sigma_x_of_qnodes);
hobj.My = hobj.fem.assema(sigma_y_of_qnodes);
hobj.Mxy = hobj.fem.assema(sigma_xy_of_qnodes);

[I, J, V, W] = hobj.fem.Assembler.assemble_grad_xy_func(...
    hobj.fem.Promoted.nodes, hobj.fem.Promoted.elems, ...
        hobj.fem.Facet.Ref, hobj.fem.Facet.RefX, hobj.fem.Facet.RefY,...
        hobj.fem.Facet.Weights, 1, 1);

hobj.P  = sparse(J, I, V);
hobj.Q  = sparse(J, I, W);


[~, ~, Vx, Wy] = hobj.fem.Assembler.assemble_grad_xy_func(...
    hobj.fem.Promoted.nodes, hobj.fem.Promoted.elems, ...
        hobj.fem.Facet.Ref, hobj.fem.Facet.RefX, hobj.fem.Facet.RefY,...
        hobj.fem.Facet.Weights, mu_sigma_y_of_qnodes - mu_sigma_x_of_qnodes, ...
        mu_sigma_x_of_qnodes - mu_sigma_y_of_qnodes);
hobj.Px = sparse(J, I, Vx);
hobj.Qy = sparse(J, I, Wy);

    function ret = sigma_x(x, y)
        ret = sigma0 *(abs(x) > 0.5).*((abs(x) - 0.5)./L - sin(2*pi*(abs(x) - 0.5)./L)/2/pi);
    end

    function ret = sigma_y(x, y)        
        ret = sigma0 *(abs(y) > 0.5).*((abs(y) - 0.5)./L - sin(2*pi*(abs(y) - 0.5)./L)/2/pi);         
    end

    function ret = mu_xy(x, y)
        ret = 1.0 + 0.2 * sin(2 * pi* x) + 0.1 * cos(2*pi* y);
    end





end


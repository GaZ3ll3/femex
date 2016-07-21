function  [ret] = ray_diffusion(fem, sigma_a, sigma_s)

%fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 128 * 128), []');

boundary = Boundary();
boundary.set_boundary('x - 1');
boundary.set_boundary('y - 1');
boundary.set_boundary('x');
boundary.set_boundary('y');

[bc1, bc2, bc3, bc4] = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 4);

boundary.setDirichlet(bc1);
boundary.setDirichlet(bc2);
boundary.setDirichlet(bc3);
boundary.setDirichlet(bc4);


sigma_a_fcn = @(x, y) (sigma_a  +  0.0.*abs(cos(2*pi*x)));
sigma_s_fcn = @(x, y) (sigma_s  +  0.0.*abs(sin(2*pi*x)));



% center = [0.6, 0.4];
% radius = 0.2;

% source_fcn = @(x,y)(((x - center(1)).^2 + (y - center(2)).^2) <= radius^2) ...
%     .*(1 + cos(pi*sqrt((x - center(1)).^2 + (y - center(2)).^2)/radius));
% source_fcn = @(x, y) (1.0 .* (x > 0.25) .* (x < 0.75) .* (y > 0.6) .* (y < 0.8));
% source_fcn = @period;

source_fcn = @ring;


% solve diffusion equation
sigma_a_qnodes = sigma_a_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
sigma_s_qnodes = sigma_s_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
sigma_t_qnodes = sigma_a_qnodes + sigma_s_qnodes;

source_q_nodes = source_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));


Q = fem.assemlbc(1, bc1) + fem.assemlbc(1, bc2) + fem.assemlbc(1, bc3) + fem.assemlbc(1, bc4);
DSA_K = fem.assems((1/3)./sigma_t_qnodes) + fem.assema(sigma_a_qnodes) + 0.5 * Q;


load = fem.asseml(source_q_nodes);

ret = DSA_K\load;

figure(2)
trisurf(fem.TriMesh', fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), ret,...
    'EdgeColor','none','LineStyle','none','FaceLighting','phong');shading interp;colormap jet;
colorbar;view(2);


end


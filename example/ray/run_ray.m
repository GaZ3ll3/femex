fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 64 * 64), []');
boundary = Boundary();
boundary.set_boundary('(x - 1) * (y - 1) * x * y');
bc = boundary.get_boundary(fem.Promoted.edges, fem.Promoted.nodes, 1);
dom = DOM(128);
tic;
dom.rayint(fem.Promoted.nodes, fem.Promoted.elems, fem.Promoted.neighbors);
toc;

sigma_a_fcn = @(x, y) (0.1  + 0.0*abs(cos(2*pi*x)));
sigma_s_fcn = @(x, y) (0.4  + 0.0*abs(sin(2*pi*x)));
source_fcn = @(x,y) ( 1 .* (x > 0.45) .* (x < 0.55) .* (y > 0.45) .* (y <  0.55));

sigma_a = sigma_a_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));
sigma_s = sigma_s_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));

sigma_t = sigma_a + sigma_s;
source = source_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:));

dom.si_init(source, sigma_t, sigma_s, fem.Promoted.nodes, fem.Promoted.elems);
% dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);

res = zeros(40, 1);


for i = 1:40
pre = dom.si_output();
tic;
dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);
toc;
post = dom.si_output();
res(i) = norm(pre - post);
end



figure(2);
trisurf(fem.TriMesh', fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), post',...
    'EdgeColor','none','LineStyle','none','FaceLighting','phong');shading interp;





function run_ray()

fem = FEM([0 0 1 0 1 1 0 1]', 1, 1/(2 * 64 * 64), []');
dom = DOM(32);

% sweeping
tic;
dom.rayint(fem.Promoted.nodes, fem.Promoted.elems, fem.Promoted.neighbors);
toc;

% set the functions
sigma_a_fcn = @(x, y) (0.1  + 0.1*abs(cos(2*pi*x)));
sigma_s_fcn = @(x, y) (10  + 0.2*abs(sin(2*pi*x)));

center = [0.6, 0.4];
radius = 0.2;

source_fcn = @(x,y)(((x - center(1)).^2 + (y - center(2)).^2) <= radius^2) ...
    .*(1 + cos(pi*sqrt((x - center(1)).^2 + (y - center(2)).^2)/radius));

sigma_a = sigma_a_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));
sigma_s = sigma_s_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));

sigma_t = sigma_a + sigma_s;
source = source_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:));


% initialize
dom.si_init(source, sigma_t, sigma_s, fem.Promoted.nodes, fem.Promoted.elems);

res = zeros(80, 1);


% source iteration
for i = 1:80
pre = dom.si_output();
dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);
post = dom.si_output();
res(i) = norm(pre - post);
end




% result
figure(1)
plot(log(res));

figure(2);
trisurf(fem.TriMesh', fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), post',...
    'EdgeColor','none','LineStyle','none','FaceLighting','phong');shading interp

end



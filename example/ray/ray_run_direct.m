function [ret] = ray_run_direct(dom, fem, sigma_a_coef, sigma_s_coef)

sigma_a_fcn = @(x, y) (sigma_a_coef  + 0.0.*abs(cos(2*pi*x)));
sigma_s_fcn = @(x, y) (sigma_s_coef  + 0.0.*abs(sin(2*pi*x)));

% center = [0.6, 0.4];
% radius = 0.2;

% source_fcn = @(x,y)(((x - center(1)).^2 + (y - center(2)).^2) <= radius^2) ...
%     .*(1 + cos(pi*sqrt((x - center(1)).^2 + (y - center(2)).^2)/radius));

% source_fcn = @(x, y) (1.0 .* (x > 0.25) .* (x < 0.75) .* (y > 0.6) .* (y < 0.8));

source_fcn = @ring;

% source_fcn = @period;


sigma_a = sigma_a_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));
sigma_s = sigma_s_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2, :));

sigma_t = sigma_a + sigma_s;

source = source_fcn(fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:));


m = dom.si_build(fem.Promoted.nodes, fem.Promoted.elems, sigma_t);
I = eye(size(m,1));

tic;
ret = (I - m' * diag(sigma_s))\(m' * (source'));
toc;
trisurf(fem.TriMesh', fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), ret,...
'EdgeColor','none','LineStyle','none','FaceLighting','phong');shading interp;



end


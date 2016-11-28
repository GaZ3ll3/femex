function [post] = ray_run( dom, fem, sigma_a_coef, sigma_s_coef )

global bc1 bc2 bc3 bc4
% set the functions
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


% solve diffusion equation
sigma_a_qnodes = sigma_a_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
sigma_s_qnodes = sigma_s_fcn(fem.Qnodes(1,:), fem.Qnodes(2,:));
sigma_t_qnodes = sigma_a_qnodes + sigma_s_qnodes;


tic;
Q = fem.assemlbc(1, bc1) + fem.assemlbc(1, bc2) + fem.assemlbc(1, bc3) + fem.assemlbc(1, bc4);
DSA_K = fem.assems((1/3)./sigma_t_qnodes) + fem.assema(sigma_a_qnodes) + 0.5 * Q;
toc;

% initialize

dom.si_init();


% import

dom.si_import(source, sigma_t, sigma_s);

% source iteration


% first run

dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);

pre = dom.si_output();
dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);
%
post = dom.si_output();
diff = post - pre;
sdiff = sigma_s' .* diff; 
tmp_load = focus_mapping(sdiff, fem.Promoted.elems, fem.Facet.Ref');
LoadVector = fem.asseml(tmp_load);
delta = DSA_K\LoadVector;
dom.si_dsa(delta);
%
post = dom.si_output;

err = norm(pre - post);

counter = 1;

fprintf('-------------------------------------------------------------------\n');
fprintf('|   iteration      |           error       |      covergence      |\n');
fprintf('-------------------------------------------------------------------\n');
fprintf('|   %6.2d         |    %12.8f       |                      |\n', counter ,  err);

while (err > 1e-6)
    counter = counter + 1;
    pre = post;
    dom.si_iter(fem.Promoted.nodes, fem.Promoted.elems);   
    %
    post = dom.si_output();
    diff = post - pre;
    sdiff = sigma_s' .* diff; 
    tmp_load = focus_mapping(sdiff, fem.Promoted.elems, fem.Facet.Ref');
    LoadVector = fem.asseml(tmp_load);
    delta = DSA_K\LoadVector;
    dom.si_dsa(delta);
    %
    post = dom.si_output;
    err_ = norm(post - pre);
    fprintf('|   %6.2d         |    %12.8f       |     %12.8f     |\n', counter ,  err_,  err/err_);
    err = err_;
end


toc;

figure(2)
trisurf(fem.TriMesh', fem.Promoted.nodes(1,:), fem.Promoted.nodes(2,:), post,...
    'EdgeColor','none','LineStyle','none','FaceLighting','phong');shading interp

colormap('jet');

end


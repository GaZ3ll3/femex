function ret = treecode(n, theta)

global tree s theta_

theta_ = theta;

x = 1/n:2/n:(n-1)/n;
y = 1/n:2/n:(n-1)/n;
z = zeros(2, n/2 * n/2);
for i = 1:n/2
    for j = 1:n/2
        z(1, n/2 * (i - 1) + j) = x(i);
        z(2, n/2 * (i - 1) + j) = y(j);
    end
end

attr = sigma_t(z(1,:), z(2,:));
num = size(x, 2)^2;

fprintf('\n1. initializing treecode root ...');
tic;
tree = Treecode(0, 0, 1, log2(n) - 1);
t1 = toc;

fprintf('done, using %f\n2. initializing attributes of all points ...', t1);
tic;
tree.set(attr);
t1 = toc ;

fprintf('done, using %f\n3. populating treecode ...', t1);
tic;
tree.split();
t1 = toc ;

s = sigma_s(z(1,:), z(2,:));
f = ring(z(1,:), z(2,:));

fprintf('done, using %f\n4. applying K(x, y) to source term ... ',t1);
tic;
% m = tree.buildmatrix(theta);
tree.preprocess(f'/(2 * pi));
p = tree.apply(theta, f'/(2*pi));
t1 = toc;


fprintf('done, using %f\n5. solving linear system using GMres ... ', t1);
tic;
ret=gmres(@forward, p, 30, 1e-12, 30,[],[], p);
t1 =toc ;
fprintf('done, using %f\n', t1);
[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;

end

function lhs = forward(rhs)
global tree theta_ s
tree.preprocess(s'.*rhs/(2 * pi));
lhs = rhs - tree.fast_apply(theta_, s'.*rhs/(2 * pi));
end

function val = sigma_s(x, y)
    val = 5.0 * ones(size(x));
    %val = 10.0 * (1 + x.^2 + y.^2);
end

function val = sigma_t(x, y)
    %val = 5.2 * ones(size(x));
    val = 0.2 + sigma_s(x, y);
end



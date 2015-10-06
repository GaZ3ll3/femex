function ret = treecode(n, theta)

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

fprintf('done, using %f\n4. building matrix of %d x %d  ...',t1, num, num);
tic;
m = tree.buildmatrix(theta);
t1 = toc;

s = sigma_s(z(1,:), z(2,:));
f = ring(z(1,:), z(2,:));

p = m * f'/(2*pi);

fprintf('done, using %f\n5. solving linear system using GMres ... ', t1);
tic;
[ret, ~, ~, ~, ~] = gmres(eye(size(m, 1)) -  m * sparse(1:num, 1:num, s) /(2 * pi), p, 10, 1e-12);
t1 =toc ;
fprintf('done, using %f\n', t1);
[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;

end

function val = sigma_t(x, y)
    val = 5.1 * ones(size(x));
end

function val = sigma_s(x, y)
    val = 5.0 * ones(size(x));
end

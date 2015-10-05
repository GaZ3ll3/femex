function [ret] = quadtree(n, sigma_t, sigma_s, theta)
x = 1/n:2/n:(n-1)/n;
y = 1/n:2/n:(n-1)/n;
z = zeros(2, n/2 * n/2);
for i = 1:n/2
    for j = 1:n/2
        z(1, n/2 * (i - 1) + j) = x(i);
        z(2, n/2 * (i - 1) + j) = y(j);
    end
end
x = [z; ones(1, n/2 * n/2)];


tic;
tree = Tree([0, 0], 1);
toc;

tic;
tree.import(x);
toc;

tic;
tree.split();
toc;

tic;
m = tree.buildmatrix(sigma_t, theta);
toc;

f = ring(x(1,:), x(2,:));
p = m' * f'/(2*pi);

ret = gmres(eye(size(m, 1)) - sigma_s/(2*pi) * m', p, 10, 1e-12);

[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;


end

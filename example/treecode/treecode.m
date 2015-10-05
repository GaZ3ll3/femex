function treecode(n, theta)
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

tic;
tree = Treecode(0, 0, 1, 6);
toc;

tic;
tree.set(attr);
toc;

tic;
tree.split();
toc;

tic;
m = tree.buildmatrix(theta);
toc;

s = sigma_s(z(1,:), z(2,:));
f = ring(z(1,:), z(2,:));
p = m * f'/(2*pi);


ret = gmres(eye(size(m, 1)) -  m * sparse(1:num, 1:num, s) /(2 * pi), p, 10, 1e-12);

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
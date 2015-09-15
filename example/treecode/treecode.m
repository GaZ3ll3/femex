function [ret] = treecode(n)
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

tree = Tree([0, 0], 1);
tree.import(x);
tree.split();
m = tree.buildmatrix();

f = ring(x(1,:), x(2,:));
p = m * f'/(2 * pi);

ret = gmres(eye(size(m, 1)) - 2.0/(2*pi) * m, p, 5, 1e-10);

[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2))


end


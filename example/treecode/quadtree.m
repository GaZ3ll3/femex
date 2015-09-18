function quadtree(n, sigma_t,theta)
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


tree = QuadTree([0, 0], 1);
tree.import(x);
tree.split();
m = tree.buildmatrix(sigma_t, theta);


end

function treecode(n)
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
tree = QuadTree([0, 0], 1);
toc;

tic;
tree.import(x);
toc;

tic;
tree.split();
toc;

end


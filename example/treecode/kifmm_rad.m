function [ret] = kifmm_rad(n, np, rank)
global m N z ncheb mu_s mu_t ker

ncheb = np;

m = rank;

x = 1/n:2/n:(n-1)/n;
y = 1/n:2/n:(n-1)/n;
z = zeros(2, n/2 * n/2);
for i = 1:n/2
    for j = 1:n/2
        z(1, n/2 * (i - 1) + j) = x(i);
        z(2, n/2 * (i - 1) + j) = y(j);
    end
end

N = size(x, 2)^2;
                                                    
charge = ring(z(1,:), z(2,:));
%charge = ones(N, 1);
% mu_s = 10.0 * (1 + z(1,:).^2 + z(2, :).^2)';
mu_s = 5.0 * ones(N, 1);

mu_t = 0.2 + mu_s;


tic;
ker = kifmm();
t = toc;

fprintf('1. Initialization kernel ... with time %f\n', t);

tic;
rhs = mu_s .*ker.calcc(ncheb, charge, z, N, m, mu_t);

t= toc;

fprintf('2. Caching necessary kernel evaluations ... with time %f\n', t);
% 

tic;
[ret,~,~,~,~] = gmres(@forward, rhs, 40, 1e-12, 40, [],[], rhs);

ret = ret./mu_s;

t = toc;
fprintf('3. GMRES takes time %f\n', t);

[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;view(2);


% ker.debug();
% export_fig result.png

end

function ret = forward(x)
global m N z ncheb mu_s mu_t ker


ret = ker.calcf(ncheb,x, z, N, m, mu_t);


ret = x - mu_s .* ret;


end

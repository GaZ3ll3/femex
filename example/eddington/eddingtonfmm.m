function [ret] = eddingtonfmm(n, nc, t)

global m N z ncheb mu_s mu_t kernels tau

tau = t;
ncheb = nc;

m = 1;

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

%mu_s = 10.0 * (1 + z(1,:).^2 + z(2, :).^2)';
lambda = (t/3) / (1 - t/3);
mu_s = 5.0 * ones(N, 1) * (1 - lambda^2);

mu_t = 0.2 + mu_s;

tic;
kernels = Eddington();
t = toc;

fprintf('1. Initialization kernel ... with time %f\n', t);

tic;
rhs_u = mu_s .* kernels.calccUU(ncheb, charge, z, N, m,mu_t);
rhs_c = mu_s .* kernels.calccUC(ncheb, charge, z, N, m, mu_t);
rhs_s = mu_s .* kernels.calccUS(ncheb, charge, z, N, m, mu_t);

[~] = kernels.calccCC(ncheb, charge, z, N, m , mu_t);
[~] = kernels.calccCS(ncheb, charge, z, N, m , mu_t);
[~] = kernels.calccSS(ncheb, charge, z, N, m , mu_t);

rhs = [rhs_u; rhs_c; rhs_s];
t = toc;

kernels.disp();
fprintf('2. Caching necessary kernel evaluations ... with time %f\n', t);

tic;
[ret, ~, ~,~,~] = gmres(@forward, rhs, 30, 1e-12, 30, [], [], rhs);

l = size(mu_s, 1);

ret(1:l) = ret(1:l)./mu_s;
ret(l+1:2*l) = ret(l+1:2*l)./mu_s;
ret(2*l+1:3*l) = ret(2*l+1:3*l) ./mu_s;

t = toc;
fprintf('3. GMRES takes time %f\n', t);

[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret(1:l), n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;view(2);

end

function ret = forward(x)

global m N z ncheb mu_s mu_t kernels tau
l = size(x, 1)/3;

ret = zeros(3*l, 1);


u = x(1:l);
c = x(l+1:2*l);
s = x(2*l+1:3*l);

u_ret = ...
    kernels.calcfUU(ncheb, u, z, N, m, mu_t) + ...
    tau * kernels.calcfUC(ncheb, c, z, N, m, mu_t) + ...
    tau * kernels.calcfUS(ncheb, s, z, N, m, mu_t);

c_ret = ...
    kernels.calcfUC(ncheb, u, z, N, m, mu_t) +...
    tau * kernels.calcfCC(ncheb, c, z, N, m, mu_t) +...
    tau * kernels.calcfCS(ncheb, s, z, N, m, mu_t);

s_ret = ...
    kernels.calcfUS(ncheb, u, z, N, m, mu_t) +...
    tau * kernels.calcfCS(ncheb, c, z, N, m, mu_t) +...
    tau * kernels.calcfSS(ncheb, s, z, N, m , mu_t);

ret(1:l) = u - mu_s .* u_ret;
ret(l + 1: 2 * l) = c - mu_s .* c_ret;
ret(2 *l + 1: 3 *l) = s - mu_s.* s_ret;



end

























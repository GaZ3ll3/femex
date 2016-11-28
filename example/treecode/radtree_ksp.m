function [ret] = radtree_ksp(n, nc)
% expriment for preconditioned gmres.

global m N Np zp z ncheb mu_s mu_sp mu_tp mu_t ker kerp

ncheb = nc;
m = 1;

[z, charge, N] = build_grid(n);

mu_s = 10.0 * (1 + z(1,:).^2 + z(2, :).^2)';
% mu_s = 5.0 * ones(N, 1);
mu_t = 0.2 + mu_s;


% preconditioner
tic;
kerp = Radfmmk();
[zp, chargep, Np] = build_grid(64);
% some approx coefficients
mu_sp = 10.0 * ones(Np, 1);
mu_tp = 0.2 + mu_sp;
rhsp= mu_sp.*kerp.calcc(4, chargep, zp, Np, m, mu_tp);
[retp,~,~,~,~] = gmres(@forwardp, rhsp, 30, 1e-5, 30, [],[], rhsp);

% retp shoul not be scaled.

% enhance retp to higher level.

grid_per_box = n/64;
ret_ap = zeros(n/2*n/2, 1);

for i = 1:n/2
    for j = 1:n/2
        ret_ap(n/2 * (i - 1) + j) = retp(...
            32 * floor((i-1)/grid_per_box) + floor((j - 1)/grid_per_box) + 1);
    end
end

t = toc;
fprintf('1. Initialization perconditioner ... with time %f\n', t);
% true operator

tic;
ker = Radfmmk();
t = toc;
fprintf('2. Initialization kernel ... with time %f\n', t);

tic;
rhs = mu_s .* ker.calcc(ncheb, charge, z, N, m, mu_t) - ...
    (ret_ap - mu_s.*ker.calcf(ncheb, ret_ap, z, N, m, mu_t));
t = toc;
fprintf('3. Caching necessary kernel evaluations ... with time %f\n', t);

tic;
[ret,~,~,~,~] = gmres(@forward, rhs, 30, 1e-12, 30, [], [], rhs);
t = toc;
fprintf('4. GMRES takes time %f\n', t);
ret = (ret + ret_ap)./mu_s;

end

function [z, charge, N]  =build_grid(n)
    x = 1/n:2/n:(n-1)/n;
    y = 1/n:2/n:(n-1)/n;
    z = zeros(2, n/2*n/2);
    for i = 1:n/2
        for j = 1:n/2
            z(1, n/2 * (i - 1) + j) = x(i);
            z(2, n/2 * (i - 1) + j) = y(j);
        end
    end

    N = size(x, 2)^2;
    charge = ring(z(1,:),z(2,:));
end

function ret = forward(x)
    global m N z ncheb mu_s mu_t ker
    ret = x - mu_s .* ker.calcf(ncheb, x, z, N, m, mu_t);
end

function ret = forwardp(x)
    global m Np zp mu_sp mu_tp kerp
    ret = x - mu_sp.* kerp.calcf(4, x, zp, Np, m, mu_tp);
end




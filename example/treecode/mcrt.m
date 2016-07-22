function [ret] = mcrt(np, n, mu_s, mu_t)

m = MCRT(np, n, mu_s, mu_t);
ret = m.simulate();

n = n * 2;

[X, Y] = meshgrid(1/n : 2/n:(n-1)/n);
surf(X,Y,reshape(ret, n/2,n/2), 'EdgeColor','None');
shading interp;colorbar; colormap jet;view(2);

end


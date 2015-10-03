function  callback_radiative(t, f, x)
global p
fprintf('%3d\t\t%0.8g\t\t%0.8g \n', t, f, norm(x - p.sigma_a)/norm(p.sigma_a));
end


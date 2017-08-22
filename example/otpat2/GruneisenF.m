function [ret] = GruneisenF(x)
 r = sqrt((x(:,1) - 0.7).^2 + (x(:,2) - 0.7).^2);
 ret = 0.5 + 0.15 * (1 + cos(r * 2 * pi/0.4)) .* (r < 0.2);
end


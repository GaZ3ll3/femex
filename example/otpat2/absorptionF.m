function [ret] = absorptionF(x)
    r = sqrt((x(:,1) - 0.3).^2 + (x(:,2) - 0.3).^2);
    ret = 0.2 + 0.2 * (1 + cos(r * 2 * pi/0.4)) .* (r < 0.2);
end


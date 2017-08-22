function [ret] = GruneisenF(x)
    ret = 0.1 + 0.2 * (x(:,1) > 0.6) .* (x(:,1) < 0.8) .* (x(:, 2) < 0.6) .* (x(:,2) > 0.2);
end


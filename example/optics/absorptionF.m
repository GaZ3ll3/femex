function [ret] = absorptionF(x)
    ret = 1.0 + (x(:,1) > 0.3) .* (x(:,1) < 0.5) .* (x(:, 2) < 0.6) .* (x(:,2) > 0.2);
end


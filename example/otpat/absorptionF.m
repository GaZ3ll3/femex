function [ret] = absorptionF(x)
    ret = 0.1 + 0.1*(x(:,1) > 0.3) .* (x(:,1) < 0.5) .* (x(:, 2) < 0.6) .* (x(:,2) > 0.2);
%     ret = 0.02 + 0.18 * external_phantom(x);
end


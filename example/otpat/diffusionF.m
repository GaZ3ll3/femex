function [ret] = diffusionF(x)
%     ret = 1 + exp(-(x(:,1)).^2 - (x(:,2)).^2);
    ret = 1.0;
    
    for i = 1:3
        for j = 1:3
            ret = ret + ...
        2.0 * (x(:,1) > (2 * i) * 0.1) .*...
        (x(:,1) < (2 * i+ 1) * 0.1) .* (x(:, 2) < ( 2 * j + 1) * 0.1) .* (x(:,2) > (2 * j) * 0.1);
        end
    end 
%     ret = 3.0 + 3.0 * external_phantom(x);
end


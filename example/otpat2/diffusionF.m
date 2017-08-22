function [ret] = diffusionF(x)
    r = sqrt((x(:,1) - 0.5).^2 + (x(:,2) - 0.5).^2);
    ret = 0.01* (1 + 0.1 * cos(r * pi/0.6).^3.* (r < 0.3)).^2 ;
%     ret = 0.01 * (1 + 0.5 * sin(pi * x(:,1)) .* sin(pi * x(:,2)));
%     ret = 1.0;
%     
%     for i = 1:3
%         for j = 1:3
%             ret = ret + ...
%         2.0 * (x(:,1) > (2 * i) * 0.1) .*...
%         (x(:,1) < (2 * i+ 1) * 0.1) .* (x(:, 2) < ( 2 * j + 1) * 0.1) .* (x(:,2) > (2 * j) * 0.1);
%         end
%     end 
%     ret = ret * 0.01;
%     ret = 3.0 + 3.0 * external_phantom(x);
end


function [ret] = period( x, y )

center = [0.7, 0.7];
center_2 = [0.3, 0.3];
radius = 0.2;

ret_1 = (((x - center(1)).^2 + (y - center(2)).^2) <= radius^2) .*(1 + cos(pi*sqrt((x - center(1)).^2 + (y - center(2)).^2)/radius));

ret_2 = (((x - center_2(1)).^2 + (y - center_2(2)).^2) <= radius^2) .*(1 + cos(pi*sqrt((x - center_2(1)).^2 + (y - center_2(2)).^2)/radius));

ret = ret_1 + ret_2;
end


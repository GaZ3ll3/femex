function ret  = ring( x, y )

center = [0.6, 0.4];

radius_1 = 0.15;
radius_2 = 0.25;


ret = 1.0 .* (((x - center(1)).^2 + (y - center(2)).^2) <= radius_2^2) .* (((x - center(1)).^2 + (y - center(2)).^2) >= radius_1^2) ;


% ret = x + y;

end


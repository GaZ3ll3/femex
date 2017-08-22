function [ret] = external_sources(d, s)
    l = @(x)( (abs(x(:,1)    ) < eps) );
    r = @(x)( (abs(x(:,1) - 1) < eps) );
    b = @(x)( (abs(x(:,2)    ) < eps) );
    t = @(x)( (abs(x(:,2) - 1) < eps) );
    
    
    ret = cell(12,  1);
%     q = sqrt(s/d);
    q = 2.0;
    theta = 2 * pi / length(ret);
    for s_id = 1:length(ret)  
        c = cos(s_id * theta);
        s = sin(s_id * theta);
        ret{s_id} = @(x)( (exp(q * (c * (x(:,2) - 0.5)  + s* (x(:,1) - 0.5))) .* (l(x) + b(x) + r(x) + t(x))));
    end
    

   
end


function [ret] = external_sources()
    l = @(x)( (abs(x(:,1)    ) < eps) );
    r = @(x)( (abs(x(:,1) - 1) < eps) );
    b = @(x)( (abs(x(:,2)    ) < eps) );
    t = @(x)( (abs(x(:,2) - 1) < eps) );
    
    vals = {l, r, b, t};
    
    ret = cell(30,  1);
    for s_id = 1:length(ret) 
        dir = mod(s_id, 4) + 1;
        frq = floor(s_id / 4);
        if dir > 2
            ret{s_id} = @(x)(cos(2 * frq * pi * x(:, 2)) .* vals{dir}(x));
        else
            ret{s_id} = @(x)(cos(2 * frq * pi * x(:, 1)) .* vals{dir}(x));
        end
    end
end


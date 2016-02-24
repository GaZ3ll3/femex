function [Z] = time_reversal(hobj, Y, dt, tl)

% hobj handle stores all geometry information
% Y stores the measurement


% extract inside information on edge

edges = hobj.fem.Promoted.edges;
len_edges = size(edges, 2);

el = 0;
for i = 1 : len_edges
    e = edges(:,i);
    if is_on_edge(e)
        el = el + 1;
    end
end

boundary = zeros(2 * el, 1);
eli = 1;
for i = 1:len_edges
    e = edges(:, i);
    if is_on_edge(e)
        boundary(eli) = e(1);
        boundary(eli + 1) = e(2);
        eli = eli + 2;
    end
end

boundary = unique(boundary);


% extract information on inside nodes

nodes = hobj.fem.Promoted.nodes;
len_nodes = size(nodes, 2);

nl = 0;
for i = 1 : len_nodes
    n = nodes(:, i);
    if is_in_domain(n)
        nl = nl + 1;
    end 
end

inside = zeros(nl, 1);
nli = 1;
for i = 1: len_nodes
    n = nodes(:, i);
    if is_in_domain(n)
        inside(nli) = i;
        nli = nli + 1;
    end
end

inner_part = setdiff(inside, boundary);

% reconstruction using time reversal
length = size(inner_part, 1);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'NormControl', 'on');

T = size(Y, 1);

tmp1 = Y(T, 1:len_nodes)';
tmp2 = Y(T, len_nodes + 1:2*len_nodes)';


% back propagate
q0 = tmp1(inner_part);
p0 = -tmp2(inner_part);

[~, Z] = ode113(@rigid, 0:dt:(T - 1)*dt, [q0;p0], options); 
 
fprintf('inf error of time reversal is %f\n', norm(Y(1, inner_part) - Z(T,1:size(Z, 2)/2), inf));


    function dy = rigid(t, y)
        
        n = floor((tl - t)/dt);
        alpha = (tl - t)/dt - n;
        n = n + 1;
        
        if alpha == 0
            rhs = Y(n, boundary);
        elseif alpha == 1
            rhs = Y(n + 1, boundary);
        else
            rhs = Y(n, boundary) *(1-alpha) + alpha * Y(n+1, boundary);
            
        end
        
        if (n ~= T)
            lhs = (Y(n + 1, boundary + len_nodes) - Y(n, boundary + len_nodes))/dt;
        else
            lhs = (Y(n , boundary + len_nodes) - Y(n - 1, boundary + len_nodes))/dt;
        end
        
        
        dy = zeros(2 * length, 1);
        u = y(1:length);
        v = y(length + 1: 2 * length);
        
        ul = v;
        vl = -hobj.M(inner_part, inner_part)\(hobj.Sm(inner_part, inner_part) * u + ...
            hobj.Sm(inner_part, boundary) * (rhs)' + hobj.M(inner_part, boundary) * lhs');
        
        dy(1:length) = ul;
        dy(length + 1: 2 * length) = vl;
        
        
    end


    function ret =  is_on_edge(e)
        node = hobj.fem.Promoted.nodes(:,e(1));
        if node(1) <= 0.51 && node(1) >= -0.51 && node(2)>= -0.51 && node(2) <= 0.51
            ret = 1;
        else
            ret = 0;
        end
    end

    function ret = is_in_domain(n)
        EPS = 1E-8;
        if n(1) <= 0.50 +EPS && n(1)>= -0.50 -EPS && n(2) >= -0.50 -EPS && n(2) <= 0.50 +EPS
            ret = 1;
        else
            ret = 0;
        end
    end


end


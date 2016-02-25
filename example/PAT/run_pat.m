

    h = pat_init(1, 1/20000);
    l = size(h.fem.Promoted.nodes,2);
    
    u0 = zeros(l , 1);
    ut0 = zeros(l, 1);


    dt = 0.005;
    tl = 1.0;

    for i = 1:l

        x = h.fem.Promoted.nodes(1, i);
        y = h.fem.Promoted.nodes(2, i);

        r = sqrt(x^2 + y^2);
        if r <= 1/4
            u0(i) = cos(4*pi*r) + 1.0;
        end

    end
    tic;
    [~, Y] = pat_solver(h, u0, ut0, dt, tl);
    toc;

    tic;
    [boundary, inner_part, Z] = time_reversal(h, Y, dt, tl);
    toc;
    
    residue = zeros(l, 1);
    residue(inner_part) = Y(1, inner_part)' - Z(end, 1:size(Z,2)/2)';
    
    disp(norm(residue, 2)/norm(Y(1, inner_part)', 2));
    
    % Z is now the first guess
    
    result = zeros(l, 1);
    result(inner_part) = Z(end, 1:size(Z,2)/2)';
    
    inc = Z(end, 1:size(Z,2)/2)';
    
    for j = 1:10
    
        next = zeros(l, 1);
        next(inner_part) = inc;

        tic;
        [~, Y1] = pat_solver(h, next, ut0, dt, tl);
        toc;

        tic;
        [~, ~, Z1] = time_reversal(h, Y1, dt, tl, boundary, inner_part);
        toc;

        inc = inc - Z1(end, 1:size(Z,2)/2)';
        
        result(inner_part) = result(inner_part) + inc;
        
        disp(norm(result(inner_part) - Y(1, inner_part)',2)/norm(Y(1, inner_part)', 2));
    end

    
    
    

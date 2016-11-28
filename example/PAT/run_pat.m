function run_pat()

    h = pat_init(1, 1/120000);
    l = size(h.fem.Promoted.nodes,2);
    
%     u0 = zeros(l , 1);
    ut0 = zeros(l, 1);


    dt = 0.002;
    tl = 1.50;

    
    u0 = pat_phantom(h, l);
    
    for i = 1:l
        
        
        x = h.fem.Promoted.nodes(1, i);
        y = h.fem.Promoted.nodes(2, i);
       
        
        
        if (x + 0.35)^2 + (y +0.35)^2 <= 0.01
            u0(i) = 0.5 * (cos(10 * pi *sqrt((x + 0.35)^2 + (y +0.35)^2 )) + 1.0);
        end
        
        if (x + 0.35)^2 + (y - 0.35)^2 <= 0.01
            u0(i) = 0.5 * (cos(10 * pi *sqrt((x + 0.35)^2 + (y -0.35)^2 )) + 1.0);
        end

        if (x - 0.35)^2 + (y - 0.35)^2 <= 0.01
            u0(i) = 0.5 * (cos(10 * pi *sqrt((x - 0.35)^2 + (y -0.35)^2 )) + 1.0);
        end
        
        if (x - 0.35)^2 + (y + 0.35)^2 <= 0.01
            u0(i) = 0.5 * (cos(10 * pi *sqrt((x - 0.35)^2 + (y +0.35)^2 )) + 1.0);
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
    
    err = norm(result(inner_part) - Y(1, inner_part)',2)/norm(Y(1, inner_part)', 2);
    
    save('data', 'result', 'Y', 'err');
    
    disp(err);
        
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
        
        err = norm(result(inner_part) - Y(1, inner_part)',2)/norm(Y(1, inner_part)', 2);
        disp(err);
        
        filename = sprintf('data-%d', j);
        save(filename, 'result', 'err');
        
    end

end
    
    
    
    
    

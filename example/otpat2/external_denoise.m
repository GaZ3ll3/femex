function w = external_denoise(fem, v)

% fem = FEM([-1 -1 1 -1 0 sqrt(3) - 1]', 1, 1/(400), []');

    n = size(fem.Promoted.nodes(1,:),2);
    nelem  = size(fem.Promoted.elems, 2);

    I = zeros(3 * nelem,1);
    J = zeros(3 * nelem,1);
    V = ones(3 * nelem,1);
    for i = 1:nelem
        I(3 * i -2) = fem.Promoted.elems(1, i);
        J(3 * i -2) = fem.Promoted.elems(2, i); 

        I(3 * i -1) = fem.Promoted.elems(2, i);
        J(3 * i -1) = fem.Promoted.elems(3, i); 

        I(3 * i ) = fem.Promoted.elems(3, i);
        J(3 * i ) = fem.Promoted.elems(1, i); 

    end

    adj = sparse(I,J,V);

    for i = 1:n
        s = sum(adj(i,:));
        adj(i,:) = adj(i,:)/s;
    end
    
    for k = 1:size(fem.Promoted.edges, 2)
        e = fem.Promoted.edges(1,k);
        adj(e, :) = 0;
    %     adj(:, e) = 0;
        adj(e,e) = 1;
    end
    

    w = adj * v;

end



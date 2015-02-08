function [T, Y] = hsolver(hobj,U0, Ut0, n)


length = size(hobj.fem.Promoted.nodes, 2);

u0 = U0;
v0 = Ut0;
w0 = zeros(length, 1);
z0 = zeros(length, 1);

% solve dy/dt = f(y), using RK4
[T, Y] = ode45(@rigid, [0, 1], [u0;v0;w0;z0]); 


    function dy = rigid(t, y)
        
        dy = zeros(4 * length, 1);
        
        u = y(1:length);
        v = y(length + 1: 2 * length);
        w = y(2 * length + 1: 3*length);
        z = y(3 * length + 1: 4*length);
        
        [u1, v1, w1, z1] = f(t, u, v, w, z);
        
        dy(1:length) = u1;
        dy(length + 1: 2*length) = v1;
        dy(2 * length + 1: 3*length) = w1;
        dy(3 * length + 1: 4*length) = z1;
        
    end





    function [u, v, w ,z] = f(t, U, V, W, Z) 
        
        u = V;
        v = -(hobj.Sm + hobj.Mxy)*U - (hobj.Mx + hobj.My)*V + ...
            hobj.P*W + hobj.Q*Z + g(t);
        w = -(hobj.Mx)*W + hobj.Px*U;
        z = -(hobj.My)*Z + hobj.Qy*U;
        
        tmp = (hobj.M)\[v ,w, z];
        v = tmp(:, 1);
        w = tmp(:, 2);
        z = tmp(:, 3);
    end

    function theta = g(t)
        
        f0 = 15;
        theta = zeros(length, 1);
        theta(n) = ...
            -exp(-pi^2 * (f0 * t - 1)^2) * ...
            (pi*pi)* 2*f0*(f0 * t- 1);
        
    end



end


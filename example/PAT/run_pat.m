h = pat_init(1, 1/20000);
l = size(h.fem.Promoted.nodes,2);
u0 = zeros(l , 1);
ut0 = zeros(l, 1);


dt = 0.005;
tl = 5;

for i = 1:l
    
    x = h.fem.Promoted.nodes(1, i);
    y = h.fem.Promoted.nodes(2, i);
    
    r = sqrt(x^2 + y^2);
    if r <= 1/4
        u0(i) = cos(4*pi * r) + 1.0 ;
    end
    
end
tic;
[T, Y] = pat_solver(h, u0, ut0, dt, tl);
toc;

tic;
Z = time_reversal(h, Y, dt, tl);
toc;
h = pml_init(1, 1/20000);
l = size(h.fem.Promoted.nodes,2);
u0 = zeros(l , 1);
ut0 = zeros(l, 1);
m = 10; n = 1; 
for i = 1:l
    if (abs(h.fem.Promoted.nodes(1,i)) + abs(h.fem.Promoted.nodes(2,i) ) < m)
    m = abs(h.fem.Promoted.nodes(1,i)) + abs(h.fem.Promoted.nodes(2,i));
    n = i;
    end
end
tic;
[T, Y] = pml_solver(h, u0, ut0, n);
toc;
% pml_show(h, Y)
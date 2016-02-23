h = pat_init(1, 1/20000);
l = size(h.fem.Promoted.nodes,2);
u0 = zeros(l , 1);
ut0 = zeros(l, 1);
m = 10; n = 1; 
for i = 1:l
    if (h.fem.Promoted.nodes(1,i)<=0.5) && (h.fem.Promoted.nodes(1,i)>=-0.5)...
            && (h.fem.Promoted.nodes(2,i)<=0.5) && (h.fem.Promoted.nodes(2,i)>=-0.5)
        u0(i) = sin(2 * pi * h.fem.Promoted.nodes(1,i)) * sin(2 * pi * h.fem.Promoted.nodes(2,i));
    end
end


tic;
[T, Y] = pat_solver(h, u0, ut0);
toc;

pat_show(h,Y);
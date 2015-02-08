function show(h, Y)

timestep = size(Y, 1);

for i = 1: timestep/10
trimesh(h.fem.TriMesh', h.fem.Promoted.nodes(1,1:h.fem.Num_nodes), h.fem.Promoted.nodes(2, 1:h.fem.Num_nodes), Y(i*10, 1:h.fem.Num_nodes));
view(2); colorbar;
caxis manual; 
caxis([-1 1]);
str = sprintf('frame at %d',i);
title(str);
pause(0.00001);
end



end


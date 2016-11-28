function pml_show(h, Y)

% filename = 'pml.gif';

figure(1);

timestep = size(Y, 1);

for i = 1: timestep
trimesh(h.fem.TriMesh', h.fem.Promoted.nodes(1,1:h.fem.Num_nodes), h.fem.Promoted.nodes(2, 1:h.fem.Num_nodes), Y(i, 1:h.fem.Num_nodes));
view(2); colorbar;
caxis manual; 
caxis([-1, 1]);
str = sprintf('frame at %d',i);
title(str);
drawnow;

% frame = getframe(1);
% im = frame2im(frame);
% [imind, cm] = rgb2ind(im, 256);
% 
% if i == 1
%    imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
% else
%    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
% end
end



end


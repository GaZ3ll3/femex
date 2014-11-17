function  fem = buildmesh(prec, min_area)
%buildmesh build mesh for finite element.

% It is a unit square.
fem = FEM([0 0 1 0 1 1 0 1]', prec, min_area);

end



function [interpolate] = mapping(func, elems, trans_ref)
    % interpolation mapping

    % allocate memory
    numberofqnodes = size(trans_ref, 1);
    interpolate = zeros(numberofqnodes, size(elems, 2));
    for i = 1: size(elems, 2)
        interpolate(:, i) = trans_ref * func(elems(:, i));
    end
end


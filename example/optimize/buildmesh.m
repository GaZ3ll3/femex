function  fem = buildmesh(prec, min_area)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('~/Documents/github/femex'));

fem = FEM([0 0 1 0 1 1 0 1]', prec, min_area);

end



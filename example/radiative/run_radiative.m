global p

p = radiative();

p.LBFGS('obj', 'grad', 'callback_radiative');

p.plot();

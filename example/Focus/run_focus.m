function run_focus()

[h, t] = focus_init(1, 1/5000, 0.2, 50);

constrains = size(h.ndofs, 1);
freedom = size(h.fem.Promoted.elems, 2);

numberofeqn = ceil(0.25 * freedom/constrains);

fprintf('%d measurements needed.\n', numberofeqn);


% in -.5 ~ .5
incidents = rand(2 , numberofeqn) - 0.5;

lhs = zeros(numberofeqn * constrains, freedom);
rhs = zeros(numberofeqn * constrains, 1);


for i = 1:numberofeqn
    
    e = focus_solve(h, incidents(:, i));
    [lhs((i - 1) * constrains + 1 : i * constrains, :), rhs( (i - 1) * constrains + 1: i * constrains )] = focus_diffusive(h, e);
    
end

cvx_quiet(true);
cvx_begin
    variable x(size(lhs, 2))
    minimize( norm( x, 2 ) )
    subject to
    lhs * x == rhs
    x >= 0
    x <= 1
cvx_end


u = zeros(size(x ,1), 1);
u(t) = 1.0;

fprintf('Error with target is %6.8f.\n',norm(x - u) );



% res = lhs \ rhs;
% 
% disp(norm(res - h.rho));


end


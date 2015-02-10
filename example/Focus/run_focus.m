function [lhs, rhs] = run_focus()

h = focus_init(1, 1/3200, 0.2, 50);

constrains = size(h.ndofs, 1);
freedom = size(h.fem.Promoted.elems, 2);


numberofeqn = 3 * ceil(freedom/constrains);

fprintf('%d measurements needed.\n', numberofeqn);


% in -.5 ~ .5
incidents = rand(2 , numberofeqn) - 0.5;

lhs = zeros(numberofeqn * constrains, freedom);
rhs = zeros(numberofeqn * constrains, 1);


for i = 1:numberofeqn
    
    e = focus_solve(h, incidents(:, i));
    [lhs((i - 1) * constrains + 1 : i * constrains, :), rhs( (i - 1) * constrains + 1: i * constrains )] = focus_diffusive(h, e);
    
end


% res = lhs \ rhs;
% 
% disp(norm(res - h.rho));


end


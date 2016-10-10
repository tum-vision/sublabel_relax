function [result] = solve_direct_rof_nd(f, lambda)
%SOLVE_DIRECT_ROF_ND Implementation of the classical ROF model with nuclear norm TV
%
% min_u (1/2) (u-f)^2 + \lambda |\nabla u|_*
%

    [ny, nx, nc] = size(f);

    u = prost.variable(nx*ny*nc);
    p = prost.variable(2*nx*ny*nc);

    problem = prost.min_problem({u}, {p});

    % u = nabla p
    problem.add_constraint(u, p, prost.block.gradient2d(nx, ny, nc, false));

    % 0.5*|u-f|^2
    problem.add_function(u, prost.function.sum_1d('square', 1, f(:), ...
                                      1, 0, 0));

    % lambda * |p|_*
    problem.add_function(p, prost.function.sum_singular_nx2(2*nc, false, 'sum_1d:abs', ...
                                                lambda, 0, 1, 0, 0));
    %% create backend and specify solver options
    backend = prost.backend.pdhg('tau0', 1, 'sigma0', 1);
    opts = prost.options('tol_abs_dual', 1e-6, ...
                         'tol_abs_primal', 1e-6, ...
                         'tol_rel_dual', 1e-6, ...
                         'tol_rel_primal', 1e-6, ...
                         'max_iters', 5500);


    %% solve problem
    prost.solve(problem, backend, opts);

    result = reshape(u.val, ny, nx, nc);
end

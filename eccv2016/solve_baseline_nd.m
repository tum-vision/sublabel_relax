function [u_proj, u_lifted] = solve_baseline_nd(vert, tri, ny, nx, lambda, data)
%SOLVE_BASELINE_ND Implementation of the standard dataterm with the
%sublabel accurate regularizer from Lellmann'13
%   vert, tri: Triangulation of the range
%   ny, nx: height width of the domain
%   lambda: regularization parameter
%   data: struct, data.f dataterm as a cost volume

    %%
    % Linear operator for data term
    [L, n] = size(vert);
    [T, ~] = size(tri);
    N = ny*nx;

    W_cell = cell(n, 1);
    for d=1:n
        W_cell{d, 1} = sparse(0, L);
    end

    for i=1:T
        M = sparse(n, 0);
        Q = sparse(0, L);
        for d=1:n
            M = cat(2, M, vert(tri(i, d), :)' - vert(tri(i, n+1), :)');
            Q = cat(1, Q, double((1:L == tri(i, d)) - (1:L == tri(i, n+1))));
        end
        M = inv(M);

        for d=1:n
            W_cell{d, 1} = cat(1, W_cell{d, 1}, M(:, d)' * Q);
        end
    end

    W_stacked = sparse(0, L);
    for d=1:n
        W_stacked = cat(1, W_stacked,  W_cell{d, 1});
    end
    W = kron(speye(2), W_stacked);
    
    u = prost.variable(N*L);
    w_tilde = prost.variable(n*2*N*T);
    
    p_tilde = prost.variable(2*N*L);
    r_tilde = prost.variable(n*2*N*T);
    

    problem = prost.min_max_problem({u, w_tilde}, ...
                                    {p_tilde, r_tilde});
    
    % \delta_Simplex(u) + <u, f>
    problem.add_function(u, prost.function.transform(prost.function.sum_ind_simplex(L, false), ...
                             1, 0, 1, data.f(:), 0));
                         
    % add functions and pairings for the regularizer
    problem.add_dual_pair(u, p_tilde, prost.block.gradient2d(nx, ny, L, false));
    problem.add_dual_pair(w_tilde, p_tilde, prost.block.sparse_kron_id(W', N));
    problem.add_dual_pair(w_tilde, r_tilde, prost.block.identity());
    problem.add_function(r_tilde, prost.function.conjugate(prost.function.sum_singular_nx2(n*2, false, 'sum_1d:abs', ...          
                                                      lambda, 0, 1, 0, 0)));

    % create backend
    backend = prost.backend.pdhg('residual_iter', 100);

    % specify solver options
    opts = prost.options('tol_abs_dual', 1e-6, ...
                         'tol_abs_primal', 1e-6, ...
                         'tol_rel_dual', 1e-6, ...
                         'tol_rel_primal', 1e-6);
    opts.max_iters = 50600;
    opts.num_cback_calls = 20;

    % solve problem
    prost.solve(problem, backend, opts);

    % project solution of lifted problem to the original domain by
    %   computing a weighted sum where at each pixel the entries of the vector
    %   u_lifted are weighted with the corresponding vertex.
    u_proj = reshape(reshape(u.val, N, L) * vert, ny, nx, n);
    u_lifted = u.val;
end



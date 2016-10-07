function [u_proj, u_lifted] = solve_sublabel_nd(vert, tri, ny, nx, mode, lambda, data)
%LINOP_SUBLABEL_ND Summary of this function goes here
%   Detailed explanation goes here

    %%
    % Linear operator for data term
    [L, n] = size(vert);
    [T, ~] = size(tri);
    N = ny*nx;
    
    AtQ_cell = cell(n, 1);
    for d=1:n
        AtQ_cell{d, 1} = sparse(0, L);
    end
    btQ = sparse(0, L);


    for i=1:T
        M = sparse(n, 0);
        Q = sparse(0, L);
        for d=1:n+1
            Q = cat(1, Q, double(1:L == tri(i, d)));
            M = cat(2, M, vert(tri(i, d), :)');
        end
        M = inv(cat(1, M, ones(1, n+1)));
        A = M(:, 1:n);
        b = M(:, n+1);

        AtQ = A'*Q;

        for d=1:n
            AtQ_cell{d, 1} = cat(1, AtQ_cell{d, 1}, AtQ(d, :));
        end

        btQ = cat(1, btQ, b'*Q);
    end

    AtQ = sparse(0, L);
    for d=1:n
        AtQ = cat(1, AtQ, AtQ_cell{d, 1});
    end
    
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
    m = prost.variable(n*N*T);
    o = prost.variable(N*T);
    
    w_tilde = prost.variable(n*2*N*T);

    v = prost.variable(N*L);
    q = prost.variable(N);
    rs = prost.variable((n+1)*N*T);
    r = prost.sub_variable(rs, n*N*T);
    s = prost.sub_variable(rs, N*T);
    
    p_tilde = prost.variable(2*N*L);
    r_tilde = prost.variable(n*2*N*T);
    
    if(strcmp(mode, 'quad') || strcmp(mode, 'quad_trunc'))
        X = sparse(0, T * n);
        for d=1:n+1
            B = sparse(T, 0);

            for k=1:n
                B = cat(2, B, vert(tri(:, d), k));
            end

            X = cat(1, X, spdiags(B, (0:n-1)*T, T, n * T));
        end
        Y = repmat(speye(T), n+1, 1);
        
        x = prost.variable(n*N*T);
        y = prost.variable(N*T);
        w = prost.variable((n+1)*N*T);
        
        ab = prost.variable((n+1)*N*T);
        a = prost.sub_variable(ab, n*N*T);
        b = prost.sub_variable(ab, N*T);
        c = prost.variable(n*N*T);
        d = prost.variable(N*T);
    end
    
    if(strcmp(mode, 'quad_trunc'))
        t = prost.variable((n+1)*N*T);
        
        problem = prost.min_max_problem({u, m, o, x, y, w, t, w_tilde}, ...
                                        {v, q, rs, ab, c, d, p_tilde, r_tilde});
    elseif(strcmp(mode, 'quad'))
        problem = prost.min_max_problem({u, m, o, x, y, w, w_tilde}, ...
                                        {v, q, rs, ab, c, d, p_tilde, r_tilde});
    elseif(strcmp(mode, 'polyhed'))
        problem = prost.min_max_problem({u, m, o, w_tilde}, ...
                                        {v, q, rs, p_tilde, r_tilde});
    else
        error('Warning: Unrecognized mode.');
    end
    

    problem.add_dual_pair(u, v, prost.block.identity());

    % r = A^T * Q * v
    problem.add_dual_pair(m, v, prost.block.sparse_kron_id(AtQ', N));
    problem.add_dual_pair(m, r, prost.block.identity(-1));

    % s = q - b^T * Q * v
    problem.add_dual_pair(o, v, prost.block.sparse_kron_id(-btQ', N));
    Z = sparse(ones(T, 1));
    problem.add_dual_pair(o, q, prost.block.sparse_kron_id(Z', N));
    problem.add_dual_pair(o, s, prost.block.identity(-1));

    if(strcmp(mode, 'quad') || strcmp(mode, 'quad_trunc'))
        % a + c = r
        problem.add_dual_pair(x, r, prost.block.identity(-1));
        problem.add_dual_pair(x, a, prost.block.identity());
        problem.add_dual_pair(x, c, prost.block.identity());

        % b + d = s
        problem.add_dual_pair(y, s, prost.block.identity(-1));
        problem.add_dual_pair(y, b, prost.block.identity());
        problem.add_dual_pair(y, d, prost.block.identity());

        % <c, t^{i_j}> <= d for all vertices t^{i_j} for all triangles Delta_i
        problem.add_dual_pair(w, c, prost.block.sparse_kron_id(X', N));
        problem.add_dual_pair(w, d, prost.block.sparse_kron_id(-Y', N));
        
        problem.add_function(w, prost.function.sum_1d('ind_leq0', 1,0,1,0,0));
        % coeffs for prox_{delta_{epi(h^*)}} with 
        % h^*(y) = 0.5 * <y,y> + <f,y> and 
        % h(x) = 0.5 * || x - f ||_2^2
        coeffs_a=0.5;
        
        coeffs_b = repmat(reshape(data.f(:, :, 1), N, 1), T, 1);
        for i=2:n
            coeffs_b=cat(1, coeffs_b, repmat(reshape(data.f(:, :, i), N, 1), T, 1));
        end
        coeffs_c=0;

        problem.add_function(ab, prost.function.sum_ind_epi_quad(n+1, false, coeffs_a, coeffs_b, coeffs_c));

    elseif(strcmp(mode, 'polyhed'))
        problem.add_function(rs, sum_ind_epi_polyhedral(n+1, false, data.pts_x, data.pts_y, ...
                                                 data.counts, data.index));
    end
    
    if(strcmp(mode, 'quad_trunc'))
        % <r, t^{i_j}> <= s for all vertices t^{i_j} for all triangles Delta_i
        problem.add_dual_pair(t, r, prost.block.sparse_kron_id(X', N));
        problem.add_dual_pair(t, s, prost.block.sparse_kron_id(-Y', N));
        
        problem.add_function(t, prost.function.sum_1d('ind_leq0', 1,0,1,-data.nu,0));
    end

    problem.add_dual_pair(u, p_tilde, prost.block.gradient2d(nx, ny, L, false));
    problem.add_dual_pair(w_tilde, p_tilde, prost.block.sparse_kron_id(W', N));
    problem.add_dual_pair(w_tilde, r_tilde, prost.block.identity());

    
    problem.add_function(q, prost.function.sum_1d('zero', 1,0,1,1,0));
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


function [u_unlifted] = baseline_lifting(cost_volume, gamma, lmb)
% this is an implementation of the standard baseline lifting method
% proposed by Pock et. al in SIIMS '10. 
    
    [ny, nx, L] = size(cost_volume);
    lmb_scaled = lmb * (gamma(2) - gamma(1));
    f = cost_volume(:);
    N = ny*nx;

    u = prost.variable(N*L); 
    q = prost.variable(N*L*3);
    
    % we need sub variables, since we have different prox operators
    % on indivudal parts of the variables.
    u1 = prost.sub_variable(u, N);
    u2 = prost.sub_variable(u, N*(L-1));
    q1 = prost.sub_variable(q, N*L*2);
    q2 = prost.sub_variable(q, N*L);

    problem = prost.min_max_problem( {u}, {q} );

    % force lowest layer of volume to 1
    problem.add_function(u1, prost.function.sum_1d('ind_eq0', 1, 1, ...
                                                   1, 0, 0));
        
    % box constraint on rest
    problem.add_function(u2, prost.function.sum_1d('ind_box01', 1, 0, 1, 0, 0));
    
    % total variation |phi_x| <= lmb
    problem.add_function(q1, prost.function.sum_norm2(2, false, 'ind_leq0', 1 / lmb_scaled, ...
                                      1, 1, 0, 0));
    
    % dataterm, phi_t >= -f
    problem.add_function(q2, prost.function.sum_1d('ind_geq0', 1, -f, 1, 0, 0));

    % linear operator is just <grad u, q> with 3d gradient
    problem.add_dual_pair(u, q, prost.block.gradient3d(nx, ny, L));
        
    %% create backend
    backend = prost.backend.pdhg(...
        'tau0', 10, ...
        'sigma0', 0.1);

    %% specify solver options
    opts = prost.options(...
        'max_iters', 5000, ...
        'num_cback_calls', 25);
    
    %% solve problem
    solution = prost.solve(problem, backend, opts);

    %% retrieve unlifted result via layer-cake formula
    u_volume = reshape(u.val, ny, nx, L);
    u_unlifted = (sum(u_volume, 3) - 1) / (L - 1);
    
end

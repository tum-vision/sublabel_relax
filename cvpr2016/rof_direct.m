%% ROF with direct convex optimization

%% read input image
im = imread('images/24004.jpg');
im = double(imresize(im, 0.5)) / 255;
[ny, nx,nc] = size(im);
lmb = 0.5;
f = im(:);

%% create problem description
u = prost.variable(nx*ny*nc);
q = prost.variable(2*nx*ny*nc);

prob = prost.min_max_problem( {u}, {q} );

prob.add_function(u, prost.function.sum_1d('square', 1, f, 2, 0, 0));
prob.add_function(q, prost.function.sum_norm2(... 
    2 * nc, false, 'ind_leq0', 1 / lmb, 1, 1, 0, 0));

prob.add_dual_pair(u, q, prost.block.gradient2d(nx,ny,nc));


%% create backend
backend = prost.backend.pdhg();

%% specify solver options
opts = prost.options(...
    'max_iters', 20000, ...
    'tol_rel_primal', 1e-7, ...
    'tol_abs_primal', 1e-7, ...
    'tol_rel_dual', 1e-7, ...
    'tol_abs_dual', 1e-7);

%% solve problem
solution = prost.solve(prob, backend, opts);


%% show result
result = solution.x;
Kmat = spmat_gradient2d(nx,ny,1);
[m, n] = size(Kmat);
grad = Kmat * result(:);
gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
en_prim = sum((result(:)-im(:)).^2) + lmb * sum(gradnorms);
fprintf('ROF_primal_energy=%f\n', en_prim);

imshow(reshape(result, ny, nx));

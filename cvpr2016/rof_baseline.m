%% ROF with baseline lifting (see Fig. 5 in the paper)

%% read input image
im = imread('images/24004.jpg');
im = double(imresize(im, 1)) / 255;
[ny, nx] = size(im);

%% setup parameters
L = 10;
gamma = linspace(0, 1, L); 
lmb = 0.5;

cost_volume = zeros(ny, nx, L);
for i=1:L
    cost_volume(:, :, i) = (gamma(i) - im).^2;
end

%% solve problem and display result
result = baseline_lifting(cost_volume, gamma, lmb);

Kmat = spmat_gradient2d(nx,ny,1);
[m, n] = size(Kmat);
grad = Kmat * result(:);
gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
en_prim = sum((result(:)-im(:)).^2) + lmb * sum(gradnorms);
fprintf('ROF_primal_energy=%f\n', en_prim);

imshow(result);

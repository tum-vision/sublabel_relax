%% Truncated quadratic dataterm + TV denoising (see Fig.7) with baseline lifting
rng(42);

% load image
im = imread('images/watercastle.jpg');
im = double(imresize(im, 1)) / 255;
[ny, nx] = size(im);
N = ny * nx;

% noise parameters
noise_sigma = 0.05; % standard deviation of gaussian noise
noise_sp = 0.25;    % percentage of salt&pepper noise

% add gaussian noise
im_noisy = im + noise_sigma * randn(ny, nx, 1);

% add salt and pepper noise
im_noisy = im_noisy(:);
perm = randperm(N);
num_sp = round(N *noise_sp * 0.5);
im_noisy(perm(1:num_sp)) = 1;
im_noisy(perm(num_sp+1:2*num_sp)) = 0;
im_noisy = min(max(im_noisy, 0), 1);
im_noisy = reshape(im_noisy, [ny, nx]);

% parameters
L = 5; % number of labels
gamma = linspace(0, 1, L)';
nu = 0.025; % a jump higher than sqrt(0.025) is considered an outlier
alpha = 25;
lmb = 1;

% create cost volume
cost_volume = zeros(ny, nx, L);
for i=1:L
    cost_volume(:, :, i) = (alpha / 2) * min((gamma(i) - im_noisy).^2, nu);
end

%% solve problem and display result
u_unlifted = baseline_lifting(cost_volume, gamma, lmb);

%% compute truncated ROF energy
Kmat = spmat_gradient2d(nx,ny,1);
[m, n] = size(Kmat);
grad = Kmat * u_unlifted(:);
gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
en_prim = 0.5 * alpha * sum(min((u_unlifted(:)-im_noisy(:)).^2, nu)) ...
          + lmb * sum(gradnorms);
fprintf('TruncROF_energy=%f\n', en_prim);

imshow(reshape(u_unlifted, [ny, nx]));

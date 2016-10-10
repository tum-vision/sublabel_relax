F_MAX = 1;

f = imread('images/210088.jpg');
f = imresize(f, 1);
f = im2double(f);


%% Truncated quadratic dataterm + TV denoising (see Fig.7) with sublabel lifting
rng(42);

[ny, nx, ~] = size(f);
N = nx*ny;

% noise parameters
noise_sigma = 0.05; % standard deviation of gaussian noise
noise_sp = 0.25;    % percentage of salt&pepper noise

% add gaussian noise
f_noisy = f + noise_sigma * randn(ny, nx, 3);

% add salt and pepper noise
f_noisy = f_noisy(:);
perm = randperm(N*3);
num_sp = round(N*3 *noise_sp * 0.5);
f_noisy(perm(1:num_sp)) = 1;
f_noisy(perm(num_sp+1:2*num_sp)) = 0;
f_noisy = min(max(f_noisy, 0), 1);
f_noisy = reshape(f_noisy, [ny, nx, 3]);

% parameters
nu = 0.025; % a jump higher than sqrt(0.025) is considered an outlier
lambda = 0.03;
n = 3;

% L=4;
% T=1;
% vert = [0 0 0; 3 0 0; 0 0 3; 0 3 0];
% tri = [1 2 3 4];

l = 4;
t=linspace(0, F_MAX, l);
[vert, tri] = triang3d_box(t);
L = size(vert, 1);
T = size(tri, 1);


dataterm = zeros(ny, nx, l, l, l);

for i=1:ny
    for j=1:nx
        dataterm(i, j, :) = min(nu, 0.5*sum((repmat(squeeze(f_noisy(i,j,:))', L, 1)-vert).^2, 2));
    end
end

data.f = dataterm;
[u_proj, ~] = solve_baseline_nd(vert, tri, ny, nx, lambda, data);

%%
% show result and compute energy
figure;
imshow(u_proj, []);

E = sum(min(nu, reshape(0.5 * sum((u_proj - f_noisy).^2, 3), N, 1)));

Nabla = spmat_gradient2d(ny, nx, 3);
Grad = reshape(Nabla * u_proj(:), ny, nx, 3, 2);
for i=1:ny
    for j=1:nx
        T = squeeze(Grad(i, j, :, :));
        E = E + lambda * sum(svd(T));
    end
end

['    energy (unlifted): ', num2str(E)]

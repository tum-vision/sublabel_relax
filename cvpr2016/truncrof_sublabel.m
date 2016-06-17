%% Truncated quadratic dataterm + TV denoising (see Fig.7) with sublabel lifting
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

if (gamma(2) - gamma(1)) > 2 * sqrt(nu)
    error('ERROR: Not enough labels for dataterm!')
end

% compute piecewise quadratic approximation of dataterm
polya = zeros(L - 1, ny, nx);
polyb = zeros(L - 1, ny, nx);
polyc = zeros(L - 1, ny, nx);

for i=1:(L-1)
    % completely piecewise constant part
    case1 = (gamma(i + 1) <= (im_noisy - sqrt(nu))) | ...
            (gamma(i) >= (im_noisy + sqrt(nu)));
    
    % completely quadratic part
    case2 = (gamma(i) >= (im_noisy - sqrt(nu))) & ...
            (gamma(i + 1) <= (im_noisy + sqrt(nu)));
    
    % on left kink
    case3 = (gamma(i) < (im_noisy - sqrt(nu))) & ...
            (gamma(i+1) > (im_noisy - sqrt(nu)));
    
    % on right kink
    case4 = (gamma(i) < (im_noisy + sqrt(nu))) & ...
            (gamma(i+1) > (im_noisy + sqrt(nu)));
    
    polya(i, case1) = 0;
    polyb(i, case1) = 0;
    polyc(i, case1) = alpha * nu / 2;
    
    polya(i, case2) = alpha / 2;
    polyb(i, case2) = -alpha * im_noisy(case2);
    polyc(i, case2) = (alpha / 2) * im_noisy(case2) .^ 2;

    polya(i, case3) = 0;
    polyb(i, case3) = ((alpha / 2) * (gamma(i + 1) - im_noisy(case3)) .^2 - (alpha * nu) / 2) / ...
        (gamma(i + 1) - gamma(i));
    polyc(i, case3) = (alpha * nu / 2) - polyb(i, case3) * gamma(i);

    polya(i, case4) = 0;
    polyb(i, case4) = ((alpha / 2) * (gamma(i) - im_noisy(case4)) .^2 - (alpha * nu) / 2) / ...
        (gamma(i) - gamma(i+1));
    polyc(i, case4) = (alpha * nu / 2) - polyb(i, case4) * gamma(i + 1);
end


%% solve problem
u_unlifted = sublabel_lifting_quad(polya, polyb, polyc, gamma, lmb);

%% compute truncated ROF energy
Kmat = spmat_gradient2d(nx,ny,1);
[m, n] = size(Kmat);
grad = Kmat * u_unlifted(:);
gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
en_prim = 0.5 * alpha * sum(min((u_unlifted(:)-im_noisy(:)).^2, nu)) ...
          + lmb * sum(gradnorms);
fprintf('TruncROF_energy=%f\n', en_prim);

imshow(reshape(u_unlifted, [ny, nx]));

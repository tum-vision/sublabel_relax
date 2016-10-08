%% ROF with sublabel lifting (see Fig. 5)

%% read input image
im = imread('images/24004.jpg');
im = double(imresize(im, 0.1)) / 255;
[ny, nx] = size(im);

%% setup parameters
L = 4;
gamma = linspace(0, 1, L)'; 
lmb = 0.1;
N = ny * nx;

%% compute piecewise quadratic approximation of dataterm
polya = zeros(L - 1, ny, nx);
polyb = zeros(L - 1, ny, nx);
polyc = zeros(L - 1, ny, nx);

for i=1:(L-1)
    polya(i,:, :)=1;
    polyb(i,:, :)=-2*im;
    polyc(i,:, :)=(im.^2);
end

%% solve problem
u_unlifted = sublabel_lifting_quad(polya, polyb, polyc, gamma, lmb);

%% compute ROF energy
Kmat = spmat_gradient2d(nx,ny,1);
[m, n] = size(Kmat);
grad = Kmat * u_unlifted(:);
gradnorms = sqrt(grad(1:n).^2 + grad(n+1:end).^2);
en_prim = sum((u_unlifted(:)-im(:)).^2) + lmb * sum(gradnorms);
fprintf('ROF_primal_energy=%f\n', en_prim);

imshow(reshape(u_unlifted, [ny, nx]));

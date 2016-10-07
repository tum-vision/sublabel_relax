f = imread('images/172032.jpg');
f = imresize(f, 0.2);
f = im2double(f);
[ny, nx, ~] = size(f);
lambda = 0.5;

%% solve unlifted problem
result = unlifted_rof(f, lambda);
figure;
imshow(result, []);
['energy (unlifted): ', num2str(energy_rof(result, f, lambda))]

l = 2;
t=linspace(0, 1, l);
[vert, tri] = triang3d_box(t);
L = size(vert, 1);
T = size(tri, 1);
data.f = f;
[u_proj, u_lifted] = solve_sublabel_nd(vert, tri, ny, nx, 'quad', lambda, data);

figure;
imshow(u_proj, []);
['energy (sublabel): ', num2str(energy_rof(u_proj, f, lambda))]



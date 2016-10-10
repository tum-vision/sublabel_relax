f = imread('images/172032.jpg');
f = imresize(f, 0.2);
f = im2double(f);
[ny, nx, ~] = size(f);
lambda = 0.5;

l = 2;
t = linspace(0, 1, l);
[vert, tri] = triang3d_box(t);
L = size(vert, 1);
T = size(tri, 1);
%% solve problem via direct optimization
result = solve_direct_rof_nd(f, lambda);
figure;
imshow(result, []);
fprintf('\n');
fprintf('energy (direct optimization): %f\n\n', energy_rof(result, f, lambda));

%%
dataterm = zeros(ny, nx, l, l, l);
for i=1:ny
    for j=1:nx
        dataterm(i, j, :) = 0.5*sum((repmat(squeeze(f(i,j,:))', L, 1)-vert).^2, 2);
    end
end

data.f = dataterm;

[u_proj, ~] = solve_baseline_nd(vert, tri, ny, nx, lambda, data);

result_baseline = u_proj;
imshow(result_baseline, []);
fprintf('\n');
fprintf('energy (baseline): %f\n\n', energy_rof(result_baseline, f, lambda));

%% solve problem via sublabel accurate lifting
data.f = f;
[u_proj, u_lifted] = solve_sublabel_nd(vert, tri, ny, nx, 'quad', lambda, data);

figure;
imshow(u_proj, []);
fprintf('\n');
fprintf('energy (direct optimization): %f\n\n', energy_rof(u_proj, f, lambda));



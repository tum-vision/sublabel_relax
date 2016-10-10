I1 = double(imresize(imread('images/circle_frame0.png'), 1));
I2 = double(imresize(imread('images/circle_frame1.png'), 1));

%I1 = double(imresize(imread('../../data/tight-convex-figs/fig11_a.png'), 1));
%I2 = double(imresize(imread('../../data/tight-convex-figs/fig11_b.png'), 1));

%I1 = double(imresize(imread('eval-data/Backyard/frame07.png'), 1));
%I2 = double(imresize(imread('eval-data/Backyard/frame14.png'), 1));

%I1 = double(imresize(imread('eval-data/Evergreen/frame10.png'), 0.5));
%I2 = double(imresize(imread('eval-data/Evergreen/frame11.png'), 0.5));

%I1 = double(imresize(imread('other-data/Urban2/frame10.png'), 1));
%I2 = double(imresize(imread('other-data/Urban2/frame11.png'), 1));

%I1 = double(imresize(imread('other-data/Beanbags/frame10.png'), 1));
%I2 = double(imresize(imread('other-data/Beanbags/frame11.png'), 1));

% lambda = 0.5, 29x29, 75 sublabel
%I1 = double(imresize(imread('other-data/Grove3/frame10.png'), 1));
%I2 = double(imresize(imread('other-data/Grove3/frame11.png'), 1));

[ny, nx, ~] = size(I1);
N= nx*ny;

l = 15;
L = l*l;
T = (l-1)*(l-1)*2;

width = 7;

lambda = 0.5;

% last parameter set to true => outputs cost volume instead
% of convexifying the energy on each triangle.
[vert, tri, cost] = ...
    compute_cost_flow_piecw_conv(I1, I2, width, l, 0, true);
cost = cost ./ max(cost(:));
vert = vert ./ (width-1);

tic;

cost = permute(cost, [1 2 4 3]);

data.f = cost;

[result, ~] = solve_baseline_nd(vert, tri, ny, nx, lambda, data);
toc;

flow = (result - 0.5) * (width - 1);

figure;
imshow(flowToColor(flow));

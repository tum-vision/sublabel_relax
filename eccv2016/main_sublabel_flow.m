%I1 = double(imresize(imread('images/circle_frame0.png'), 1));
%I2 = double(imresize(imread('images/circle_frame1.png'), 1));

%I1 = double(imresize(imread('images/eval-data/Backyard/frame07.png'), 1));
%I2 = double(imresize(imread('images/eval-data/Backyard/frame14.png'), 1));

%I1 = double(imresize(imread('images/eval-data/Evergreen/frame10.png'), 0.5));
%I2 = double(imresize(imread('images/eval-data/Evergreen/frame11.png'), 0.5));

%I1 = double(imresize(imread('images/other-data/Urban2/frame10.png'), 1));
%I2 = double(imresize(imread('images/other-data/Urban2/frame11.png'), 1));

%I1 = double(imresize(imread('images/other-data/Beanbags/frame10.png'), 1));
%I2 = double(imresize(imread('images/other-data/Beanbags/frame11.png'), 1));

% lambda = 0.5, 29x29, 75 sublabel
I1 = double(imresize(imread('images/other-data/Grove3/frame10.png'), 1));
I2 = double(imresize(imread('images/other-data/Grove3/frame11.png'), 1));

[ny, nx, ~] = size(I1);
N= nx*ny;


%% solve problem
l = 5;
L = l*l;
T = (l-1)*(l-1)*2;

width = 29;
total_sublabels = 75;

[vert, tri, index, counts, pts_x, pts_y] = ...
    compute_cost_flow_piecw_conv(I1, I2, width, l, floor(total_sublabels / (l-1)));

vert = vert / (width-1);
pts_x = pts_x / (width-1);
pts_y = pts_y / (max(pts_y));

lambda = 0.5;
n=2;

data.pts_x = pts_x;
data.pts_y = pts_y;
data.counts = double(counts);
data.index = double(index);
[u_proj, u_lifted] = solve_sublabel_nd(vert, tri, ny, nx, 'polyhed', lambda, data);

ucounts = squeeze(sum(reshape(u_lifted, [ny, nx, L]) > 0.001, 3));
figure; imagesc(ucounts); colorbar;


flow = (u_proj - 0.5) * (width - 1);

figure;
imshow(flowToColor(flow));
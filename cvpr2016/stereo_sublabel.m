%% Stereo matching with sublabel accurate lifting (see Fig. 1 and 6 in the paper)

%% compute matching cost
ndisps = 135;

if ~exist('stereo_cost')
    stereo_cost = compute_stereo_cost(...
        'images/motorcycle_im0.png', ...
        'images/motorcycle_im1.png', ndisps, 2, 2);
end

%% setup parameters
L = 4;
gamma = linspace(0, 1, L)'; 
lmb = 0.15;

%% solve problem and display result
tic;
[u_unlifted] = sublabel_lifting_convex(stereo_cost, gamma, lmb);
toc;

[ny, nx, ~] = size(stereo_cost);  
u_unlifted = min(max(u_unlifted, 0), 1);
num_colors = 65536;
cmap = jet(num_colors);
cmap_index = 1 + round(u_unlifted * (num_colors - 1));
image_rgb = reshape(cmap(cmap_index,:),ny,nx,3);
imshow(image_rgb);

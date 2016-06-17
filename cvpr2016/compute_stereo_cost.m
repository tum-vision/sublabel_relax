function [cost_volume] = compute_stereo_cost(path_im0, path_im1, ndisps, dscl1, dscl2)
% computes a stereo matching cost between two images given by
% path_im0 and path_im1. ndisps denotes the number of disparities,
% dscl1 is a factor by which to downscale the original image, and
% dscl2 is a factor by which to downscale the dataterm. dscl2 is
% essentially a "patch-size" for the dataterm computation.
    
    % read stereo images
    im_0 = imresize(double(imread(path_im0)) / 255, (1/dscl1));
    im_1 = imresize(double(imread(path_im1)) / 255, (1/dscl1));
    [ny, nx, nc] = size(im_0);

    % compute derivatives
    im_0_dx = [im_0(:, 2:end, :) im_0(:, end, :)] - im_0;
    im_0_dy = [im_0(2:end, :, :); im_0(end, :, :)] - im_0;
    im_1_dx = [im_1(:, 2:end, :) im_1(:, end, :)] - im_1;
    im_1_dy = [im_1(2:end, :, :); im_1(end, :, :)] - im_1;
    
    dt_alpha = 0.1;
    dt_beta = 0.1;

    % compute full dataterm
    idx_L = 1;
    f = zeros(ny, nx, ndisps);
    for d=1:ndisps
        sum_dx = sum(abs([repmat(im_1_dx(:, 1, :), 1, d-1) im_1_dx(:, 1:(nx-d+1), :)] - ...
                         im_0_dx), 3);
        sum_dy = sum(abs([repmat(im_1_dy(:, 1, :), 1, d-1) im_1_dy(:, 1:(nx-d+1), :)] - ...
                         im_0_dy), 3);
        f(:, :, idx_L) = min(dt_alpha, sum_dx) + min(dt_beta, sum_dy);
        idx_L = idx_L + 1;
    end

    % downscale dataterm
    ny2 = floor(ny / dscl2);
    nx2 = floor(nx / dscl2);
    cost_volume = zeros(ny2, nx2, ndisps);

    for i=1:dscl2
        for j=1:dscl2
            cost_volume = cost_volume + f(i:dscl2:end, j:dscl2:end, :);
        end
    end
    cost_volume = cost_volume / (dscl2 * dscl2);

end

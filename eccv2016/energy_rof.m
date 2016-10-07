function [E] = energy_rof(u, f, lambda)
%ENERGY_ROF_3D Summary of this function goes here
%   Detailed explanation goes here
    [ny, nx, nc] = size(f);
    E = 0.5 * sum((u(:) - f(:)).^2);
    Nabla = spmat_gradient2d(ny, nx, nc);
    Grad = reshape(Nabla * u(:), ny, nx, nc, 2);
    for i=1:ny
        for j=1:nx
            T = squeeze(Grad(i, j, :, :));
            E = E + lambda * sum(svd(T));
        end
    end
end


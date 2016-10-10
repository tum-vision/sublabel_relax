function [E] = energy_rof(u, f, lambda)
%ENERGY_ROF Computes the energy of the classical ROF model with nuclear norm TV
% min_u (1/2) (u-f)^2 + \lambda |\nabla u|_* 
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


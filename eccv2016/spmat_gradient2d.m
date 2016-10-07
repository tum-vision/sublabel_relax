function [Nabla] = spmat_gradient2d(ny, nx, nc)
%SPMAT_GRADIENT2D Summary of this function goes here
%   Detailed explanation goes here
    e = ones(ny, 1);
    A = spdiags([-e e], 0:1, ny, ny);
    A(end, end) = 0;
    Dy = kron(speye(nx), A);
    e = ones(nx, 1);
    B = spdiags([-e e], 0:1, nx, nx)';
    B(end, end) = 0;
    Dx = kron(B', speye(ny));
    Nabla = cat(1, kron(speye(nc), Dx), kron(speye(nc), Dy));
end


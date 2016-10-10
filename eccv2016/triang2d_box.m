function [vert, tri] = triang2d_box(t)
%TRIANG2D_BOX Computes a triangulation of a 2d box
%   t is a vector of length l that contains the labels in x, y direction 
%   all elements in t x t correspond to a vertex of the triangulation
    
    l = size(t, 2);
    L = l*l;
    vert = [reshape(repmat(t, l, 1), L, 1), repmat(t, 1, l)']; % vertices, (L^2) * 2 matrix
    T = ((l - 1).^2) * 2; % number of triangles
    tri = zeros(T, 3);

    % make triangle mesh
    k=1;
    for i=1:T/2+l-1
        if(mod(i, l) ~= 0)
            tri(k, :) = [i, i+1, i+l+1];
            k=k+1;
            tri(k, :) = [i+l+1, i+l, i];
            k=k+1;
        end
    end
end


function [vert, tri] = triang3d_box(t)
%TRIANG3D_BOX Computes a triangulation of a 3d box
%   t is a vector of length l that contains the labels in x, y direction 
%   all elements in t x t x t correspond to a vertex of the triangulation
  
    l = size(t, 2);
    L = l*l*l;
    vert = [repmat([reshape(repmat(t, l, 1), l^2, 1), repmat(t, 1, l)'], l, 1), reshape(repmat(t, l^2, 1), L, 1)];

    T = (l-1)^3*6;
    tri = zeros(T, 4);

    a=1;
    for i=1:l-1
        for j=1:l-1
            for k=1:l-1
                b = (i-1)*l*l + (j-1)*l + k;

                x1 = b;
                x2 = b+1;
                x3 = b+l;
                x4 = b+l+1;
                x5 = b+l^2;
                x6 = b+l^2+1;
                x7 = b+l^2+l;
                x8 = b+l^2+l+1;

                tri(a, :) = [x1, x2, x4, x6];
                a=a+1;
                tri(a, :) = [x4, x3, x1, x6];
                a=a+1;
                tri(a, :) = [x1, x3, x5, x6];
                a=a+1;
                tri(a, :) = [x4, x3, x8, x6];
                a=a+1;
                tri(a, :) = [x3, x7, x5, x6];
                a=a+1;
                tri(a, :) = [x3, x7, x8, x6];
                a=a+1;
            end
        end
    end
end


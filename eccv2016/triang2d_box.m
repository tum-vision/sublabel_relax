function [vert, tri] = triang2d_box(t)
%GET_TRIANGULATION_BOX Summary of this function goes here
%   Detailed explanation goes here
    
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


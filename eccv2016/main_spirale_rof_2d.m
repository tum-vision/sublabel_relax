N = 50; % number of points on the spirale
l = 3; % number of labels in each direction
t = linspace(-1, 1, l); % 1d label space
lambda = 0.5; % parameter in front of TV

%% 
% construct label-sapce
[vert, tri] = triang2d_box(t);
L = size(vert, 1);
T = size(tri, 1);

%% 
% compute spiral
spiral = zeros(N, 2);
for i=1:N
    angle = double(i) / N;
    radius = 0.9 - 0.5 * (i / N);
    px = radius * cos(angle * 4 * pi);% + 0.05 * randn;
    py = radius * sin(angle * 4 * pi);% + 0.05 * randn;
    spiral(i, :) = [ px py ];
end

%% 
% compute baseline dataterm
dataterm = zeros(N, 1, l, l);
for i=1:N
    for j=1:l
        for k=1:l
            dy = t(j) - spiral(i,2);
            dx = t(k) - spiral(i,1);
            dataterm(i, 1, j, k) = 0.5 * (dx.^2 + dy.^2);
        end
    end
end

%%
%compute precise solution

data.f = reshape(spiral, N, 1, 2);

[u_proj, u_lifted] = solve_sublabel_nd(vert, tri, N, 1, 'quad', lambda, data);

result = u_proj;

%%
% compute unlifted spiral
u = unlifted_rof(reshape(spiral, [N, 1, 2]), lambda);
result_unlifted = u; %reshape(u, [N, 1, 2]);

%% 
% plot result
figure;
hold on;

%labels
scatter(vert(:,1), vert(:,2), 1, [0.5,0.5,0.5], 'filled');

%plot triangles
for i=1:size(tri,1)
    for j=1:3
        j1 = j;
        j2 = mod(j,3) + 1;
        h = plot([vert(tri(i,j1), 1) vert(tri(i,j2), 1)], ...
                 [vert(tri(i,j1), 2) vert(tri(i,j2), 2)], ...
                 'k', 'LineWidth', 0.025);
        set(h, 'color', [0.75, 0.75, 0.75]);
    end
end
%%
% plot dataterm
scatter(spiral(:,1), spiral(:,2), 2, 'red', 'filled');
for i=1:(N-1)
    plot([spiral(i,1), spiral(i+1,1)], [spiral(i,2), spiral(i+1,2)], ...
         'r', 'LineWidth', 0.25);
end

% plot unlifted result
scatter(result_unlifted(:, 1, 1), result_unlifted(:, 1, 2), 2, ...
        'markerfacecolor', [0.0, 1.0, 0.0], 'markeredgecolor', [0.0, 1.0, 0.0]);

for i=1:(N-1)
    h = plot([result_unlifted(i,1,1), result_unlifted(i+1,1,1)], ...
         [result_unlifted(i,1,2), result_unlifted(i+1,1,2)], 'k', ...
             'LineWidth', 0.25);
    set(h, 'color', [0.0, 1.0, 0.0]);
end

%plot result
scatter(result(:, 1, 1), result(:, 1, 2), 2, 'blue', 'filled');
for i=1:(N-1)
    plot([result(i,1,1), result(i+1,1,1)], ...
         [result(i,1,2), result(i+1,1,2)], ...
         'b', 'LineWidth', 0.25);
end
%%
set(gca, 'ytick', t);
set(gca, 'xtick', t);

xlabel = cell(L, 1);
ylabel = cell(L, 1);
for i=1:N
    xlabel{i} = '';
    ylabel{i} = '';
end
xlabel{1} = int2str(t(1));
xlabel{floor(l/2)+1} = int2str((t(1)+t(end))/2);
xlabel{l} = int2str(t(end));
ylabel{1} = int2str(t(1));
ylabel{floor(l/2)+1} = int2str((t(1)+t(end))/2);
ylabel{l} = int2str(t(end));

set(gca, 'xticklabel', xlabel);
set(gca, 'yticklabel', ylabel);

grid on;
axis equal;
axis tight;
axis([t(1), t(end), t(1), t(end)]); 





function [fig_handle] = plot_net_pores2(nodes, fibers, pores, xdim, fiber_rad)


% plot_net(nodes, fibers)
%
% in:
%
% nodes            N x 3 coordinates for N nodes
% fibers           N x 2 start-end nodes for N fibers
%
% last rev:
%
% tue nov 6 2012 mfh


figure;

plot3(nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',0.2, ...
    'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.8],'MarkerSize',3);

%axis_val = max( nodes(1:end) ) + 0.1;     % find a max/min for plot
%axis([-axis_val axis_val -axis_val axis_val -axis_val axis_val]);
%axis([-1 +1 -1 +1 -1 +1]);

hold on;

for n = 1 : size(fibers, 1) % count rows
    
    x(1) = nodes(fibers(n,1), 1); % node 1 x coord
    y(1) = nodes(fibers(n,1), 2); % node 1 y coord
    z(1) = nodes(fibers(n,1), 3); % node 1 z coord
    
    x(2) = nodes(fibers(n,2), 1); % node 2 x coord
    y(2) = nodes(fibers(n,2), 2); % node 2 y coord
    z(2) = nodes(fibers(n,2), 3); % node 2 z coord
    
    plot3(x, y, z, 'b');
    
    hold on;
    
end

% Plot pores
theta = 0:pi/100:2*pi;

% Set up colorbar variables
pore_color_index = [0: 1e-9: 10e-8];
%pore_color_index = sort(pore_color_index, 'descend');  %this %inverts the scale
pore_color_index = round(pore_color_index, 9);  %this helps matlab find a match in "intersect()" below.

CMap = hot(length(pore_color_index));
CMap_fib=colormap(CMap);  %This sets the number of colors
colorbar();

ax = gca;
ax.CLim = [0 10e-8];

% Convert pore radius to real units from computational
pore_radii = pores(:,4).*xdim - fiber_rad;

for i = 1:size(pores,1)
    [X,Y,Z] = sphere;
    r = pores(i,4);
    X2 = X * r;
    Y2 = Y * r;
    Z2 = Z * r;
    
    % Plot pore center
    a = pores(i,1);
    b = pores(i,2);
    c = pores(i,3);
    plot3(a, b, c, '*k'); %, 'MarkerFaceColor','k');
    hold on
    
    % Find color for circle
    % Round your pore sizes so they match colorbar indices
    [C, is, col_index] = intersect(round(pore_radii(i), 9), pore_color_index, 'stable');
    %col_index is the index of the "fib_color_index" that shares
    %the radius value in pore_color_index. 
    if isempty(col_index)
     col_index = length(pore_color_index);  %in the case where the stretch exceeds our color bounds, give it the max
    end

%     %and now it assigns colors to the stretches:
%     set(pcircle, 'Color', CMap_fib(col_index,:));
    Color = CMap_fib(col_index,:);
    %pcircle = patch(circle_pts(1,:), circle_pts(2,:), circle_pts(3,:),C);
%     for j = 1:size(X2, 1)
%         for k = 1: size(X2, 2)
%             C(j, k, 1) = Color(1);
%             C(j, k, 2) = Color(2);
%             C(j, k, 3) = Color(3);
%         end
%     end
    
    hs = surf(X2+pores(i,1), Y2+pores(i,2), Z2+pores(i,3)); %, C);
    set(hs, 'FaceColor', Color, 'FaceAlpha',0.4, 'FaceLighting','gouraud','EdgeColor','none') %,'EdgeColor','none'
    hold on

end

h = colorbar;
set(get(h,'title'),'string','Pore \newline size');

set(gcf, 'color', 'white');
axis equal;

axis_lim = max( nodes(1:end) );

axis( [-0.5 +0.5 -0.5 +0.5 -0.5 +0.5] );

xlabel('x')'; ylabel('y'); zlabel('z');

fig_handle = getframe(gcf);
end

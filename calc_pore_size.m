%%% Calc_pore_size.m
%%% Written by Lauren Bersie-Larson, 11-25-2020
%%% This function takes in an existing network, calculates all the pore
%%% sizes present within it by sampling the network randomly 50,000 times, 
%%% finding the largest sphere that can fit at said sampling point without 
%%% intersecting a fiber, based on the methods of Stylianopoulos et al.:
%%% Stylianopoulos, T., Diop-Frimpong, B., Munn, L. L., & Jain, R. K. (2010) 
%%% Diffusion Anisotropy in Collagen Gels and Tumors: The Effect of Fiber 
%%% Network Orientation. Biophysical Journal, 99(10), 3119–3128. 
%%% https://doi.org/10.1016/J.BPJ.2010.08.065
%%% Once all pore sizes are calculated, the code then reports the average, 
%%% max, and min pore size of the network. 
%%% Plotting options to visualize these pores are also included. 

%%% Part of this code is a modified version of the code for the primitive 
%%% test "Closest Point on Line Segment to Point", and "Distance of Point 
%%% to Segment", as found on pages 127-130 of   
%%% "Real-Time Collision Detection", by Christer Ericson:
%%% Ericson, C. (2005). Chapter 5 - Basic Primitive Tests. In C. Ericson (Ed.), 
%%% Real-Time Collision Detection (pp. 125–233). San Francisco: Morgan Kaufmann. 
%%% https://doi.org/https://doi.org/10.1016/B978-1-55860-732-3.50010-3

function [avg_pore_size, min_size, max_size, pores] = calc_pore_size(nodes, fibers, x, fiber_rad)
tic
% Step 1: Find network properties
boundaries = [-0.5 +0.5 -0.5 +0.5 -0.5 +0.5];
bnd_node_nums = find_boundary_nodes(nodes,boundaries);

node1 = [nodes(fibers(:, 1), 1), nodes(fibers(:, 1), 2), nodes(fibers(:, 1), 3)];
node2 = [nodes(fibers(:, 2), 1), nodes(fibers(:, 2), 2), nodes(fibers(:, 2), 3)]; 

% Step 2: Sample network by randomly generating points within
num_pts = 50000; 
pores = zeros(num_pts, 4);
for i = 1:num_pts
    point = rand(3,1)- 0.5;
    point = point.';
    dist = zeros(size(fibers, 1), 1);
    
    % Calculate distances from point to every fiber line segment. 
    % The algorithm for this was adapted from "Real-Time Collision
    % Detection" by Christer Ericson
    % Calculations are vectorized to speed up process
    ab = node2 - node1; %node1 - node2; 
    ac = repmat(point, size(fibers,1), 1) - node1;
    bc = repmat(point, size(fibers,1), 1) - node2;
    
    % Take dot products, treating rows as the vectors
    e = dot(ac, ab, 2); 
    f = dot(ac, ac, 2);
    g = dot(ab, ab, 2);
    h = dot(bc, bc, 2);
    j = dot(ac, ac, 2) - (e .* e )./ g;
    
    % Handle cases where c (point) projects outside ab
    indices1 = find(e <= 0);
    indices2 = find(e >= g);
    
    dist(indices1) = f(indices1); % Where e <=0, use f calculation
    dist(indices2) = h(indices2); % Where e >=g, use h calculation
    
    % Otherwise, c (point) projects onto ab
    indices3 = [1:size(fibers,1)];
    indices3 = setdiff(indices3, union(indices1, indices2)); % Otherwise for all other indices
    
    dist(indices3) = j(indices3);
    
    dist = sqrt(dist); 
    
    % Make sure sphere isn't contacting fiber
    inradius = min(dist((dist.*x)>= fiber_rad));

    % Record pore size in matrix
    pores(i, :) = [point, inradius];
    
%     % Plot individual pore in network- uncomment if you want this
%     fig_handle = plot_net_pores2(nodes, fibers, pores, x, fiber_rad)
%     figure_file = ['Network_pore_',num2str(i),'.fig'];
%     im{i} = frame2im(fig_handle);
%     close
    clear point dist ab ac bc e f g h j indices1 indices2 indices3
end

% % Write all individual pore figures to a gif- uncomment if you want this
% filename = 'network_pores.gif'; % Specify the output file name
% for idx = 1:num_pts
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
%     end
% end

% Step 6: Average all radii to find average pore size
% Convert all pore sizes to real units through x
avg_pore_size = mean(pores(:,4).*x - fiber_rad)
min_size = min(pores(:,4).*x - fiber_rad)
max_size = max(pores(:,4).*x - fiber_rad)

% Step 7: Optional visualization
% plot_net_pores2(nodes, fibers, pores, x, fiber_rad)

toc
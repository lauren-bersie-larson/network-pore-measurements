%%% Calc_pore_dist.m
%%% Written by Lauren Bersie-Larson, 12-7-2020
%%% This function takes in an existing network and pore size threshold, and
%%% determines the distance a sphere the size of the pore size threshold
%%% can travel before intersecting a fiber or the end of the network. This
%%% process is done by randomly sampling the network 20,000 times. 

%%% Part of this code is a modified version of the code for the primitive 
%%% test "Closest Point on Line Segment to Point", and "Distance of Point 
%%% to Segment", as found on pages 127-130 of   
%%% "Real-Time Collision Detection", by Christer Ericson:
%%% Ericson, C. (2005). Chapter 5 - Basic Primitive Tests. In C. Ericson (Ed.), 
%%% Real-Time Collision Detection (pp. 125–233). San Francisco: Morgan Kaufmann. 
%%% https://doi.org/https://doi.org/10.1016/B978-1-55860-732-3.50010-3

%%% INPUTS: 
%%% Nodes= Nx3 matrix of the X,Y,Z coordinates of the nodes(pts) in network
%%% Fibers = MX2 matrix of the fibers in the network. The two entries
%%% correspond to the two nodes that comprise the fiber
%%% x = the dimension of the network, in meters. To input 20 micron, use 20e-6, etc. 
%%% fiber_rad = radius of the fibers, in meters.
%%% Threshold = the pore size you are using, in meters. 
%%% Direction = Either 1, 2, or 3, corresponding to X, Y, or Z,
%%% respectively. 

%%% OUTPUTS:
%%% Avg_dist = the averaged distance traveled by a sphere of size threshold, 
%%% before hitting a fiber or end of the network, 
%%% in the direction specified. This is a scalar value.
%%% dist_traveled = A vector (sized 20,000) of all the distances traveled 
%%% by spheres of size threshold, before hitting a fiber or end of the network, 
%%% in the direction specified.


function [avg_dist, dist_traveled] = calc_pore_dist(nodes, fibers, x, fiber_rad, threshold, direction)
tic
% Step 1: Find network properties
boundaries = [-0.5 +0.5 -0.5 +0.5 -0.5 +0.5];
bnd_node_nums = find_boundary_nodes(nodes,boundaries);

node1 = [nodes(fibers(:, 1), 1), nodes(fibers(:, 1), 2), nodes(fibers(:, 1), 3)];
node2 = [nodes(fibers(:, 2), 1), nodes(fibers(:, 2), 2), nodes(fibers(:, 2), 3)]; 

% Step 2: Sample network by randomly generating points within
num_pts = 20000; 
pores = zeros(num_pts, 4);
dist_traveled = zeros(1, num_pts);

for i = 1:num_pts
    % Generate random point within network
    point = rand(3,1)- 0.5;
    point = point.';
    dist = zeros(size(fibers, 1), 1); 
    
    inradius = 10; % High value to enter while loop
    loop_count = 0; 
    
    while (inradius*x - fiber_rad) > threshold
        % Update sphere position 
        point(direction) = point(direction) + loop_count*0.1*(threshold/x); 
        
        % Make sure this new position isn't outside the network
        if point(direction) < boundaries(direction*2-1) || point(direction) > boundaries(direction*2)
            break
        end
        
        % Calculate distances from point to every fiber line segment. 
        % The algorithm for this was adapted from "Real-Time Collision Detection"
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

        inradius = min(dist);

        % See if dimensionalized distance away from nearest fiber is within
        % the radius threshold to continue on path or not
        if (inradius*x - fiber_rad) >= threshold 
            % Record pore size in matrix
            pores(i, :) = [point, (threshold/x)];
            
                        
            loop_count = loop_count + 1;
            
            % Plot individual pore in network- uncomment if you want this
%             fig_handle = plot_net_pore_dist(nodes, fibers, pores, x, fiber_rad)
%             figure_file = ['Network_pore_',num2str(loop_count),'.fig'];
%             im{i, loop_count} = frame2im(fig_handle);

        else
            break % Done
        end
        
        close
          
        clear dist ab ac bc e f g h j indices1 indices2 indices3
    end
    
    if loop_count > 1
        dist_traveled(i) = loop_count-1;
    else 
        dist_traveled(i) = loop_count; 
    end
end

% % Write all individual pore figures to a gif- uncomment if you want this
% filename = 'network_pore_dist_10.gif'; % Specify the output file name
% for idx1 = 1:size(im, 1)
%     for idx2 = 1:size(im, 2)
%         if ~isempty(im{idx1, idx2})
%             [A,map] = rgb2ind(im{idx1, idx2},256);
%             if idx1 == 1 && idx2 == 1
%                 imwrite(A,map,filename,'gif','WriteMode','overwrite','LoopCount',0,'DelayTime',0.25);
%             else
%                 imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.25);
%             end
%         end
%     end
% end

% Dimensionalize resulting distance traveled
dist_traveled = dist_traveled.*(0.1*(threshold/x)); 
avg_dist = mean(dist_traveled); 

end
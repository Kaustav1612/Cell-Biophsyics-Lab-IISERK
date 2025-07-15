function [cluster_centers,cluster_assignments,features,final_object_positions] = cluster_object_space(roi,bandwidth,type,filtered_img_norm)
        [features,Intensity_map]=feature_extraction_space(roi,filtered_img_norm);
       
        % all_feature_vector is the entire dataset of the objects with rows
        % having the data for each individual objects and columns having the
        % features
        

        % Define Mean Shift Parameters
        % Weight Assign using Gaussian Kernel (Size of Kernel) to create a density map 
        epsilon = bandwidth*0.01; % Convergence Checker Parameter for shifting of the cluster center until it reaches the attraction basin
        learning_parameter = 0.1; % Rate of Learning
        merge_distance = 0.05*bandwidth; %Thresholding using the density map
        
        [final_object_positions] = shift(features,epsilon,bandwidth,learning_parameter,type,Intensity_map,filtered_img_norm);
        
        [cluster_centers,cluster_assignments] = find_clusters(final_object_positions,merge_distance);
%         figure();
%         gscatter(features(:,1),features(:,2),cluster_assignments);
%         hold on;
%         plot(cluster_centers(:,1),cluster_centers(:,2),'kx','MarkerSize', 15, 'LineWidth', 3);
%         directory='F:\Mean_Shift_Cluster\Cluster_Images\Plots\Cell_01_T1';
%         filename=sprintf('ROI_%d.fig',k);
%         savefig(fullfile(directory,filename));
%         close;

end
    
    % Shifting the position of an object by assigning weights to all the other objects with respect to a single object
    % Then getting the final shifted object positions
    
function [final_object_positions] = shift(features, epsilon, bandwidth, learning_parameter, type,Intensity_map,filtered_img_norm)
    final_object_positions = features;  % Initialize final positions as initial positions
    num_points = size(features, 1);
    max_iterations = 2000;  % Maximum number of iterations to prevent infinite loops

    for i = 1:num_points
        current_point = final_object_positions(i, :);
        current_init_pos = features(i,:);
        other_points_pos = features(~features(i,:));

        converged = false;
        iter = 0;
        
        while ~converged && iter < max_iterations
            mean_shift_vector = zeros(1, size(features, 2));
            
            % Calculate mean shift vector
            distances = pdist2(features, current_point);
%             distance_map = zeros(size(filtered_img_norm));
%             distance_matrix = pdist2(features, current_point);
%             distance_map(sub2ind(size(distance_map), features(:,1), features(:,2))) = distance_matrix;

            % weights = isotropic_gaussian_kernel(distances, bandwidth)./(Intensity_map);
            weights = sigmo_gaussian(bandwidth,distances,Intensity_map,current_init_pos,features);
            mean_shift_vector = sum(weights.* (features-current_point), 1);
            
            % Update the current point and check for convergence
            if norm(mean_shift_vector) > epsilon
                current_point = current_point + learning_parameter * mean_shift_vector;
                converged = false;
            else
                converged = true;
            end
            
            iter = iter + 1;
            if iter == max_iterations-1
                fprintf("Not converged \n");
            end
        end

        final_object_positions(i, :) = current_point;
        
        if type == 2
            if (i / num_points) > 0.25 && (i / num_points) < 0.2501
                fprintf("25 percent of Cluster Detection reached \n ");
            elseif (i / num_points) > 0.50 && (i / num_points) < 0.5001
                fprintf("50 percent of Cluster Detection reached \n");
            elseif (i / num_points) > 0.75 && (i / num_points) < 0.7501
                fprintf("75 percent of Cluster Detection reached \n");
            end          
        end
    end
    end
    
function [cluster_centers, cluster_assignments] = find_clusters(final_object_positions, merge_distance)
    [cluster_centers, flag] = merge_close_centers(final_object_positions, merge_distance);

    % Assign each point to the nearest cluster center
    cluster_assignments = zeros(size(final_object_positions, 1), 1);
    distances = pdist2(cluster_centers, final_object_positions);  % Compute all distances once

    for k = 1:size(final_object_positions, 1)
        [~, nearest_cluster_center_index] = min(distances(:, k));
        
        if flag(nearest_cluster_center_index) == 1
            cluster_assignments(k) = nearest_cluster_center_index;
        else
            cluster_assignments(k) = 0;
        end
    end
end

function [merged_centers, flag] = merge_close_centers(final_object_positions, merge_distance)
    num_centers = size(final_object_positions, 1);
    merged = false(num_centers, 1);
    flag = [];
    
    for i = 1:num_centers
        if merged(i)
            continue;
        end
        
        object_size = 0;
        
        for j = 1:num_centers  
            if norm(final_object_positions(i, :) - final_object_positions(j, :)) < merge_distance && i~=j
                object_size = object_size + 1;
                final_object_positions(j, :) = (final_object_positions(i, :) + final_object_positions(j, :)) / 2;
                merged(j) = true;
            end
        end
        
        if object_size > 5
            flag = [flag, 1];
        else 
            flag = [flag, 0];
        end
    end
    
    merged_centers = final_object_positions(~merged, :);
end


% Calculation of weight of objects using the isotropic gaussian kernel
    
    function weight = isotropic_gaussian_kernel(euclidean_distance,bandwidth)
        weight = (1/(bandwidth*sqrt(2*pi)))*exp(-((euclidean_distance).^2)/(2*bandwidth^2));
    end

    
% Calculation of weight of objects using the uniform kernel

    function weight = uniform_kernel(euclidean_distance,bandwidth)
        if euclidean_distance<=bandwidth
            weight = 1;
        else 
            weight = 0;
        end    
    end
  
 % Calculation of weight of objects using a combination of sigmoidal and gaussian function    
   function weight = sigmo_gaussian(bandwidth,distances,intensity_map,current_init_pos,features)
       intensity_weight = (exp(-intensity_map(features))/(1+exp(-intensity_map(current_init_pos(1),current_init_pos(2)))));
       weight = intensity_weight .*((1/(bandwidth*sqrt(2*pi)))*exp(-((distances).^2)/(2*bandwidth^2)));
    end

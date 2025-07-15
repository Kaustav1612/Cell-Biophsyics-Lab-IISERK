function [cluster_centers,class_assignment,classify_features,all_haralick_features,all_haralick_features_scaled] = cluster_object_shape(roi,bandwidth,type,filtered_img_norm,roi_clusteredLabeledImage,cluster_assignments,n_objects)
        [classify_features,all_haralick_features_scaled,all_haralick_features]=feature_extraction_shape(filtered_img_norm,cluster_assignments,n_objects,roi_clusteredLabeledImage);
        % all_feature_vector is the entire dataset of the objects with rows
        % having the data for each individual objects and columns having the
        % features
        

        % Define Mean Shift Parameters
        % Weight Assign using Gaussian Kernel (Size of Kernel) to create a density map 
        epsilon = bandwidth*0.1; % Convergence Checker Parameter for shifting of the cluster center until it reaches the attraction basin
        learning_parameter = 0.1; % Rate of Learning
        merge_distance = 0.01*bandwidth; %Thresholding using the density map
        
        [final_object_positions] = shift(classify_features,epsilon,bandwidth,learning_parameter,type);
        
        [cluster_centers,class_assignment] = find_clusters(final_object_positions,merge_distance);

   end
    
    % Shifting the position of an object by assigning weights to all the other objects with respect to a single object
    % Then getting the final shifted object positions
    
    function [final_object_positions] = shift(features, epsilon, bandwidth, learning_parameter, type)
    final_object_positions = features;  % Initialize final positions as initial positions
    num_points = size(features, 1);
    max_iterations = 1000;  % Maximum number of iterations to prevent infinite loops

    for i = 1:num_points
        current_point = final_object_positions(i, :);
        converged = false;
        iter = 0;
        
        while ~converged && iter < max_iterations
            mean_shift_vector = zeros(1, size(features, 2));
            
            % Calculate mean shift vector
            distances = pdist2(features, current_point);  % Calculate all distances once
            weights = isotropic_gaussian_kernel(distances, bandwidth);
            mean_shift_vector = sum((weights .* (features - current_point)), 1);
            
            % Update the current point and check for convergence
            if norm(mean_shift_vector) > epsilon
                current_point = current_point + learning_parameter * mean_shift_vector;
                converged = false;
            else
                converged = true;
            end
            
            iter = iter + 1;
        end
        
        final_object_positions(i, :) = current_point;
        
    end
    end

function [cluster_centers, class_assignments] = find_clusters(final_object_positions, merge_distance)
    [cluster_centers, flag] = merge_close_centers(final_object_positions, merge_distance);

    % Assign each point to the nearest cluster center
    class_assignments = zeros(size(final_object_positions, 1), 1);
    distances = pdist2(cluster_centers, final_object_positions);  % Compute all distances once

    for k = 1:size(final_object_positions, 1)
        [~, nearest_cluster_center_index] = min(distances(:, k));
        
        if flag(nearest_cluster_center_index) == 1
            class_assignments(k) = nearest_cluster_center_index;
        else
            class_assignments(k) = 0;
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
        
        for j = i+1:num_centers  
            if norm(final_object_positions(i, :) - final_object_positions(j, :)) < merge_distance
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
    


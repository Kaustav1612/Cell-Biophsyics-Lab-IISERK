function [object_properties, cluster_properties,roi_cluster_properties,distances] = enlist_properties_updated(cluster_centers, features, roi_cluster_label_img, roi_LabeledImage, roi, stats, componentIndices)
    num_clusters = length(componentIndices);
    num_features = size(features, 1);
    
    % Preallocate struct arrays for object properties
    object_properties = struct('Area', cell(1, num_clusters),'Radius', cell(1, num_clusters), 'Centroid', cell(1, num_clusters), ...
        'BoundingBox', cell(1, num_clusters), 'ConvexArea', cell(1, num_clusters), ...
        'Eccentricity', cell(1, num_clusters), 'Circularity', cell(1, num_clusters), ...
        'EquivDiameter', cell(1, num_clusters), 'EulerNumber', cell(1, num_clusters), ...
        'Extent', cell(1, num_clusters), 'Extrema', cell(1, num_clusters), ...
        'FilledArea', cell(1, num_clusters), 'FilledImage', cell(1, num_clusters), ...
        'Image', cell(1, num_clusters), 'MaxFeretDiameter', cell(1, num_clusters), ...
        'MaxFeretAngle', cell(1, num_clusters), 'MaxFeretCoordinates', cell(1, num_clusters), ...
        'MinFeretDiameter', cell(1, num_clusters), 'MinFeretAngle', cell(1, num_clusters), ...
        'MinFeretCoordinates', cell(1, num_clusters), 'MajorAxisLength', cell(1, num_clusters), ...
        'MinorAxisLength', cell(1, num_clusters), 'Orientation', cell(1, num_clusters), ...
        'Perimeter', cell(1, num_clusters), 'PixelIdxList', cell(1, num_clusters), ...
        'PixelList', cell(1, num_clusters), 'Solidity', cell(1, num_clusters), ...
        'SubarrayIdx', cell(1, num_clusters), 'MeanIntensity', cell(1, num_clusters), ...
        'MaxIntensity', cell(1, num_clusters), 'MinIntensity', cell(1, num_clusters), ...
        'WeightedCentroid', cell(1, num_clusters));

    % Loop through each assignment to find objects belonging to the current cluster
    distances =[];
    for j = 1:num_clusters
        object_properties(j).Area = stats(j).Area;
        object_properties(j).Radius = ((stats(j).EquivDiameter).*15)./2;
        object_properties(j).Centroid = stats(j).Centroid;
        object_properties(j).BoundingBox = stats(j).BoundingBox;
        object_properties(j).ConvexArea = stats(j).ConvexArea;
        object_properties(j).Eccentricity = stats(j).Eccentricity;
        object_properties(j).Circularity = stats(j).Circularity;
        object_properties(j).EquivDiameter = stats(j).EquivDiameter;
        object_properties(j).EulerNumber = stats(j).EulerNumber;
        object_properties(j).Extent = stats(j).Extent;
        object_properties(j).Extrema = stats(j).Extrema;
        object_properties(j).FilledArea = stats(j).FilledArea;
        % object_properties(j).FilledImage = stats(j).FilledImage;
        % object_properties(j).Image = stats(j).Image;
        object_properties(j).MaxFeretDiameter = stats(j).MaxFeretDiameter;
        object_properties(j).MaxFeretAngle = stats(j).MaxFeretAngle;
        object_properties(j).MaxFeretCoordinates = stats(j).MaxFeretCoordinates;
        object_properties(j).MinFeretDiameter = stats(j).MinFeretDiameter;
        object_properties(j).MinFeretAngle = stats(j).MinFeretAngle;
        object_properties(j).MinFeretCoordinates = stats(j).MinFeretCoordinates;
        object_properties(j).MajorAxisLength = stats(j).MajorAxisLength;
        object_properties(j).MinorAxisLength = stats(j).MinorAxisLength;
        object_properties(j).Orientation = stats(j).Orientation;
        object_properties(j).Perimeter = stats(j).Perimeter;
        object_properties(j).PixelIdxList = stats(j).PixelIdxList;
        object_properties(j).PixelList = stats(j).PixelList;
        object_properties(j).Solidity = stats(j).Solidity;
        object_properties(j).SubarrayIdx = stats(j).SubarrayIdx;
        object_properties(j).MeanIntensity = stats(j).MeanIntensity;
        object_properties(j).MaxIntensity = stats(j).MaxIntensity;
        object_properties(j).MinIntensity = stats(j).MinIntensity;
        object_properties(j).WeightedCentroid = stats(j).WeightedCentroid;
    end
    
    
    number_clusters = size(cluster_centers,1);
    % Preallocate struct arrays for cluster properties
     cluster_properties = struct('AverageIntensity', cell(1, number_clusters), ...
        'MaxIntraClusterDistance', cell(1, number_clusters), 'MomentOfInertia', cell(1, number_clusters), ...
        'Area', cell(1, number_clusters), ...
        'ClusterRadius', cell(1, number_clusters),'NearestClusterDistance',cell(1,1),'NearClusterNumber',cell(1,1));

    % Get linear indices for the ROI cluster labels
    roi_indices = sub2ind(size(roi_cluster_label_img), features(:,1), features(:,2));
    feature_cluster_labels = roi_cluster_label_img(roi_indices);
    feature_label_image = roi_LabeledImage(roi_indices);

    % Initialize variables
    excluded_indices = NaN(number_clusters, 1);
    near_cluster_number = zeros(number_clusters, 1);
    

    % Loop through each cluster
    for i = 1:number_clusters
        distance = [];
        % Get the current cluster's feature indices
        cluster_indices = feature_cluster_labels == i;
        current_features = features(cluster_indices, :);

        % Get features and labels from other clusters
        other_indices = feature_cluster_labels ~= i;
        other_features = features(other_indices, :);
        other_labels = feature_cluster_labels(other_indices);  % Get labels of other features

        % Calculate distances to other clusters
        if ~isempty(current_features) && ~isempty(other_features)
            cluster_distances = pdist2(current_features, other_features).* 15;  % Scale by 15

            % Track which other clusters are within 100 distance
            nearby_clust_dist = [];
            clusters_within_threshold = [];
            min_distance_threshold = 100;
            % Find valid distances within the threshold
           for k = 1:size(current_features, 1)
                % Get valid distances for current feature
                valid_distances = cluster_distances(k, :);
                valid_indices = find(valid_distances <=  min_distance_threshold);

                % Get the labels of the features within the valid distance
                nearby_cluster_labels = other_labels(valid_indices);

             
                % Store unique nearby clusters within the threshold
                clusters_within_threshold = unique([clusters_within_threshold; nearby_cluster_labels]);
                
                
                 
           end

            % Store the number of unique clusters within 150 distance
         
            near_cluster_number(i) = numel(clusters_within_threshold);
            if ~isempty(clusters_within_threshold)
                    % Find minimum distance for each cluster label
                    for m = 1:numel(clusters_within_threshold)
                        % Get the minimum distance for this label
                        current_label = clusters_within_threshold(m);
                        label_distances = valid_distances(other_labels == current_label);
                        min_dist = min(label_distances); % Minimum distance for this label

                        % Store the minimum distance
                        distance = [distance, min_dist];
                    end
            end
          
            
        end


        % Calculate intra-cluster distances
        intra_cluster_distances = pdist2(current_features, current_features) * 15;

        % Calculate intensity, area, and moment of inertia
        intensity = roi(roi_indices(cluster_indices));
        area = sum(cluster_indices);
        moi = sum((vecnorm(current_features - cluster_centers(i,:), 2, 2) * 15).^2 .* intensity);


        % Assign properties
        if ~isnan(mean(intensity))
            cluster_properties(i).AverageIntensity = mean(intensity);
            cluster_properties(i).MaxIntraClusterDistance = max(intra_cluster_distances(:));
            cluster_properties(i).MomentOfInertia = moi;
            cluster_properties(i).Area = area;
            cluster_properties(i).ClusterRadius = sqrt((225 * area) / pi);
            cluster_properties(i).NearClusterNumber = near_cluster_number(i)';
            cluster_properties(i).NearestClusterDistance = distance';
        end
    end
    roi_cluster_properties = struct('AverageIntensity', cell(1, 1), ...
        'MaxIntraClusterDistance', cell(1, 1), 'MomentOfInertia', cell(1, 1), ...
        'Area', cell(1, 1), 'NearestClusterDistance', cell(1, 1), ...
        'ClusterRadius', cell(1, 1),'NearClusterNumber',cell(1,1));


        roi_cluster_properties.AverageIntensity = median([cluster_properties.AverageIntensity]);
        roi_cluster_properties.MaxIntraClusterDistance = median([cluster_properties.MaxIntraClusterDistance]);
        roi_cluster_properties.MomentOfInertia = median([cluster_properties.MomentOfInertia]);
        roi_cluster_properties.Area = median([cluster_properties.Area]);
        roi_cluster_properties.NearestClusterDistance = median(distance);
        roi_cluster_properties.ClusterRadius = median([cluster_properties.ClusterRadius]);
        roi_cluster_properties.NearClusterNumber = mean(near_cluster_number);
      
end

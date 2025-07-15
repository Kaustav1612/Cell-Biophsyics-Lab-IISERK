

function [shape_class_features,all_haralick_features_scaled,all_haralick_features] = feature_extraction_shape(filtered_img_norm,cluster_assignments,n_objects,roi_clusteredLabeledImage)
    
        classify_features = zeros(n_objects, 14);  % Pre-allocate for 14 Haralick features
        fixed_size = [50,50];

        for k = 1 : size(cluster_assignments)
            [px] = find(roi_clusteredLabeledImage==k);
            padded_image = zeros(size(filtered_img_norm));
            padded_image(ind2sub(size(filtered_img_norm),px))=filtered_img_norm(ind2sub(size(filtered_img_norm),px));
            stats = regionprops(padded_image>0,"BoundingBox","PixelIdxList","Area");
            componentIndices=[];
            for i = 1:length(stats)
                currentArea = stats(i).Area;
                % Check if the area is within the specified thresholds
                if currentArea >= 5  
                    componentIndices = [componentIndices, i];
                elseif currentArea >= 2000
                    continue;
                end
            end
            
            for j = 1:length(componentIndices)
                bounding_box = stats(componentIndices(j)).BoundingBox;
                min_row = floor(bounding_box(2));
                max_row = min_row + bounding_box(4);
                min_col = floor(bounding_box(1));
                max_col = min_col + bounding_box(3);


                % Extract sub-image using bounding box
                if min_row ==0
                    min_row = min_row+1;
                elseif min_col==0
                    min_col = min_col+1;
                end
                    if max_col > size(filtered_img_norm,2) || max_row > size(filtered_img_norm,1)
                        continue;
                    else
                        if max_col == 0 || max_row == 0
                            continue;
                        else
                            sub_image = filtered_img_norm(min_row:max_row, min_col:max_col);
                            [rows,cols]=size(sub_image);
                            % Check if sub_image is empty
                            if ~isempty(sub_image)
                                if rows < fixed_size(1) || cols < fixed_size(2)
                                    padded_image = zeros(fixed_size);
                                    padded_image(1:rows, 1:cols) = sub_image;
                                else
                                    padded_image = sub_image;
                                end
                            end
                            

                            x = haralickTextureFeatures(padded_image);

                                % Store the features
                                classify_features(k, :) = x;
                %                 figure();
                %                 imshow(padded_image,[]);
                %                 close;
                        end
                    end
                end
          end
        
               
         class_features=[
         classify_features(:,3),...
         classify_features(:,5),...
         classify_features(:,6),...
         classify_features(:,8),...
         classify_features(:,11),...
         classify_features(:,14)];

     scaled_class_features = zeros(size(class_features));
     
     scaled_class_features(:, 3) = robust_scaling(classify_features(:, 3));
     scaled_class_features(:, 5) = robust_scaling(classify_features(:, 5));
     scaled_class_features(:, 6) = robust_scaling(classify_features(:, 6));
     scaled_class_features(:, 8) = robust_scaling(classify_features(:, 8));
     scaled_class_features(:, 11) = robust_scaling(classify_features(:, 11));
     scaled_class_features(:, 14) = robust_scaling(classify_features(:, 14));

     shape_class_features = [scaled_class_features(:, 3).*scaled_class_features(:, 5), ...
                        scaled_class_features(:, 6).*scaled_class_features(:, 8), ...
                        scaled_class_features(:, 11).*scaled_class_features(:, 14), ...
                        ];
                    
      all_haralick_nonscaled = [
         classify_features(:,1),...
         classify_features(:,2),...
         classify_features(:,3),...
         classify_features(:,4),...
         classify_features(:,5),...
         classify_features(:,6),...
         classify_features(:,7),...
         classify_features(:,8),...
         classify_features(:,9),...
         classify_features(:,10),...
         classify_features(:,11),...
         classify_features(:,12),...
         classify_features(:,13),...
         classify_features(:,14),...
          ];
      
      all_haralick_features_scaled =  zeros(size(all_haralick_nonscaled));
      
     all_haralick_features_scaled(:, 1) = robust_scaling(classify_features(:, 1));
     all_haralick_features_scaled(:, 2) = robust_scaling(classify_features(:, 2));
     all_haralick_features_scaled(:, 3) = robust_scaling(classify_features(:, 3));
     all_haralick_features_scaled(:, 4) = robust_scaling(classify_features(:, 4));
     all_haralick_features_scaled(:, 5) = robust_scaling(classify_features(:, 5));
     all_haralick_features_scaled(:, 6) = robust_scaling(classify_features(:, 6));
     all_haralick_features_scaled(:, 7) = robust_scaling(classify_features(:, 7));
     all_haralick_features_scaled(:, 8) = robust_scaling(classify_features(:, 8));
     all_haralick_features_scaled(:, 9) = robust_scaling(classify_features(:, 9));
     all_haralick_features_scaled(:, 10) = robust_scaling(classify_features(:, 10));
     all_haralick_features_scaled(:, 11) = robust_scaling(classify_features(:, 11));
     all_haralick_features_scaled(:, 12) = robust_scaling(classify_features(:, 12));
     all_haralick_features_scaled(:, 13) = robust_scaling(classify_features(:, 13));
     all_haralick_features_scaled(:, 14) = robust_scaling(classify_features(:, 14));
     
     all_haralick_features = zeros(size(all_haralick_nonscaled));
     
     all_haralick_features(:, 1) = (classify_features(:, 1));
     all_haralick_features(:, 2) = (classify_features(:, 2));
     all_haralick_features(:, 3) = (classify_features(:, 3));
     all_haralick_features(:, 4) = (classify_features(:, 4));
     all_haralick_features(:, 5) = (classify_features(:, 5));
     all_haralick_features(:, 6) = (classify_features(:, 6));
     all_haralick_features(:, 7) = (classify_features(:, 7));
     all_haralick_features(:, 8) = (classify_features(:, 8));
     all_haralick_features(:, 9) = (classify_features(:, 9));
     all_haralick_features(:, 10) = (classify_features(:, 10));
     all_haralick_features(:, 11) = (classify_features(:, 11));
     all_haralick_features(:, 12) = (classify_features(:, 12));
     all_haralick_features(:, 13) = (classify_features(:, 13));
     all_haralick_features(:, 14) = (classify_features(:, 14));
    
end

function scaled_data = z_score_scaling(data)
    % Standardize data to have a mean of 0 and a standard deviation of 1
    mean_data = mean(data);
    std_data = std(data);
    scaled_data = (data - mean_data) / std_data;
end

function scaled_features = robust_scaling(features)
%ROBUST_SCALING Scales features using robust scaling.
%
%   scaled_features = ROBUST_SCALING(features) scales the features using robust
%   scaling, which is less sensitive to outliers.
%
%   Inputs:
%   - features: A matrix of features to be scaled.
%
%   Outputs:
%   - scaled_features: A matrix of scaled features.

% Calculate the median and interquartile range (IQR)
median_values = median(features, 1);
iqr = quantile(features, 0.75, 1) - quantile(features, 0.25, 1);

% Scale the features
scaled_features = (features - median_values) ./ iqr;
end
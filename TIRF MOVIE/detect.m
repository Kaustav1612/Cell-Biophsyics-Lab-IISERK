%% Detect Function

function [componentIndices,stats_label,labeled_image] = detect(normalized_roi,lowerThreshold,upperThreshold)
    

    L = bwlabel(normalized_roi);
   
    % Create a labeled image initialized with the original image
    labeled_image = zeros(size(normalized_roi));


    % Get properties using a single regionprops call
   
    stats_label = regionprops(L,normalized_roi, 'all');
    
    
    % Initialize component indices and areas arrays
    componentIndices = [];
    componentAreas = [];
    
    % Loop through all connected components
    for i = 1:length(stats_label)
        currentArea = stats_label(i).Area;
        % Check if the area is within the specified thresholds
        if currentArea >= lowerThreshold && currentArea<upperThreshold
            % Access object pixels using PixelIdxList
            object_pixels = stats_label(i).PixelIdxList;
            labeled_image(object_pixels) = i;  % Label with object index
            % Add the index of the component to the componentIndices array
            componentIndices = [componentIndices, i];
            componentAreas = [componentAreas, currentArea];
        end
    end
    
    % Number of objects detected
    n_objects = length(componentIndices);
end
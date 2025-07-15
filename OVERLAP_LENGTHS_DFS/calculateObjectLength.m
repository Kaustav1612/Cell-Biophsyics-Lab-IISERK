function all_lengths = calculateObjectLength(image, componentIndices, stats_label)
    all_lengths = []; % Initialize array to store lengths of all objects
    fixed_size = [30, 30]; % Fixed size for padded sub-image

    % Loop through all detected objects
    num_objects = numel(componentIndices);
    for i = 1:num_objects
        % Extract pixel indices for the current object
         % Extract pixel indices for the current object
    pixelList = stats_label(componentIndices(i)).PixelIdxList;
    
    % Convert the pixel indices to row and column subscripts
    [rows, cols] = ind2sub(size(image), pixelList);
    
    % Determine the bounding box coordinates using the minimum and maximum row and column values
    min_row = min(rows);
    max_row = max(rows);
    min_col = min(cols);
    max_col = max(cols);

    % Calculate the size of the bounding box
    bbox_rows = max_row - min_row + 1;
    bbox_cols = max_col - min_col + 1;

    % Create a zero-padded image of fixed size
    padded_image = zeros(fixed_size);

    % Initialize an empty sub-image within the bounding box
    sub_image = zeros(bbox_rows, bbox_cols);
    
    % Adjust the row and column indices relative to the bounding box
    adjusted_rows = rows - min_row + 1;
    adjusted_cols = cols - min_col + 1;

    % Populate the sub-image with the object pixels
    linear_indices_within_bbox = sub2ind([bbox_rows, bbox_cols], adjusted_rows, adjusted_cols);
    sub_image(linear_indices_within_bbox) = image(pixelList);
    

    % Check if sub_image needs to be padded to match the fixed size
    [sub_image_rows, sub_image_cols] = size(sub_image);
    if sub_image_rows < fixed_size(1) || sub_image_cols < fixed_size(2)
        % Zero-pad the sub-image to the fixed size
        padded_image(1:sub_image_rows, 1:sub_image_cols) = sub_image;
    else
        % If sub_image is already the fixed size or larger, use it directly
        padded_image = sub_image;
    end

        % Skeletonize the image to reduce thickness
        skel_image = bwmorph(padded_image, 'skel', Inf);
        padded_skel_image = padarray(skel_image, [0, 0], 0, 'both');
        branchpoints = bwmorph(padded_skel_image, 'branchpoints');
        endpoints=bwmorph(padded_skel_image, 'endpoints');
%         figure(4);
%         imshow(branchpoints,[]);
%         title('Branch Points');
%         figure(5);
%         imshow(endpoints,[]);
%         title('Endpoints');
%         figure(6);
%         imshow(skel_image,[]);


        % Find a starting pixel within the skeletonized object
        [row_endpoints, col_endpoints] = find(endpoints==1);
        [row_branchpoints, col_branchpoints] = find(endpoints==1);

        % Initialize visited matrix and find starting pixel
        visited = false(size(skel_image));
        start_points = [row_endpoints, col_endpoints];
        branch_points = [row_branchpoints, col_branchpoints];
        object_lengths = [];
        for i = 1: size(start_points)
            object_lengths = [];
            startPixel=[row_endpoints(i), col_endpoints(i)];
            visited = NaN(size(skel_image));
            indices = find(skel_image);
            visited(indices) = 0;
            
            [object_lengths, ~, visited] = dfs_traversal_singleBranch(skel_image, startPixel, visited, object_lengths,branch_points,start_points);
        end
              
        all_lengths = [all_lengths, object_lengths];
    end
end

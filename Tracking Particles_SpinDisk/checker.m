    % Clear workspace and close figures
    clc; clear; close all;
    
    % Initialize parameters
    params = struct();
    params.MinArea = 12;          % Minimum object area (pixels)
    params.MaxArea = 1000;        % Maximum object area (pixels)
    params.FrameStart = 1;        % First frame to process
    params.FrameEnd = 56;         % Last frame to process
    params.OverlapThreshold = 0.5;% Threshold for object matching
    nframes = params.FrameEnd-params.FrameStart+1;
    % Main processing pipeline
   
    [image_stack, mask,marked_frame,data_dir] = load_and_mask_images(params);
    [all_stats, detected_movie,intensity_thresh] = detect_objects(image_stack,marked_frame, params);
    new_stats = update_all_stats(all_stats,params);   
    all_objects = track_object_lifetimes(new_stats,nframes,params);

    %% --- Generate Animated GIF from Detected Movie ---
    gifFile = fullfile(data_dir, 'binary_detection.gif');
    for i = 1:nframes
        [A, map] = gray2ind(mat2gray(detected_movie(:, :, i)), 256);
        if i == 1
            imwrite(A, map, gifFile, 'LoopCount', Inf, 'DelayTime', 1);
        else
            imwrite(A, map, gifFile, 'WriteMode', 'append', 'DelayTime', 1);
        end
    end


gifFile_annotated = fullfile(data_dir,'annotated_detection.gif');
% Precompute frame-to-object map
frame_object_map = cell(nframes, 1);

for t = 1:numel(all_objects)
    obj = all_objects(t);
    for frame_idx = obj.start_frame : obj.end_frame
        idx = frame_idx - params.FrameStart + 1;
        time_relative = frame_idx - obj.start_frame + 1;

        if idx > 0 && idx <= nframes && time_relative <= length(obj.XC)
            x = obj.XC(time_relative);
            y = obj.YC(time_relative);

            % Store the object ID at appropriate index
            if isempty(frame_object_map{idx})
                frame_object_map{idx} = [abs(obj.id-size(all_objects,2))];
            else
                frame_object_map{idx}(end+1) = abs(obj.id-size(all_objects,2));
            end

        end
    end
end


for i = 1:nframes
    current_frame = mat2gray(detected_movie(:, :, i));
    rgb_frame = repmat(current_frame, [1, 1, 3]);

    % Overlay object IDs based on precomputed mapping
    stats = new_stats{i};
    id_map = frame_object_map{i};
    for k = 1:numel(stats)
        if ~isempty(id_map) && length(id_map) >= k && id_map(k) ~= 0
            id = id_map(k);
            pos = round(stats(k).Centroid);
            rgb_frame = insertText(rgb_frame, pos, num2str(id), 'TextColor', 'white', 'BoxOpacity', 0.1, 'FontSize', 10);
        end
    end

    % Convert to indexed image for GIF
    if i == 1
        global_map = gray(256); 
        
        [A, ~] = rgb2ind(rgb_frame, global_map, 'nodither'); % Use 'nodither' for sharper text
        imwrite(A, global_map, gifFile_annotated, 'LoopCount', Inf, 'DelayTime', 0.5);
    else
        [A, ~] = rgb2ind(rgb_frame, global_map, 'nodither'); % Use the same global_map
        imwrite(A, global_map, gifFile_annotated, 'WriteMode', 'append', 'DelayTime', 0.5);
    end
end

%% Image Loading and Masking
function [image_stack, cell_mask,marked_frame,data_dir] = load_and_mask_images(params)
    % Get directory containing TIFF files
    data_dir = uigetdir('', 'Select directory containing TIFF files');
    if data_dir == 0
        error('No directory selected.');
    end
    
    % Load image stack
    tiff_files = dir(fullfile(data_dir, '*.tif'));
    if isempty(tiff_files)
        error('No TIFF files found in the directory.');
    end
    
    % Read first frame to get dimensions
    first_frame = imread(fullfile(data_dir, tiff_files(1).name));
    [height, width] = size(first_frame);
    
    % Initialize image stack
    nframes = min(params.FrameEnd, length(tiff_files)) - params.FrameStart + 1;
    image_stack = zeros(height, width, nframes, 'like', first_frame);
    
    % Load all frames
    for i = params.FrameStart:nframes
        image_stack(:,:,i) = imread(fullfile(data_dir, tiff_files(params.FrameStart + i - 1).name));
    end

    image_stack = freehand_rect_crop(image_stack,params);
    % Create cell mask
    cell_mask = create_cell_mask(image_stack(:,:,params.FrameStart));

    % Only need if mask is complicated to draw
    for k = 1 : (params.FrameEnd-params.FrameStart+1)
        frame = image_stack(:,:,k);
        frame(~cell_mask) = 0;
        freehand_masked_stack(:,:,k) = frame;
    end
    image_stack = freehand_masked_stack;
    
    % Mark events if needed
    if questdlg('Mark events in first frame?', 'Event Marking', 'Yes', 'No', 'No') == "Yes"
        marked_frame = mark_events(image_stack(:,:,params.FrameStart), cell_mask);
    end
end

function mask = create_cell_mask(frame)
    fig = figure();
    imshow(frame, []);
    title('Draw the cell boundary');
    
    redrawing = true;
    while redrawing
        h = drawfreehand('Color', 'r', 'LineWidth', 2);
        choice = questdlg('Keep this mask?', 'Mask Confirmation', 'Yes', 'No', 'Yes');
        
        if strcmp(choice, 'Yes')
            mask = createMask(h);
            redrawing = false;
        else
            delete(h);
        end
    end
    close(fig);
end

function marked_frame = mark_events(frame, cell_mask)
    fig = figure();
    frame(~cell_mask)=0;
    imshow(frame, []);
    title('Mark events (e.g., vesicles)');
    marked_frame = frame;
    
    while true
        choice = questdlg('Mark an event?', 'Event Marking', 'Yes', 'No', 'No');
        if strcmp(choice, 'No')
            break;
        end
        
        h = drawfreehand('Color', 'g', 'LineWidth', 1);
        event_mask = createMask(h);
        marked_frame(event_mask) = 0;  % Set events to background
    end
    close(fig);
end

%% Object Detection
function [all_stats, detected_movie,intensity_thresh] = detect_objects(image_stack, marked_frame, params)
    nframes = size(image_stack, 3);
    detected_movie = zeros(size(image_stack), 'like', image_stack);
    all_stats = cell(nframes, 1);
    
    % Calculate intensity threshold from cell region
    intensity_thresh = max(marked_frame(:));
    
    for i = 1:nframes
        % Threshold and clean up binary image
        binary_frame = image_stack(:,:,i) > intensity_thresh;
        binary_frame = bwareaopen(binary_frame, params.MinArea);
        binary_frame = imfill(binary_frame, 'holes');
        
        % Connected components analysis
        cc = bwconncomp(binary_frame);
        stats = regionprops(cc, 'Area', 'Centroid', 'PixelIdxList', ...
            'Eccentricity', 'EquivDiameter', 'MajorAxisLength', 'MinorAxisLength');
        
        % Filter objects by size and calculate additional properties
        valid_objects = [];
        for j = 1:numel(stats)
            if stats(j).Area >= params.MinArea && stats(j).Area <= params.MaxArea
                % Calculate intensity features
                pixels = stats(j).PixelIdxList;
                intensity_values = image_stack(pixels + (i-1)*numel(image_stack(:,:,1)));
                
                stats(j).MeanIntensity = mean(intensity_values);
                stats(j).MaxIntensity = max(intensity_values);
                stats(j).MinIntensity = min(intensity_values);
                stats(j).StdIntensity = std(double(intensity_values));
                
                % Store spatial features
                stats(j).XC = stats(j).Centroid(1);
                stats(j).YC = stats(j).Centroid(2);
                
                % Add to valid objects
                valid_objects = [valid_objects; stats(j)];
                
                % Create labeled image
                detected_movie(pixels + (i-1)*numel(image_stack(:,:,1))) = ...
                stats(j).MeanIntensity;
            end
        end
        
        all_stats{i} = valid_objects;
        
        % Display progress
        fprintf('Processed frame %d/%d - Found %d objects\n', i, nframes, numel(valid_objects));
    end
end

function cropped_stack = freehand_rect_crop(image_stack,params)

    % Display first frame for ROI selection
    fig = figure('Name', 'Select Rectangular ROI', 'NumberTitle', 'off');
    imshow(image_stack(:,:,params.FrameStart), []);
    title('Draw rectangular crop region');
    
    % Initialize crop confirmation flag
    crop_confirmed = false;
    
    while ~crop_confirmed
        % Create interactive rectangle ROI
        h = drawrectangle('Color', 'r', 'LineWidth', 2, ...
                         'DrawingArea', [1, 1, size(image_stack,2), size(image_stack,1)]);
        
        % Ask for confirmation
        choice = questdlg('Keep this crop region?', ...
                         'Crop Confirmation', ...
                         'Yes','No','No');
        
        % Process choice
        if strcmp(choice, 'Yes')
            crop_confirmed = true;
        else
            delete(h);  % Delete the current ROI to allow redrawing
        end
    end
    
    % Get ROI position [x y width height]
    roi_position = round(h.Position);
    delete(h);
    close(fig);
    
    % Validate ROI
    if any(roi_position(3:4) < 5)
        error('Selected ROI too small - minimum 5x5 pixels required');
    end
    
    % Apply crop to all frames
    cropped_stack = image_stack(...
        roi_position(2):roi_position(2)+roi_position(4)-1, ...
        roi_position(1):roi_position(1)+roi_position(3)-1, ...
        :);
    
    % Display confirmation
    fprintf('Cropped stack from %dx%d to %dx%d pixels\n', ...
            size(image_stack,1), size(image_stack,2), ...
            size(cropped_stack,1), size(cropped_stack,2));
end

% This code helps you analyze a movie of Tf. Any one of the can be used as
% a referrence see  variable 'frametostart' on line 34
% objects are detected. centring on centroids, 13 pixel linescan in
% performed and kymograph built for each such puncta. (For Changing linescan width change variable 
% 'line_scan_length' on line 138)
% Min_Area_Thresh is the threshold of pixels in objects for considering them as
% puncta.

% This code is modularised with all the function are written as separate
% MATLAB function scripts please keep all the MATLAB Scripts in a single directory

clearvars -except filter_mobile;
close all;
clc;


% Give your correct data directory
olddir=cd('F:\Uday_data\02_03_24(CytoD_live_img2)\0.1um_cytoD\Cell 1\TNFR1_1');
directory = "F:\Uday_data\mobility_analysis\Cyto_live_img_2\0.1 microM\01\";

frametostart=1;
A=dir('*.tif');
nframes = 300; 
% nframes=30;
for i=1:nframes
    filename=A(i).name;
    image_stack (:,:,i)=imread(filename);
end
hold on;
redrawing = true;
while redrawing
  imshow(image_stack (:,:,1),[95 130]);  % Display the black and white image
  title('Draw the cell bounary');
  % Allow the user to select only the cell of importance
  mask_handle = drawfreehand('Color','r','LineWidth',2);
  
  % Option to redraw or exit
  choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
  
  % Update flag based on user choice
  redrawing = strcmp(choice, "Yes");
  
  % Clear the previous mask for redraw (if needed)
  if redrawing
    delete(mask_handle);
  end
end
% Create a binary mask from the freehand boundary
mask = mask_handle.createMask;
close;
cd(olddir);  
image_1 =  image_stack (:,:,1);
% ---- Background Selection ----
continue_drawing = true;
while continue_drawing
    imshow(image_1,[95 130]);  % Display the black and white image
    title('Draw the Background');
    h = drawrectangle('InteractionsAllowed', 'all', 'Color', 'r');
    
    % Extract top-left corner and dimensions
    top_left = h.Position(1:2);
    width = h.Position(3);
    height = h.Position(4);
    
    % Ensure the shape is a square
    side_length = min(width, height);
    square_coords = [top_left;
        top_left + [side_length, 0];
        top_left + [side_length, side_length];
        top_left + [0, side_length];
        top_left];  % Close the square
    
    
    % Confirm the ROI or redraw
    choice = questdlg('Do you want to keep this ROI or redraw?', ...
        'Confirm ROI', ...
        'Keep', 'Redraw', 'Keep');
    if strcmp(choice, 'Redraw')
        delete(h); % Remove the drawn rectangle
        delete(square_plot); % Remove the displayed red square
        continue;  % Restart the loop
    else
        continue_drawing=false;
    end
end
close;

max_row = round(square_coords(3, 2));
min_row = round(square_coords(2, 2));
max_col = round(square_coords(2, 1));
min_col = round(square_coords(1, 1));


bg_img_sted = image_1(min_row:max_row, min_col:max_col);


thresh = mean(bg_img_sted(:));
%% GET ALL CENTROIDS OR READ ALL THE CENTROIDS

for i =1 : nframes 
    current_frame = image_stack (:,:,i);
 
    
    % Set the threshold calibration which is the number of frames after
    % which sensitivity for adaptthresh is recalibrated based on the number of focused clustera
    
    threshold_calibration = 150;
    
      if rem(i,threshold_calibration)==1        
        Min_Area_Thresh = 10;
       Max_Area_Thresh = 100; 
        threshIntensity = thresh;
        correlation_threshold = 0.65;
        final_image =  zeros(size(current_frame));
        threshold = [];
        for k = 6 : 12    
            sigma = k;
            line_scan_length = k-2;
            [final_image,gauss_blur,threshold] = ObjectDetector_un(current_frame,correlation_threshold,line_scan_length,mask,sigma,final_image,threshIntensity,Min_Area_Thresh,Max_Area_Thresh);
        end
        final_image(~mask) = 0;
        eroded_image = imerode(final_image, strel('disk',1));
        output_img=imdilate(eroded_image,strel('disk',1));
        
        % Load the connected objects image (binary image)
        connected_objects = output_img; % Replace with your image file
        ground_truth_mask = imbinarize(connected_objects); % Ensure it's binary

        % Parameters for optimization
        sensitivity_range = 0.1:0.05:0.9; % Sensitivity range to test 
        best_iou = 0;

        % Iterate over sensitivity values
        for sensitivity = sensitivity_range
            % Apply adaptive thresholding
            adaptive_mask = imbinarize(current_frame, 'adaptive', 'Sensitivity', sensitivity);

            % Calculate Intersection over Union (IoU)
            intersection = sum(adaptive_mask & ground_truth_mask, 'all');
            union = sum(adaptive_mask | ground_truth_mask, 'all');
            iou = intersection / union;

            % Track the best sensitivity value
            if iou > best_iou
                best_iou = iou;
                best_sensitivity = sensitivity;
            end
        end
        
        optimal_thresholded_frame = imbinarize(current_frame, 'adaptive', 'Sensitivity', best_sensitivity);
        eroded_image = imerode(optimal_thresholded_frame, strel('disk',1));
        dilated_image=imdilate(eroded_image,strel('disk',1));
        dilated_image(~mask)=0;
        binary_stack(:,:,i) = dilated_image;   
    else
        
        % Visualize the best thresholded mask
        optimal_thresholded_frame = imbinarize(current_frame, 'adaptive', 'Sensitivity', best_sensitivity);
        eroded_image = imerode(optimal_thresholded_frame, strel('disk',1));
        dilated_image=imdilate(eroded_image,strel('disk',1));
        dilated_image(~mask)=0;
        % Store processed binary image
        binary_stack(:,:,i) = dilated_image;
      end
      
      % Connected component analysis
      connected_components = bwconncomp(dilated_image);
      stats = regionprops(connected_components, 'all');
      
      % Compute intensity and store additional properties
      for u = 1:numel(stats)
          pixel_indices = stats(u).PixelIdxList;
          intensity = mean(current_frame(pixel_indices));
          stats(u).Intensity = intensity;
          centroid = stats(u).Centroid;
          stats(u).XC = centroid(1);
          stats(u).YC = centroid(2);
          
          % Extract all other properties
          stats(u).Area = stats(u).Area;
          stats(u).Circularity = stats(u).Circularity;
          stats(u).EquivDiameter = stats(u).EquivDiameter;
          stats(u).BoundingBox = stats(u).BoundingBox;
          stats(u).EulerNumber = stats(u).EulerNumber;
          stats(u).MinorAxisLength = stats(u).MinorAxisLength;
          stats(u).Extent = stats(u).Extent;
          stats(u).Orientation = stats(u).Orientation;
          stats(u).Extrema = stats(u).Extrema;
          stats(u).Perimeter = stats(u).Perimeter;
          stats(u).ConvexArea = stats(u).ConvexArea;
          stats(u).FilledArea = stats(u).FilledArea;
          stats(u).PixelIdxList = stats(u).PixelIdxList;
          stats(u).ConvexHull = stats(u).ConvexHull;
          stats(u).FilledImage = stats(u).FilledImage;
          stats(u).PixelList = stats(u).PixelList;
          stats(u).ConvexImage = stats(u).ConvexImage;
          stats(u).Image = stats(u).Image;
          stats(u).Solidity = stats(u).Solidity;
          stats(u).Eccentricity = stats(u).Eccentricity;
          stats(u).MajorAxisLength = stats(u).MajorAxisLength;
          stats(u).SubarrayIdx = stats(u).SubarrayIdx;
      end
      
      labeled_image = zeros(size(binary_stack(:,:,i)));
      labeled_img_binary =  zeros(size(binary_stack(:,:,i)));
      labeled_image_norm = zeros(size(binary_stack(:,:,i)));
      for j = 1 : numel(stats)
          % Check if the area is within the specified thresholds
            if stats(j).Area >= Min_Area_Thresh 
                % Access object pixels using PixelIdxList
                object_pixels = stats(j).PixelIdxList;
                labeled_image_norm(object_pixels) = current_frame(object_pixels);
                labeled_image(object_pixels) = j;  % Label with object index
                labeled_img_binary(object_pixels)=1;
            end
      end
      labeled_img=bwlabeln(labeled_image);
      labeled_norm_img = bwlabeln(labeled_image_norm);
      
      
      % Store final stats
      final_stats = regionprops(labeled_norm_img, 'all');
      for u = 1:numel(final_stats)
          pixel_indices = final_stats(u).PixelIdxList;
          intensity = mean(current_frame(pixel_indices));
          final_stats(u).Intensity = intensity;
          centroid = final_stats(u).Centroid;
          final_stats(u).XC = centroid(1);
          final_stats(u).YC = centroid(2);
          final_stats(u).EquivDiameter=final_stats(u).EquivDiameter;
          
          final_stats(u).Area = final_stats(u).Area;
          final_stats(u).Circularity = final_stats(u).Circularity;
          final_stats(u).EquivDiameter = final_stats(u).EquivDiameter;
          final_stats(u).BoundingBox = final_stats(u).BoundingBox;
          final_stats(u).EulerNumber = final_stats(u).EulerNumber;
          final_stats(u).MinorAxisLength = final_stats(u).MinorAxisLength;
          final_stats(u).Extent = final_stats(u).Extent;
          final_stats(u).Orientation = final_stats(u).Orientation;
          final_stats(u).Extrema = final_stats(u).Extrema;
          final_stats(u).Perimeter = final_stats(u).Perimeter;
          final_stats(u).ConvexArea = final_stats(u).ConvexArea;
          final_stats(u).FilledArea = final_stats(u).FilledArea;
          final_stats(u).PixelIdxList = final_stats(u).PixelIdxList;
          final_stats(u).ConvexHull = final_stats(u).ConvexHull;
          final_stats(u).FilledImage = final_stats(u).FilledImage;
          final_stats(u).PixelList = final_stats(u).PixelList;
          final_stats(u).ConvexImage = final_stats(u).ConvexImage;
          final_stats(u).Image = final_stats(u).Image;
          final_stats(u).Solidity = final_stats(u).Solidity;
          final_stats(u).Eccentricity = final_stats(u).Eccentricity;
          final_stats(u).MajorAxisLength = final_stats(u).MajorAxisLength;
          final_stats(u).SubarrayIdx = final_stats(u).SubarrayIdx;
      end
      
        % Save stats for current frame in an all_stats array
        all_stats{i} = final_stats;
        labeled_movie_norm(:,:,i)=labeled_image_norm;
        detected_movie(:,:,i)= labeled_image_norm>0;
        labeled_movie(:,:,i)=labeled_img;
end
all_stats=all_stats';

    
figure();
for i = 1:nframes
    imshow(detected_movie(:,:,i),[]);
    pause(0.03); % Adjust speed
end


%% Choose an ROI
image_stack_processed = detected_movie;
run = 1;
nslot=nframes;
current_frame = detected_movie(:,:,1);

%% Tracking parameters (must be consistent)
max_speed_per_frame = 2;      % pixels/frame (consistent with 5px max_step)
cost_cutoff = 2;              % Maximum allowed matching cost
min_traj_length = 5;          % Minimum frames for valid trajectory
max_gap = 5;                  % Maximum allowed gap (frames) to miss
max_step_distance = 3;        % Maximum distance between successive points
alpha = 0.7;                  % Position vs. features weight

%% Build trajectories and populate all_intensity_tracked and all_area_tracked

% Initialize mapping from particles to track IDs
particle_to_track_map = cell(nframes, 1);
next_track_id = 1;

% Extract particle features
% Initialize cell arrays for each property
centroids = cell(nframes, 1);
intensity_details = cell(nframes, 1);
area_details = cell(nframes, 1);
circularity_details = cell(nframes, 1);
equivdiameter_details = cell(nframes, 1);
boundingbox_details = cell(nframes, 1);
eulernumber_details = cell(nframes, 1);
minoraxislength_details = cell(nframes, 1);
extent_details = cell(nframes, 1);
orientation_details = cell(nframes, 1);
perimeter_details = cell(nframes, 1);
convexarea_details = cell(nframes, 1);
filledarea_details = cell(nframes, 1);
pixelidxlist_details = cell(nframes, 1);
convexhull_details = cell(nframes, 1);
filledimage_details = cell(nframes, 1);
pixellist_details = cell(nframes, 1);
conveximage_details = cell(nframes, 1);
image_details = cell(nframes, 1);
solidity_details = cell(nframes, 1);
eccentricity_details = cell(nframes, 1);
maxferetproperties_details = cell(nframes, 1);
subarrayidx_details = cell(nframes, 1);

for frame = 1:nframes
    if ~isempty(all_stats{frame})
        centroids{frame} = vertcat(all_stats{frame}.Centroid);
        intensity_details{frame} = vertcat(all_stats{frame}.Intensity);
        area_details{frame} = vertcat(all_stats{frame}.Area);
        circularity_details{frame} = vertcat(all_stats{frame}.Circularity);
        equivdiameter_details{frame} = vertcat(all_stats{frame}.EquivDiameter);
        boundingbox_details{frame} = vertcat(all_stats{frame}.BoundingBox);
        eulernumber_details{frame} = vertcat(all_stats{frame}.EulerNumber);
        minoraxislength_details{frame} = vertcat(all_stats{frame}.MinorAxisLength);
        extent_details{frame} = vertcat(all_stats{frame}.Extent);
        orientation_details{frame} = vertcat(all_stats{frame}.Orientation);
        perimeter_details{frame} = vertcat(all_stats{frame}.Perimeter);
        convexarea_details{frame} = vertcat(all_stats{frame}.ConvexArea);
        filledarea_details{frame} = vertcat(all_stats{frame}.FilledArea);
        pixelidxlist_details{frame} = {all_stats{frame}.PixelIdxList}';
        convexhull_details{frame} = {all_stats{frame}.ConvexHull}';
        filledimage_details{frame} = {all_stats{frame}.FilledImage}';
        pixellist_details{frame} = {all_stats{frame}.PixelList}';
        conveximage_details{frame} = {all_stats{frame}.ConvexImage}';
        image_details{frame} = {all_stats{frame}.Image}';
        solidity_details{frame} = vertcat(all_stats{frame}.Solidity);
        eccentricity_details{frame} = vertcat(all_stats{frame}.Eccentricity);
        majoraxislength_details{frame} = vertcat(all_stats{frame}.MajorAxisLength);
        subarrayidx_details{frame} = {all_stats{frame}.SubarrayIdx}';
    else
        centroids{frame} = [];
        intensity_details{frame} = [];
        area_details{frame} = [];
        circularity_details{frame} = [];
        equivdiameter_details{frame} = [];
        boundingbox_details{frame} = [];
        eulernumber_details{frame} = [];
        minoraxislength_details{frame} = [];
        extent_details{frame} = [];
        orientation_details{frame} = [];
        perimeter_details{frame} = [];
        convexarea_details{frame} = [];
        filledarea_details{frame} = [];
        pixelidxlist_details{frame} = [];
        convexhull_details{frame} = [];
        filledimage_details{frame} = [];
        pixellist_details{frame} = [];
        conveximage_details{frame} = [];
        image_details{frame} = [];
        solidity_details{frame} = [];
        eccentricity_details{frame} = [];
        majoraxislength_details{frame} = [];
        subarrayidx_details{frame} = [];
    end
end
% Frame-to-frame matching using the linear assignment algorithm munkres(standard program)
% The assignment of the munkres for each frame is matched with sucessive
% frame to form a matched_pairs variable of size nframes-1

matched_pairs = cell(nframes-1, 1);
for t = 1:nframes-1
    pts_t = centroids{t};
    pts_t1 = centroids{t+1};
    if isempty(pts_t) || isempty(pts_t1)
        continue;
    end
    
    % Create cost matrix and perform matching
    cost_matrix = enhanced_cost_matrix(all_stats{t}, all_stats{t+1}, max_speed_per_frame, alpha);
    [assignments, ~] = munkres(cost_matrix);
    
    % Validate matches
    valid_assignments = zeros(size(assignments));
    for k = 1:length(assignments)
        if assignments(k) > 0
            distance = norm(pts_t(k,:) - pts_t1(assignments(k),:));
            if cost_matrix(k, assignments(k)) <= cost_cutoff && distance <= max_speed_per_frame
                valid_assignments(k) = assignments(k);
            end
        end
    end
    matched_pairs{t} = valid_assignments;
end

% Initialize tracks in first frame
if ~isempty(all_stats{1})
    num_particles = length(all_stats{1});
    particle_to_track_map{1} = (next_track_id:next_track_id+num_particles-1)';
    next_track_id = next_track_id + num_particles;
end

% The track is built using the matched pairs variable, which stores the
% assigned particle for the next frame

for t = 1:nframes-1
    if ~isempty(matched_pairs{t})
        num_particles_next = length(all_stats{t+1});
        particle_to_track_map{t+1} = zeros(num_particles_next, 1);
        
        % Assign existing track IDs to matched particles
        for i = 1:length(matched_pairs{t})
            j = matched_pairs{t}(i);
            if j > 0 && particle_to_track_map{t}(i) > 0
                particle_to_track_map{t+1}(j) = particle_to_track_map{t}(i);
            end
        end
        
        % Assign new IDs to unmatched particles where there is no match in
        % sucessive frame and the particle is assigned 0 a new track is
        % generated
        unmatched = find(particle_to_track_map{t+1} == 0);
        particle_to_track_map{t+1}(unmatched) = next_track_id:next_track_id+length(unmatched)-1;
        next_track_id = next_track_id + length(unmatched);
    end
end
% Total number of tracks possible
num_tracks = next_track_id - 1;

%% Enhanced trajectory building with intensity and area interpolation

% Initialize trajectory storage
% For tracking regions across multiple frames
trajectories = struct(...
    'trajectory_id', {}, ...      % Unique ID for each tracked region
    'particle_id', {}, ...         % Changed from 'region_id'
    'frames', {}, ...              % Frame numbers where region appears
    'positions', {}, ...           % Changed from 'centroids'
    'intensity', {}, ...           % Intensity values at each frame
    'area', {}, ...                % Area at each frame
    'circularity', {}, ...         % Circularity at each frame
    'equivdiameter', {}, ...       % Equivalent diameter at each frame
    'boundingbox', {}, ...         % Bounding boxes at each frame (Nx4 array)
    'perimeter', {}, ...           % Perimeter at each frame
    'eccentricity', {}, ...        % Eccentricity at each frame
    'solidity', {}, ...            % Solidity at each frame
    'extent', {}, ...              % Extent at each frame
    'orientation', {}, ...         % Orientation at each frame
    'majoraxislength', {}, ...     % Major axis length at each frame
    'minoraxislength', {}, ...     % Minor axis length at each frame
    'convexarea', {}, ...          % Convex area at each frame
    'filledarea', {}, ...          % Filled area at each frame
    'eulernumber', {}, ...         % Euler number at each frame
    'convexhull', {}, ...          % Convex hull at each frame
    'filledimage', {}, ...         % Filled image at each frame
    'conveximage', {}, ...         % Convex image at each frame
    'image', {}, ...               % Region image at each frame
    'pixelidxlist', {}, ...        % Pixel indices at each frame
    'pixellist', {}, ...           % Pixel coordinates at each frame
    'subarrayidx', {}, ...         % Subarray indices at each frame
    'length', {}, ...              % Number of frames region appears
    'max_gap', {}, ...             % Maximum gap in tracking
    'valid_frames', {});           % Valid frame indices

% Build complete trajectories with interpolated values of positions,area
% and intensity kind of a guess to track accurately 

for track_id = 1:num_tracks
    % Find all frames where this object appears across the movie
    valid_frames = [];
    particle_indices = [];
    for t = 1:nframes
        if ~isempty(particle_to_track_map{t})
            idx = find(particle_to_track_map{t} == track_id);
            if ~isempty(idx)
                valid_frames(end+1) = t;
                particle_indices(end+1) = idx;
            end
        end
    end
    
    % Skip if track is too short
    if length(valid_frames) < min_traj_length
        continue;
    end
    
    % Extract actual measurements of all properties for valid frames
    measured_centroids = zeros(length(valid_frames), 2);
    measured_intensity = zeros(1, length(valid_frames));
    measured_area = zeros(1, length(valid_frames));
    measured_circularity = zeros(1, length(valid_frames));
    measured_equivdiameter = zeros(1, length(valid_frames));
    measured_boundingbox = zeros(length(valid_frames), 4);
    measured_eulernumber = zeros(1, length(valid_frames));
    measured_minoraxislength = zeros(1, length(valid_frames));
    measured_extent = zeros(1, length(valid_frames));
    measured_orientation = zeros(1, length(valid_frames));
    measured_perimeter = zeros(1, length(valid_frames));
    measured_convexarea = zeros(1, length(valid_frames));
    measured_filledarea = zeros(1, length(valid_frames));
    measured_solidity = zeros(1, length(valid_frames));
    measured_eccentricity = zeros(1, length(valid_frames));
    measured_majoraxislength = zeros(1, length(valid_frames));
    
    % For properties that are arrays or cells
    measured_convexhull = cell(1, length(valid_frames));
    measured_filledimage = cell(1, length(valid_frames));
    measured_pixellist = cell(1, length(valid_frames));
    measured_conveximage = cell(1, length(valid_frames));
    measured_image = cell(1, length(valid_frames));
    measured_pixelidxlist = cell(1, length(valid_frames));
    measured_subarrayidx = cell(1, length(valid_frames));
    
    for i = 1:length(valid_frames)
        t = valid_frames(i);
        idx = particle_indices(i);
        
        % Numeric properties
        measured_centroids(i,:) = centroids{t}(idx,:);
        measured_intensity(i) = intensity_details{t}(idx);
        measured_area(i) = area_details{t}(idx);
        measured_circularity(i) = circularity_details{t}(idx);
        measured_equivdiameter(i) = equivdiameter_details{t}(idx);
        measured_boundingbox(i,:) = boundingbox_details{t}(idx,:);
        measured_eulernumber(i) = eulernumber_details{t}(idx);
        measured_minoraxislength(i) = minoraxislength_details{t}(idx);
        measured_extent(i) = extent_details{t}(idx);
        measured_orientation(i) = orientation_details{t}(idx);
        measured_perimeter(i) = perimeter_details{t}(idx);
        measured_convexarea(i) = convexarea_details{t}(idx);
        measured_filledarea(i) = filledarea_details{t}(idx);
        measured_solidity(i) = solidity_details{t}(idx);
        measured_eccentricity(i) = eccentricity_details{t}(idx);
        measured_majoraxislength(i) = majoraxislength_details{t}(idx);
        
        % Cell/array properties
        measured_convexhull{i} = convexhull_details{t}{idx};
        measured_filledimage{i} = filledimage_details{t}{idx};
        measured_pixellist{i} = pixellist_details{t}{idx};
        measured_conveximage{i} = conveximage_details{t}{idx};
        measured_image{i} = image_details{t}{idx};
        measured_pixelidxlist{i} = pixelidxlist_details{t}{idx};
        measured_subarrayidx{i} = subarrayidx_details{t}{idx};
    end
    
    % Create full timeline from first to last appearance
    full_frames = min(valid_frames):max(valid_frames);
    full_centroids = NaN(length(full_frames), 2);
    full_intensity = NaN(1, length(full_frames));
    full_area = NaN(1, length(full_frames));
    full_circularity = NaN(1, length(full_frames));
    full_equivdiameter = NaN(1, length(full_frames));
    full_boundingbox = NaN(length(full_frames), 4);
    full_eulernumber = NaN(1, length(full_frames));
    full_minoraxislength = NaN(1, length(full_frames));
    full_extent = NaN(1, length(full_frames));
    full_orientation = NaN(1, length(full_frames));
    full_perimeter = NaN(1, length(full_frames));
    full_convexarea = NaN(1, length(full_frames));
    full_filledarea = NaN(1, length(full_frames));
    full_solidity = NaN(1, length(full_frames));
    full_eccentricity = NaN(1, length(full_frames));
    full_majoraxislength = NaN(1, length(full_frames));
    
    % Cell arrays for full timeline
    full_convexhull = cell(1, length(full_frames));
    full_filledimage = cell(1, length(full_frames));
    full_pixellist = cell(1, length(full_frames));
    full_conveximage = cell(1, length(full_frames));
    full_image = cell(1, length(full_frames));
    full_pixelidxlist = cell(1, length(full_frames));
    full_subarrayidx = cell(1, length(full_frames));

    
    % Fill in measured values
    [~, loc] = ismember(valid_frames, full_frames);
    full_centroids(loc,:) = measured_centroids;
    full_intensity(loc) = measured_intensity;
    full_area(loc) = measured_area;
    full_circularity(loc) = measured_circularity;
    full_equivdiameter(loc) = measured_equivdiameter;
    full_boundingbox(loc,:) = measured_boundingbox;
    full_eulernumber(loc) = measured_eulernumber;
    full_minoraxislength(loc) = measured_minoraxislength;
    full_extent(loc) = measured_extent;
    full_orientation(loc) = measured_orientation;
    full_perimeter(loc) = measured_perimeter;
    full_convexarea(loc) = measured_convexarea;
    full_filledarea(loc) = measured_filledarea;
    full_solidity(loc) = measured_solidity;
    full_eccentricity(loc) = measured_eccentricity;
    full_majoraxislength(loc) = measured_majoraxislength;
    
    % Fill cell arrays
    full_convexhull(loc) = measured_convexhull;
    full_filledimage(loc) = measured_filledimage;
    full_pixellist(loc) = measured_pixellist;
    full_conveximage(loc) = measured_conveximage;
    full_image(loc) = measured_image;
    full_pixelidxlist(loc) = measured_pixelidxlist;
    full_subarrayidx(loc) = measured_subarrayidx;
    
    % Interpolate missing values for numeric properties
    for dim = 1:2
        full_centroids(:,dim) = interp1(full_frames(loc), measured_centroids(:,dim), full_frames, 'linear');
    end
    full_intensity = interp1(valid_frames, measured_intensity, full_frames, 'nearest');
    full_area = interp1(valid_frames, measured_area, full_frames, 'nearest');
    full_circularity = interp1(valid_frames, measured_circularity, full_frames, 'linear');
    full_equivdiameter = interp1(valid_frames, measured_equivdiameter, full_frames, 'linear');
    full_eulernumber = interp1(valid_frames, measured_eulernumber, full_frames, 'nearest');
    full_minoraxislength = interp1(valid_frames, measured_minoraxislength, full_frames, 'linear');
    full_extent = interp1(valid_frames, measured_extent, full_frames, 'linear');
    full_orientation = interp1(valid_frames, measured_orientation, full_frames, 'linear');
    full_perimeter = interp1(valid_frames, measured_perimeter, full_frames, 'linear');
    full_convexarea = interp1(valid_frames, measured_convexarea, full_frames, 'linear');
    full_filledarea = interp1(valid_frames, measured_filledarea, full_frames, 'linear');
    full_solidity = interp1(valid_frames, measured_solidity, full_frames, 'linear');
    full_eccentricity = interp1(valid_frames, measured_eccentricity, full_frames, 'linear');
    full_majoraxislength = interp1(valid_frames, measured_majoraxislength, full_frames, 'linear');
    
    % For bounding box, interpolate each dimension separately
    for dim = 1:4
        full_boundingbox(:,dim) = interp1(full_frames(loc), measured_boundingbox(:,dim), full_frames, 'linear');
    end
    
    % Calculate gap statistics
    gap_lengths = diff(valid_frames) - 1;
    max_gap = max([0, gap_lengths(gap_lengths > 0)]);
    
    % Store trajectory with ALL properties
    trajectories(end+1) = struct(...
        'trajectory_id', length(trajectories)+1, ...
        'particle_id', track_id, ...
        'frames', full_frames, ...
        'positions', full_centroids, ...
        'intensity', full_intensity, ...
        'area', full_area, ...
        'circularity', full_circularity, ...
        'equivdiameter', full_equivdiameter, ...
        'boundingbox', full_boundingbox, ...
        'eulernumber', full_eulernumber, ...
        'minoraxislength', full_minoraxislength, ...
        'extent', full_extent, ...
        'orientation', full_orientation, ...
        'perimeter', full_perimeter, ...
        'convexarea', full_convexarea, ...
        'filledarea', full_filledarea, ...
        'solidity', full_solidity, ...
        'eccentricity', full_eccentricity, ...
        'majoraxislength', full_majoraxislength, ...
        'convexhull', {full_convexhull}, ...
        'filledimage', {full_filledimage}, ...
        'pixellist', {full_pixellist}, ...
        'conveximage', {full_conveximage}, ...
        'image', {full_image}, ...
        'pixelidxlist', {full_pixelidxlist}, ...
        'subarrayidx', {full_subarrayidx}, ...
        'length', length(full_frames), ...
        'max_gap', max_gap, ...
        'valid_frames', valid_frames);
end

%% Filter trajectories by length and gap size
valid_trajectories = trajectories([trajectories.length] >= min_traj_length & ...
                                 [trajectories.max_gap] <= max_gap);
%% Parameters for merging trajectories
max_merge_gap = 5;          % Maximum allowed frame gap between trajectories
max_merge_distance = 15;     % Maximum spatial distance for merging
min_merged_length = 75;     % Minimum length after merging to keep

% Function to calculate distance between trajectory endpoints
endpoint_distance = @(t1,t2) min([
    norm(t1.positions(end,:) - t2.positions(1,:)),   % end of t1 to start of t2
    norm(t1.positions(1,:) - t2.positions(end,:))    % start of t1 to end of t2
]);

% Initialize merge tracking
merge_iteration = 0;
max_iterations = 5;       % Safety limit to prevent infinite loops

% Initialize new_trajectories as a structure array with the same fields
if ~isempty(trajectories)
    new_trajectories = repmat( struct(...
    'trajectory_id', length(trajectories)+1, ...
    'particle_id', track_id, ...
    'frames', full_frames, ...
    'positions', full_centroids, ...
    'intensity', full_intensity, ...
    'area', full_area, ...
    'circularity', full_circularity, ...
    'equivdiameter', full_equivdiameter, ...
    'boundingbox', full_boundingbox, ...
    'eulernumber', full_eulernumber, ...
    'minoraxislength', full_minoraxislength, ...
    'extent', full_extent, ...
    'orientation', full_orientation, ...
    'perimeter', full_perimeter, ...
    'convexarea', full_convexarea, ...
    'filledarea', full_filledarea, ...
    'solidity', full_solidity, ...
    'eccentricity', full_eccentricity, ...
    'majoraxislength', full_majoraxislength, ...
    'convexhull', {full_convexhull}, ...
    'filledimage', {full_filledimage}, ...
    'pixellist', {full_pixellist}, ...
    'conveximage', {full_conveximage}, ...
    'image', {full_image}, ...
    'pixelidxlist', {full_pixelidxlist}, ...
    'subarrayidx', {full_subarrayidx}, ...
    'length', length(full_frames), ...
    'max_gap', max_gap, ...
    'valid_frames', valid_frames), 0, 1);
else
    new_trajectories = struct.empty(0,1);
end

% Main merging loop
while merge_iteration < max_iterations
    merge_iteration = merge_iteration + 1;
    fprintf('\n=== Merge Iteration %d ===\n', merge_iteration);
    
    % Initialize variables for this iteration
    n_trajectories = length(trajectories);
    merged_indices = false(1, n_trajectories);
    new_trajectories = repmat( struct(...
    'trajectory_id', length(trajectories)+1, ...
    'particle_id', track_id, ...
    'frames', full_frames, ...
    'positions', full_centroids, ...
    'intensity', full_intensity, ...
    'area', full_area, ...
    'circularity', full_circularity, ...
    'equivdiameter', full_equivdiameter, ...
    'boundingbox', full_boundingbox, ...
    'eulernumber', full_eulernumber, ...
    'minoraxislength', full_minoraxislength, ...
    'extent', full_extent, ...
    'orientation', full_orientation, ...
    'perimeter', full_perimeter, ...
    'convexarea', full_convexarea, ...
    'filledarea', full_filledarea, ...
    'solidity', full_solidity, ...
    'eccentricity', full_eccentricity, ...
    'majoraxislength', full_majoraxislength, ...
    'convexhull', {full_convexhull}, ...
    'filledimage', {full_filledimage}, ...
    'pixellist', {full_pixellist}, ...
    'conveximage', {full_conveximage}, ...
    'image', {full_image}, ...
    'pixelidxlist', {full_pixelidxlist}, ...
    'subarrayidx', {full_subarrayidx}, ...
    'length', length(full_frames), ...
    'max_gap', max_gap, ...
    'valid_frames', valid_frames), 0, 1);
    merge_count = 0;
    
    % Try to merge all possible trajectory pairs
    for i = 1:n_trajectories
        if merged_indices(i), continue; end
        
        best_merge = [];
        best_j = 0;
        best_dist = inf;
        best_gap = inf;
        
        % Find the best candidate to merge with current trajectory
        for j = 1:n_trajectories
            if i == j || merged_indices(j), continue; end
            
            % Calculate distance and gap between trajectories
            dist = endpoint_distance(trajectories(i), trajectories(j));
            gap = min([
                trajectories(j).frames(1) - trajectories(i).frames(end), 
                trajectories(i).frames(1) - trajectories(j).frames(end)
            ]);
            
            % Check merge conditions
            if gap > 0 && gap <= max_merge_gap && dist <= max_merge_distance
                disp('There is gap');
                if dist < best_dist || (dist == best_dist && gap < best_gap)
                    best_dist = dist;
                    best_gap = gap;
                    best_j = j;
                end
            end
        end
        
        % Perform the merge if a candidate was found
        if best_j > 0
            % Determine merge direction
            if trajectories(best_j).frames(1) > trajectories(i).frames(end)
                % t1 then t2
                new_pos = [trajectories(i).positions; trajectories(best_j).positions];
                new_frames = [trajectories(i).frames, trajectories(best_j).frames];
                new_intensity = [trajectories(i).intensity, trajectories(best_j).intensity];
                new_area = [trajectories(i).area, trajectories(best_j).area];
                new_valid_frames = [trajectories(i).valid_frames, trajectories(best_j).valid_frames];
            else
                % t2 then t1
                new_pos = [trajectories(best_j).positions; trajectories(i).positions];
                new_frames = [trajectories(best_j).frames, trajectories(i).frames];
                new_intensity = [trajectories(best_j).intensity, trajectories(i).intensity];
                new_area = [trajectories(best_j).area, trajectories(i).area];
                new_valid_frames = [trajectories(best_j).valid_frames, trajectories(i).valid_frames];
            end
            
            % Create merged trajectory
merged_traj = struct(...
    'trajectory_id', length(new_trajectories)+1, ...
    'particle_id', max([trajectories(i).particle_id, trajectories(best_j).particle_id]), ...
    'frames', new_frames, ...
    'positions', new_pos, ...
    'intensity', new_intensity, ...
    'area', new_area, ...
    'circularity', [trajectories(i).circularity, trajectories(best_j).circularity], ...
    'equivdiameter', [trajectories(i).equivdiameter, trajectories(best_j).equivdiameter], ...
    'boundingbox', [trajectories(i).boundingbox; trajectories(best_j).boundingbox], ...
    'eulernumber', [trajectories(i).eulernumber, trajectories(best_j).eulernumber], ...
    'minoraxislength', [trajectories(i).minoraxislength, trajectories(best_j).minoraxislength], ...
    'extent', [trajectories(i).extent, trajectories(best_j).extent], ...
    'orientation', [trajectories(i).orientation, trajectories(best_j).orientation], ...
    'perimeter', [trajectories(i).perimeter, trajectories(best_j).perimeter], ...
    'convexarea', [trajectories(i).convexarea, trajectories(best_j).convexarea], ...
    'filledarea', [trajectories(i).filledarea, trajectories(best_j).filledarea], ...
    'solidity', [trajectories(i).solidity, trajectories(best_j).solidity], ...
    'eccentricity', [trajectories(i).eccentricity, trajectories(best_j).eccentricity], ...
    'majoraxislength', [trajectories(i).majoraxislength, trajectories(best_j).majoraxislength], ...
    'convexhull', [trajectories(i).convexhull, trajectories(best_j).convexhull], ...
    'filledimage', [trajectories(i).filledimage, trajectories(best_j).filledimage], ...
    'pixellist', [trajectories(i).pixellist, trajectories(best_j).pixellist], ...
    'conveximage', [trajectories(i).conveximage, trajectories(best_j).conveximage], ...
    'image', [trajectories(i).image, trajectories(best_j).image], ...
    'pixelidxlist', [trajectories(i).pixelidxlist, trajectories(best_j).pixelidxlist], ...
    'subarrayidx', [trajectories(i).subarrayidx, trajectories(best_j).subarrayidx], ...
    'length', length(new_frames), ...
    'max_gap', max([trajectories(i).max_gap, trajectories(best_j).max_gap, gap-1]), ...
    'valid_frames', new_valid_frames);
            
            fprintf('Merged Trajectories %d and %d (dist=%.2f, gap=%d)\n',...
                i, best_j, best_dist, best_gap);
            
            % Add to new trajectories
            new_trajectories(end+1) = merged_traj;
            merged_indices([i, best_j]) = true;
            merge_count = merge_count + 1;
        else
            % Keep unmerged trajectory (with updated ID)
            trajectories(i).trajectory_id =  length(new_trajectories)+1;
            new_trajectories(end+1) = trajectories(i);
            merged_indices(i) = true;
        end
    end
    
    % Update trajectories for next iteration
    trajectories = new_trajectories;
    
    % Stop if no more merges were performed
    if merge_count == 0
        fprintf('No more merges possible. Ending merge loop.\n');
        break;
    end
end

% Final filtering by minimum length
trajectories = trajectories([trajectories.length] >= min_merged_length);

% Re-number trajectories
for i = 1:length(trajectories)
    trajectories(i).trajectory_id = i;
end

% Diagnostic output
fprintf('\n=== Final Trajectory Statistics ===\n');
fprintf('Total trajectories after merging: %d\n', length(trajectories));
if ~isempty(trajectories)
    lengths = [trajectories.length];
    max_gaps = [trajectories.max_gap];
    
    fprintf('Length stats: Min=%d, Mean=%.1f, Max=%d\n', ...
        min(lengths), mean(lengths), max(lengths));
    fprintf('Max gap stats: Min=%d, Mean=%.1f, Max=%d\n', ...
        min(max_gaps), mean(max_gaps), max(max_gaps));
end

[drift_corr_traj] = drift_removal(trajectories,nframes);

%% Enhanced Visualization
figure;
imshow(image_stack(:,:,1), []);
hold on;
title('Trajectories with Unique Colors');
xlabel('X position (px)');
ylabel('Y position (px)');

% Create a colormap with enough distinct colors
num_trajectories = length(drift_corr_traj);
colors = lines(num_trajectories);  % 'lines' is a MATLAB colormap with distinct colors

for i = 1:num_trajectories
    traj = drift_corr_traj(i);
    
    % Plot trajectory with unique color
    h = plot(traj.positions(:,1), traj.positions(:,2), ...
        'Color', colors(i,:), 'LineWidth', 1);
    
    % Plot start/end markers
    plot(traj.positions(1,1), traj.positions(1,2), 'o', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
    plot(traj.positions(end,1), traj.positions(end,2), 's', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    % Store handle and label for legend (optional)
    legend_handles(i) = h;
    legend_labels{i} = ['Trajectory ' num2str(i)];
end


legend(legend_handles, legend_labels, 'Location', 'best');

grid on;
filename=sprintf('All_Trajectories.fig');
savefig(fullfile(directory,filename));

figure;
imshow(image_stack(:,:,1), []);
hold on;
title('Trajectories with Unique Colors');
xlabel('X position (px)');
ylabel('Y position (px)');

% Create a colormap with enough distinct colors
num_trajectories = length(drift_corr_traj);
colors = lines(num_trajectories);  % 'lines' is a MATLAB colormap with distinct colors

for i = 1:num_trajectories
    traj = trajectories(i);
    
    % Plot trajectory with unique color
    h = plot(traj.positions(:,1), traj.positions(:,2), ...
        'Color', colors(i,:), 'LineWidth', 1);
    
    % Plot start/end markers
    plot(traj.positions(1,1), traj.positions(1,2), 'o', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
    plot(traj.positions(end,1), traj.positions(end,2), 's', ...
        'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
    % Store handle and label for legend (optional)
    legend_handles(i) = h;
    legend_labels{i} = ['Trajectory ' num2str(i)];
end


legend(legend_handles, legend_labels, 'Location', 'best');

grid on;
% filename=sprintf('Non_Drift_Corrected_Trajectories.fig');
% savefig(fullfile(directory,filename));
close;

%% GET ALL THE OBJECTS AND GENERATE KYMOGRAPH AND STICH THEM TOGETHER

x_start = 1;
y_start = 1;
x_end = size(detected_movie, 2);
y_end = size(detected_movie, 1);
    
sizest=size(all_stats{1,1});
np=sizest(1);
line_scan_length = 6;
 


num_particles_detected=zeros(1,nframes);
max_val=0;
max_index=1;
for i =1 :nframes
    num_particles_detected(i)=size(all_stats{i}, 1);
    S=all_stats{i};
    if length(all_stats{i})>max_val
        max_val = length(all_stats{i});
        max_index=i;
    end    
end


% Extract XC, YC, and Intensity 
x_coords = int16([all_stats{max_index}.XC]);
y_coords = int16([all_stats{max_index}.YC]);
Intparticle = [all_stats{max_index}.Intensity];
 
%% Following a single detected object

max_particles=size(trajectories,2);



num_trajectories = length(trajectories);
max_nframes = 0;
if num_trajectories > 0
    for traj_idx = 1:num_trajectories
        if ~isempty(trajectories(traj_idx).positions)
            current_nframes = size(trajectories(traj_idx).positions, 1);
            if current_nframes > max_nframes
                max_nframes = current_nframes;
            end
        end
    end
else
    max_nframes = 0;
end

X_Centroids = NaN(num_trajectories, max_nframes);
Y_Centroids = NaN(num_trajectories, max_nframes);

for traj_idx = 1:num_trajectories
    % Check if the trajectory has position data.
    if ~isempty(trajectories(traj_idx).positions)
        % Get the length of the current trajectory's data.
        current_traj_nframes = size(trajectories(traj_idx).positions, 1);
        
        % Extract the X coordinates (first column) and assign to the correct row.
        % The assignment only populates up to the length of the current trajectory.
        X_Centroids(traj_idx, 1:current_traj_nframes) = trajectories(traj_idx).positions(:, 1)';
        
        % Extract the Y coordinates (second column) and assign to the correct row.
        Y_Centroids(traj_idx, 1:current_traj_nframes) = trajectories(traj_idx).positions(:, 2)';
    end
end

total_max_movement = 10;
max_movement_allowed_per_frame = 1;
params = struct();
params.MinArea = Min_Area_Thresh;       
params.MaxArea = Max_Area_Thresh;      
params.FrameStart = 1;       
params.FrameEnd = nframes;       


merge_stats = merge_tracker(all_stats,params); 
all_objects = track_object_lifetimes(merge_stats,nframes,params);
minimum_particle_size=6;

wrong_indices = find(num_particles_detected>max_particles);
if length(wrong_indices)>0
    num_particles_detected(wrong_indices)=max_particles;
end
[kymograph_position,kymograph_intensity,particle_kymograph,dist_from_centroid] = ...
    generate_kymograph(labeled_movie_norm,binary_stack,x_coords,y_coords,all_stats,nframes,minimum_particle_size,...
    max_particles,num_particles_detected,current_frame,X_Centroids,Y_Centroids);

clear y kymograph_stitched kymograph_stitched_intensity;

[kymograph_stitched,kymograph_stitched_intensity,total_lifetime_normalized] = ...
    get_lifetime(num_particles_detected,kymograph_position,kymograph_intensity,nframes,max_particles);

%%
[sRg_x,sRg_y,R_g_x,R_g_y,r_x,r_y] = rg_trajectory(trajectories);
idx_mobile = find(sRg_x>filter_mobile);
idx_fixed = find(sRg_x<filter_mobile);

% 1. Extract the mobile trajectories
all_trajectory_mobile = trajectories(idx_mobile);

% 2. Extract the 9fixed trajectories
all_trajectory_fixed = trajectories(idx_fixed);

% Optional: Print a summary to the console
fprintf('Segregation Complete:\n');
fprintf(' - Mobile Clusters: %d\n', length(all_trajectory_mobile));
fprintf(' - Fixed Clusters:  %d\n', length(all_trajectory_fixed));


%% PLOTS 


figure(2); 
imshow(labeled_movie(:,:,max_index), []); 
colormap jet; 
colorbar; 
title("Labelled endosomes");
filename=sprintf('Labeled Objects.fig');
savefig(fullfile(directory,filename));
close;

figure(3);
imshow(kymograph_stitched);
colormap jet;
title("DETECTION KYMOGRAPH X");
image_data = getframe(gcf).cdata;
filename=sprintf('Detection Kymograph.fig');
savefig(fullfile(directory,filename));
close;

% Display intensity kymograph
figure(8);
imshow(kymograph_stitched_intensity, [0 200]);
colormap jet;
colorbar;
title("Intensity Kymograph");
filename=sprintf('Intensity Kymograph.fig');
savefig(fullfile(directory,filename));
close;



cd(directory);
gifFile = 'Detection Movie.gif';

% Check if the GIF file exists and its format
if exist(gifFile, 'file')
    % Check if the GIF file is in GIF89a format
    fid = fopen(gifFile, 'r');
    format = fread(fid, 6, '*char')';
    fclose(fid);
    
    % If it's not GIF89a, delete the file to create a new GIF
    if ~strcmp(format, 'GIF89a')
        warning('GIF file is not in GIF89a format. Creating a new file.');
        delete(gifFile);
    end
end

% Write the first frame (or create a new GIF89a file)
imwrite(detected_movie(:,:,1), gifFile, 'LoopCount', Inf, 'DelayTime', 0.1);

% Append subsequent frames
for i = 2:size(detected_movie, 3)
    imwrite(detected_movie(:,:,i), gifFile, 'WriteMode', 'append', 'DelayTime', 0.1);
end

%% SAVING KEY VARIABLES 


cd(directory);
% Save results
results_file_name = 'all_information';
save(results_file_name, 'kymograph_stitched_intensity', 'kymograph_stitched','trajectories','merge_stats','all_objects','all_stats','directory','all_trajectory_fixed','all_trajectory_mobile');
cd(olddir);


function [assignment,cost] = munkres(costMat)

assignment = zeros(1,size(costMat,1));
cost = 0;
validMat = costMat == costMat & costMat < Inf;
bigM = 10^(ceil(log10(sum(costMat(validMat))))+1);
costMat(~validMat) = bigM;
% costMat(costMat~=costMat)=Inf;
% validMat = costMat<Inf;
validCol = any(validMat,1);
validRow = any(validMat,2);
nRows = sum(validRow);
nCols = sum(validCol);
n = max(nRows,nCols);
if ~n
    return
end
maxv=10*max(costMat(validMat));
dMat = zeros(n) + maxv;
dMat(1:nRows,1:nCols) = costMat(validRow,validCol);
%*************************************************
% Munkres' Assignment Algorithm starts here
%*************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STEP 1: Subtract the row minimum from each row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minR = min(dMat,[],2);
minC = min(bsxfun(@minus, dMat, minR));
%**************************************************************************  
%   STEP 2: Find a zero of dMat. If there are no starred zeros in its
%           column or row start the zero. Repeat for each zero
%**************************************************************************
zP = dMat == bsxfun(@plus, minC, minR);
starZ = zeros(n,1);
while any(zP(:))
    [r,c]=find(zP,1);
    starZ(r)=c;
    zP(r,:)=false;
    zP(:,c)=false;
end
while 1
%**************************************************************************
%   STEP 3: Cover each column with a starred zero. If all the columns are
%           covered then the matching is maximum
%**************************************************************************
    if all(starZ>0)
        break
    end
    coverColumn = false(1,n);
    coverColumn(starZ(starZ>0))=true;
    coverRow = false(n,1);
    primeZ = zeros(n,1);
    [rIdx, cIdx] = find(dMat(~coverRow,~coverColumn)==bsxfun(@plus,minR(~coverRow),minC(~coverColumn)));
    while 1
        %**************************************************************************
        %   STEP 4: Find a noncovered zero and prime it.  If there is no starred
        %           zero in the row containing this primed zero, Go to Step 5.  
        %           Otherwise, cover this row and uncover the column containing 
        %           the starred zero. Continue in this manner until there are no 
        %           uncovered zeros left. Save the smallest uncovered value and 
        %           Go to Step 6.
        %**************************************************************************
        cR = find(~coverRow);
        cC = find(~coverColumn);
        rIdx = cR(rIdx);
        cIdx = cC(cIdx);
        Step = 6;
        while ~isempty(cIdx)
            uZr = rIdx(1);
            uZc = cIdx(1);
            primeZ(uZr) = uZc;
            stz = starZ(uZr);
            if ~stz
                Step = 5;
                break;
            end
            coverRow(uZr) = true;
            coverColumn(stz) = false;
            z = rIdx==uZr;
            rIdx(z) = [];
            cIdx(z) = [];
            cR = find(~coverRow);
            z = dMat(~coverRow,stz) == minR(~coverRow) + minC(stz);
            rIdx = [rIdx(:);cR(z)];
            cIdx = [cIdx(:);stz(ones(sum(z),1))];
        end
        if Step == 6
            % *************************************************************************
            % STEP 6: Add the minimum uncovered value to every element of each covered
            %         row, and subtract it from every element of each uncovered column.
            %         Return to Step 4 without altering any stars, primes, or covered lines.
            %**************************************************************************
            [minval,rIdx,cIdx]=outerplus(dMat(~coverRow,~coverColumn),minR(~coverRow),minC(~coverColumn));            
            minC(~coverColumn) = minC(~coverColumn) + minval;
            minR(coverRow) = minR(coverRow) - minval;
        else
            break
        end
    end
    %**************************************************************************
    % STEP 5:
    %  Construct a series of alternating primed and starred zeros as
    %  follows:
    %  Let Z0 represent the uncovered primed zero found in Step 4.
    %  Let Z1 denote the starred zero in the column of Z0 (if any).
    %  Let Z2 denote the primed zero in the row of Z1 (there will always
    %  be one).  Continue until the series terminates at a primed zero
    %  that has no starred zero in its column.  Unstar each starred
    %  zero of the series, star each primed zero of the series, erase
    %  all primes and uncover every line in the matrix.  Return to Step 3.
    %**************************************************************************
    rowZ1 = find(starZ==uZc);
    starZ(uZr)=uZc;
    while rowZ1>0
        starZ(rowZ1)=0;
        uZc = primeZ(rowZ1);
        uZr = rowZ1;
        rowZ1 = find(starZ==uZc);
        starZ(uZr)=uZc;
    end
end
% Cost of assignment
rowIdx = find(validRow);
colIdx = find(validCol);
starZ = starZ(1:nRows);
vIdx = starZ <= nCols;
assignment(rowIdx(vIdx)) = colIdx(starZ(vIdx));
pass = assignment(assignment>0);
pass(~diag(validMat(assignment>0,pass))) = 0;
assignment(assignment>0) = pass;
cost = trace(costMat(assignment>0,assignment(assignment>0)));
end

function [minval,rIdx,cIdx]=outerplus(M,x,y)
ny=size(M,2);
minval=inf;
for c=1:ny
    M(:,c)=M(:,c)-(x+y(c));
    minval = min(minval,min(M(:,c)));
end
[rIdx,cIdx]=find(M==minval);
end



%% Enhanced cost matrix function



function cost = enhanced_cost_matrix(stats1, stats2, max_speed, alpha)
    % Position cost (normalized Euclidean distance)
    pos1 = vertcat(stats1.Centroid);
    pos2 = vertcat(stats2.Centroid);
    dist_cost = pdist2(pos1, pos2)/max_speed;
    
    % Intensity cost (if available)
    if isfield(stats1, 'Intensity')
        int1 = [stats1.Intensity]';
        int2 = [stats2.Intensity]';
        int_cost = abs(int1 - int2')/max([int1; int2]);
    else
        int_cost = zeros(size(dist_cost));
    end
    
    % Area cost (if available)
    if isfield(stats1, 'Area')
        area1 = [stats1.Area]';
        area2 = [stats2.Area]';
        area_cost = abs(area1-area2')./max([area1; area2]);
    else
        area_cost = zeros(size(dist_cost));
    end
    
    % Combined cost
    cost = alpha*dist_cost + (1-alpha)*(0.3*int_cost + 0.7*area_cost);
end
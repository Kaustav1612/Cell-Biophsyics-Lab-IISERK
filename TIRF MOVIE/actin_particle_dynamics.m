% This code helps you analyze a movie of Tf. Any one of the can be used as
% a referrence see  variable 'frametostart' on line 34
% objects are detected. centring on centroids, 13 pixel linescan in
% performed and kymograph built for each such puncta. (For Changing linescan width change variable
% 'line_scan_length' on line 138)
% Min_Area_Thresh is the threshold of pixels in objects for considering them as
% puncta.

% This code is modularised with all the function are written as separate
% MATLAB function scripts please keep all the MATLAB Scripts in a single directory

clearvars -except displacement_table;
close all;
clc;
olddir = pwd;
%%
num_channels = 2;
channel_1_dir='I:\Uday_data\02_03_24(CytoD_live_img2)\1um_cytoD\Cell 4\TNFR1_1';
channel_2_dir= 'I:\Uday_data\02_03_24(CytoD_live_img2)\1um_cytoD\Cell 4\Actin_1';
directory = 'K:\new_codes\mobility_analysis\Uday\CytoD(1micM)\004';


channel_2_drift_X = displacement_table.XC_Ch1-displacement_table.XC_Ch2;
channel_2_drift_Y = displacement_table.YC_Ch1-displacement_table.YC_Ch2;

mean_channel_drift_X = mean(channel_2_drift_X);
mean_channel_drift_Y = mean(channel_2_drift_Y);

cd(channel_2_dir);

clear A;

A=dir('*.tif');
nframes=size(A,1); 
for i=1:nframes
    filename=A(i).name;
    image_stack_channel_2_uncorrected(:,:,i)=imread(filename);
end
cd(olddir);
num_frames = size(image_stack_channel_2_uncorrected,3);
image = image_stack_channel_2_uncorrected(:,:,1);
frame_size = [size(image,1),size(image,2)];

for f = 1:num_frames
    % Read frame
    current_frame = image_stack_channel_2_uncorrected(:,:,f);

    % Apply translation correction
    corrected_frame = imtranslate(current_frame, [-mean_channel_drift_X, -mean_channel_drift_Y], 'OutputView', 'same');

    corrected_channel_2(:,:,f) = corrected_frame;
end
close all;

% Results of the pixel shift correction 



clear A;
for k = 1 :num_channels
    if k ==1
        cd(channel_1_dir);
        A=dir('*.tif');
        nframes=size(A,1);
        
        for i=1:nframes
            filename=A(i).name;
            image_stack (:,:,i)=imread(filename);
            image_stack_channel_1 = image_stack (:,:,i);
        end
        hold on;
        redrawing = true;
        while redrawing
            imshow(image_stack (:,:,1),[0 170]);  % Display the black and white image
            title('Draw the cell bounary for TNFR1');
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



        for i =1 : nframes
            current_frame = image_stack (:,:,i);


            % Set the threshold calibration which is the number of frames after
            % which sensitivity for adaptthresh is recalibrated based on the number of focused clustera

            threshold_calibration = 100;

            if rem(i,threshold_calibration)==1
                Min_Area_Thresh = 10;
                Max_Area_Thresh = 500;
                threshIntensity = 10;
                correlation_threshold = 0.8;
                final_image =  zeros(size(current_frame));
                threshold = [];
                for k = 6 : 10
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
            stats = regionprops(connected_components, 'Area', 'Eccentricity', 'PixelIdxList', 'Centroid','EquivDiameter');

            % Compute intensity and store additional properties
            for u = 1:numel(stats)
                pixel_indices = stats(u).PixelIdxList;
                intensity = mean(current_frame(pixel_indices));
                stats(u).Intensity = intensity;
                centroid = stats(u).Centroid;
                stats(u).XC = centroid(1);
                stats(u).YC = centroid(2);
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
            final_stats = regionprops(labeled_norm_img, 'Area', 'Eccentricity', 'PixelIdxList', 'Centroid','EquivDiameter');
            for u = 1:numel(final_stats)
                pixel_indices = final_stats(u).PixelIdxList;
                intensity = mean(current_frame(pixel_indices));
                final_stats(u).Intensity = intensity;
                centroid = final_stats(u).Centroid;
                final_stats(u).XC = centroid(1);
                final_stats(u).YC = centroid(2);
                final_stats(u).EquivDiameter=final_stats(u).EquivDiameter;
            end

            % Save stats for current frame in an all_stats array
            all_stats{i} = final_stats;
            labeled_movie_norm_channel_1(:,:,i)=labeled_image_norm;
            detected_movie_channel_1(:,:,i)= labeled_image_norm;
            labeled_movie(:,:,i)=labeled_img;
            
        end

        all_stats=all_stats';
        all_stats_channel_1 = all_stats;
    end

    %%
    
    clear A;
    close all;
    if k==2
        image_stack_channel_2 = corrected_channel_2;
            
        for i =1 : nframes
            current_frame = image_stack_channel_2(:,:,i);

            % Filtering image using gaussian blur
            gauss_blur=imgaussfilt(current_frame,8,"FilterSize",[13 13]);
            filtered_image = current_frame-gauss_blur;

            % Normalisation of the filtered image
            min_value = min(filtered_image(:));  % Get minimum value of all elements
            max_value = max(filtered_image(:));  % Get maximum value of all elements
            normalized_image = double(filtered_image - min_value) / double(max_value - min_value);

            %Binarising using threshold
            threshold = 0.05;
            binary_image = imbinarize(normalized_image, threshold);

            % Erosion and Dilation
            eroded_image = imerode(binary_image, strel('disk',1));
            dilated_image=imdilate(eroded_image,strel('disk',1));
            BW = dilated_image;
            BW_1=dilated_image;

            image = BW .* mask;
            masked_normalized_image = normalized_image.*mask;

            figure();
            imshow(image,[0 150]);
            title('Masked Image');
            close;

            [componentIndices,stats_label,labeled_image] = detect(image,15,masked_normalized_image);
            skel_image = bwmorph(image, 'skel', Inf);


            labeled_movie_norm_channel_2(:,:,i)=masked_normalized_image;
            detected_movie_channel_2(:,:,i) = image;
            skel_detected_movie_channel_2(:,:,i)= skel_image;

            figure();
            imshow(image,[0 150]);
            close;
        end
    end


    figure();
    for i = 1:nframes
        imshow(labeled_movie_norm_channel_1(:,:,i),[]);
        pause(0.03); % Adjust speed
    end

end
%%
close all;
figure;
imshowpair(image_stack_channel_1(:,:,1),image_stack_channel_2_uncorrected(:,:,1));
title('Uncorrected Channel');
figure;
imshowpair(image_stack_channel_1(:,:,1),image_stack_channel_2(:,:,1));
title('Corrected Channel');



%% Choose an ROI
image_stack_processed = detected_movie_channel_1;
run = 1;
nslot=nframes;
current_frame = detected_movie_channel_1(:,:,1);

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
centroids = cell(nframes, 1);
intensity_details = cell(nframes, 1);
area_details = cell(nframes, 1);
for frame = 1:nframes
    if ~isempty(all_stats_channel_1{frame})
        centroids{frame} = vertcat(all_stats_channel_1{frame}.Centroid);
        intensity_details{frame} = vertcat(all_stats_channel_1{frame}.Intensity);
        area_details{frame} = vertcat(all_stats_channel_1{frame}.Area);
    else
        centroids{frame} = [];
        intensity_details{frame} = [];
        area_details{frame} = [];
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
    cost_matrix = enhanced_cost_matrix(all_stats_channel_1{t}, all_stats_channel_1{t+1}, max_speed_per_frame, alpha);
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
if ~isempty(all_stats_channel_1{1})
    num_particles = length(all_stats_channel_1{1});
    particle_to_track_map{1} = (next_track_id:next_track_id+num_particles-1)';
    next_track_id = next_track_id + num_particles;
end

% The track is built using the matched pairs variable, which stores the
% assigned particle for the next frame

for t = 1:nframes-1
    if ~isempty(matched_pairs{t})
        num_particles_next = length(all_stats_channel_1{t+1});
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
trajectories = struct(...
    'trajectory_id', {}, ...
    'particle_id', {}, ...
    'frames', {}, ...
    'positions', {}, ...
    'intensity', {}, ...
    'area', {}, ...
    'length', {}, ...
    'max_gap', {}, ...
    'valid_frames', {});

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

    % Extract actual measurements of the positions, area and intensity
    measured_positions = zeros(length(valid_frames), 2);
    measured_intensity = zeros(1, length(valid_frames));
    measured_area = zeros(1, length(valid_frames));

    for i = 1:length(valid_frames)
        t = valid_frames(i);
        idx = particle_indices(i);
        measured_positions(i,:) = centroids{t}(idx,:);
        measured_intensity(i) = intensity_details{t}(idx);
        measured_area(i) = area_details{t}(idx);
    end

    % Create full timeline from first to last appearance
    full_frames = min(valid_frames):max(valid_frames);
    full_positions = NaN(length(full_frames), 2);
    full_intensity = NaN(1, length(full_frames));
    full_area = NaN(1, length(full_frames));

    % Fill in measured values
    [~, loc] = ismember(valid_frames, full_frames);
    full_positions(loc,:) = measured_positions;
    full_intensity(loc) = measured_intensity;
    full_area(loc) = measured_area;



    %     % Interpolate missing values
    for dim = 1:2
        full_positions(:,dim) = interp1(full_frames(loc), measured_positions(:,dim), full_frames, 'linear');
    end
    full_intensity = interp1(valid_frames, measured_intensity, full_frames, 'nearest');
    full_area = interp1(valid_frames, measured_area, full_frames, 'nearest');

    % Calculate gap statistics
    gap_lengths = diff(valid_frames) - 1;
    max_gap = max([0, gap_lengths(gap_lengths > 0)]);

    % Store trajectory
    trajectories(end+1) = struct(...
        'trajectory_id', length(trajectories)+1, ...
        'particle_id', track_id, ...
        'frames', full_frames, ...
        'positions', full_positions, ...
        'intensity', full_intensity, ...
        'area', full_area, ...
        'length', length(full_frames), ...
        'max_gap', max_gap, ...
        'valid_frames', valid_frames);
end



%% Filter trajectories by length and gap size
valid_trajectories = trajectories([trajectories.length] >= min_traj_length & ...
    [trajectories.max_gap] <= max_gap);
%% Parameters for merging trajectories
max_merge_gap = 4;          % Maximum allowed frame gap between trajectories
max_merge_distance = 10;     % Maximum spatial distance for merging
min_merged_length = 100;     % Minimum length after merging to keep

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
    new_trajectories = repmat(struct(...
        'trajectory_id', [], ...
        'particle_id', [], ...
        'frames', [], ...
        'positions', [], ...
        'intensity', [], ...
        'area', [], ...
        'length', [], ...
        'max_gap', [], ...
        'valid_frames', []), 0, 1);
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
    new_trajectories = repmat(struct(...
        'trajectory_id', [], ...
        'particle_id', [], ...
        'frames', [], ...
        'positions', [], ...
        'intensity', [], ...
        'area', [], ...
        'length', [], ...
        'max_gap', [], ...
        'valid_frames', []), 0, 1);
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
imshowpair(detected_movie_channel_1(:,:,1), skel_detected_movie_channel_2(:,:,1));
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
x_end = size(detected_movie_channel_1, 2);
y_end = size(detected_movie_channel_1, 1);

sizest=size(all_stats_channel_1{1,1});
np=sizest(1);
line_scan_length = 6;



num_particles_detected=zeros(1,nframes);
max_val=0;
max_index=1;
for i =1 :nframes
    num_particles_detected(i)=size(all_stats_channel_1{i}, 1);
    S=all_stats_channel_1{i};
    if length(all_stats_channel_1{i})>max_val
        max_val = length(all_stats_channel_1{i});
        max_index=i;
    end
end


% Extract XC, YC, and Intensity
x_coords = int16([all_stats{max_index}.XC]);
y_coords = int16([all_stats{max_index}.YC]);
Intparticle = [all_stats{max_index}.Intensity];


max_particles = 0;
for i = 1:nframes
    [max_particles] = max(max_particles, size(all_stats{i}, 1));
end
num_particles_frames = NaN(max_particles,1);
for i = 1:nframes
    num_particles_frames(i)= size(all_stats{i}, 1);
end
[~,index_max]= max(num_particles_frames);
%% Following a single detected object
if max_particles>size(trajectories,2)
    max_particles=size(trajectories,2);
end


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
[new_stats] = update_all_stats(all_stats,index_max,nframes,max_particles,max_movement_allowed_per_frame,total_max_movement);
minimum_particle_size=6;

wrong_indices = find(num_particles_detected>max_particles);
if length(wrong_indices)>0
    num_particles_detected(wrong_indices)=max_particles;
end
[kymograph_position,kymograph_intensity,particle_kymograph,dist_from_centroid] = ...
    generate_kymograph(labeled_movie_norm_channel_1,binary_stack,x_coords,y_coords,new_stats,nframes,minimum_particle_size,...
    max_particles,num_particles_detected,current_frame,X_Centroids,Y_Centroids);

[Centroids_3D,frames,~,~,particle_presence] = follow_track(new_stats, nframes,max_particles);


clear y kymograph_stitched kymograph_stitched_intensity;

[kymograph_stitched,kymograph_stitched_intensity,total_lifetime_normalized] = ...
    get_lifetime(num_particles_detected,kymograph_position,kymograph_intensity,nframes,max_particles);

%% Define thresholds for stuck and mobile particles

mobile_threshold_low = 8;
mobile_threshold_high = 50;
stuck_threshold_low = 0.8 * nframes;
stuck_threshold_high = nframes;

[num_stuck_particles,num_mobile_particles,stuck_image_full,mobile_image_full,stuck_image,mobile_image] = find_stuck(stuck_threshold_low,stuck_threshold_high,mobile_threshold_low,mobile_threshold_high,labeled_movie,num_particles_detected,total_lifetime_normalized,max_index);

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

% Plot size of detected particles
figure(4);
plot(total_lifetime_normalized);
title('Lifetime Variation');
hold on;
close;


% Display stuck and mobile particles
figure(7);
imshowpair(stuck_image, mobile_image);

x_coords_stuck = X_Centroids(total_lifetime_normalized >= stuck_threshold_low & total_lifetime_normalized <= stuck_threshold_high);
y_coords_stuck = Y_Centroids(total_lifetime_normalized >= stuck_threshold_low & total_lifetime_normalized <= stuck_threshold_high);

% Display intensity kymograph
figure(8);
imshow(kymograph_stitched_intensity, [10 200]);
colormap jet;
colorbar;
title("Intensity Kymograph");
filename=sprintf('Intensity Kymograph.fig');
savefig(fullfile(directory,filename));
close;

figure(11);
hold on;
colors = lines(max_particles); % Generate distinguishable colors for each particle
X = Centroids_3D(:, :, 1);  % X centroids
Y = Centroids_3D(:, :, 2);  % Y centroids
Z = Centroids_3D(:, :, 3);  % Frame numbers
for j = 1:max_particles
    % Plot each particle's trajectory in 3D
    plot3(X(j, :), Y(j, :), Z(j, :), 'LineWidth', 2,'Color', colors(j, :), 'DisplayName', sprintf('Particle %d', j));
end
title('Trajectory of particles across frames');
xlabel('X Centroid');
ylabel('Y Centroid');
zlabel('Frame Number (Z)');
legend('show'); % Show legend for the particles
grid on;
filename=sprintf('Trajectoryacrossframes_D1C1.fig');
savefig(fullfile(directory,filename));
hold off;
close;
%%
close all;
pixel = 0.065;
X_Centroids = X_Centroids*pixel;
Y_Centroids = Y_Centroids*pixel;
del_t=0.2;
TAMSD =NaN(size(nframes,max_particles));
for j = 1 : max_particles
    N = size(find(particle_presence(:,j)==1),1);
    for n  = 1 : N
        sum_X=0;
        sum_Y=0;
        for t = 1:n
            sum_X = sum_X+(X_Centroids(j,t) - X_Centroids(j, min(find(particle_presence(:,j)==1))))^2;
            sum_Y = sum_Y+(Y_Centroids(j,t) - Y_Centroids(j, min(find(particle_presence(:,j)==1))))^2;
        end
        if sum_X+sum_Y ~= 0
            TAMSD(n,j)=(sum_X+sum_Y)/((N-n+1)*del_t);

        elseif  sum_X+sum_Y == 0
            TAMSD(n,j)=NaN;
        end
    end
end

% Initialize slope vector
slope_vector = [];

% Time vector (must match TAMSD's time dimension)
a = (1:size(TAMSD, 1))' * del_t;  % Column vector of time points

% Loop through each particle
for j = 1:max_particles
    % Get valid indices where TAMSD > 0 (avoid log(0) or negative values)
    valid_idx = find(TAMSD(:, j) > 0);

    % Skip if insufficient data points
    if length(valid_idx) <= 50
        continue;
    end
    a=a*10^3;
    % Extract valid time and TAMSD values
    valid_time = a(valid_idx);
    valid_TAMSD = TAMSD(valid_idx, j);

    % Linear fit in log10 space: log10(TAMSD) = slope * t + intercept
    p = polyfit(log(valid_time), log10(valid_TAMSD), 1);
    slope = p(1);
    intercept = p(2);

    % Skip if slope is invalid (NaN or zero)
    if isnan(slope) || slope == 0
        continue;
    end

    % Store slope if all conditions are met
    slope_vector = [slope_vector; slope];

    % Generate fitted curve for plotting
    x_fit = 10.^linspace(min(valid_time), max(valid_time), 100);
    y_fit = 10.^(slope * x_fit + intercept);  % Convert back from log10

    % Create log-log plot
    figure(j);


    % Plot raw data (log-log scale)
    loglog(valid_time, valid_TAMSD, 'b-','LineWidth', 1.5,  'DisplayName', 'TAMSD Data');
    hold on;
    % Plot fitted line (log-log scale)
    loglog(x_fit, y_fit, 'ro','MarkerSize', 6,'DisplayName', sprintf('Fit: slope=%.6f', slope));

    % Add labels and title
    title(sprintf('TAMSD for Particle %d (Slope = %.3f)', j, slope));
    xlabel('Time, t (ms)');
    ylabel('TAMSD, \langle \Delta r^2(\tau) \rangle');
    legend('Location', 'best');
    grid on;

    % Save figure
    filename = sprintf('TAMSD_Particle_%d.fig', j);
    savefig(fullfile(directory, filename));
    close;
end

MSD = NaN(max_particles,1);
D = NaN(max_particles,1);
particle_sizes = [];
for j =1 : max_particles
    if length(find(TAMSD(:,j)>0)) > 50
        sum_X=0;
        sum_Y=0;
        N=0;
        X_Particle = X_Centroids(j,:);
        Y_Particle = Y_Centroids(j,:);
        for k =1 :nframes
            if ~isnan(X_Centroids(j,k))
                N = N+1;
                sum_X = sum_X+(X_Centroids(j, k) - mean((~isnan(X_Centroids(j,:)))))^2;
                sum_Y = sum_Y + (Y_Centroids(j, k) -mean((~isnan(Y_Centroids(j,:)))))^2;
            end
        end
        MSD(j,1) = (sum_X+sum_Y)/N;
        D(j,1) = (MSD (j,1))/(4*N*del_t); % micro meter^2/s
    end
end

%%

% [G',G''] = GSER();
sizes = zeros(size(trajectories,2),1);
pixel = 0.065;
for particle = 1 : size(trajectories,2)
    rad=[];
    for k = 2 : size(new_stats,2)
        if particle < size(new_stats{1,k},1)
            if new_stats{1,k}(particle).Area > 0
                rad = [rad,sqrt(new_stats{1,k}(particle).Area/3.14)];
            else
                continue
            end
        end

    end
    sizes(particle) = mean(rad(:))*pixel*1e-6;
end

if size(TAMSD,2) < num_trajectories
    num_trajectories = size(TAMSD,2);
end
T=310;
[Gp, Gpp, Gabs,Gp_mean,Gpp_mean,omega_unique] = mason_weitz_MSD(TAMSD, sizes,T,num_trajectories);
%%
dt=0.2;
results = meanback_infer(trajectories, dt, T);

% Plot all trajectories' MBR curves
figure; hold on
for i = 1:results.ntraj
    if ~isempty(results.MBR{i})
        plot(results.MBR_t{i}, results.MBR{i}, 'Color',[0.7 0.7 0.7]);
    end
end

% (Assuming all MBR_t have the same length; if not, interpolate them first)
tcommon = results.MBR_t{1};  % common lag times

% --- Stack all trajectories ---
ntraj = results.ntraj;
MBRmat = nan(length(tcommon), ntraj);

for i = 1:ntraj
    if isfield(results,'MBR') && ~isempty(results.MBR{i})
        if length(results.MBR{i})==length(tcommon)
            MBRmat(:,i)=results.MBR{i};
        else
            % interpolate onto tcommon if needed:
            MBRmat(:,i)=interp1(results.MBR_t{i},results.MBR{i},tcommon,'linear','extrap');
        end
    end
end

% --- Mean MBR ---
meanMBR = nanmean(MBRmat,2);   % average over all trajectories

% --- Plot mean MBR ---
figure;
plot(tcommon, meanMBR,'k','LineWidth',2)
xlabel('Lag time (s)')
ylabel('Mean Back Relaxation (MBR)')
title('Mean MBR across all trajectories')
grid on
savefig(fullfile(directory, 'Mean_MBR_All_Particles.fig'));
print(fullfile(directory, 'Mean_MBR_All_Particles.png'), '-dpng', '-r300');
close;

K = results.k_per_dim./results.k_radial;
D_trap = results.D_trap;
D = results.D_short;

kB = 1.380649e-23;
for i = 1: numel(tcommon)
    model_MBR(i) = 0.5*(1-(mean(D_trap./D))).*(1-exp(1-(-mean(K.*D)*i)/kB*T));
end

kB = 1.380649e-23;
k_eff = nanmean(results.k_per_dim);          % average stiffness
D_eff = nanmean(results.D_short);            % short-time diffusion
gamma_eff = kB*T / D_eff;                    % friction
tau_rel = gamma_eff / k_eff;                 % relaxation time

% choose amplitude from your first MBR point or just set A=1
A = abs(mean(results.MBR{1}(1)));            % for example

tvec = tcommon;
MBR_model = -A * exp(-tvec / tau_rel);

plot(tvec, meanMBR, 'k','LineWidth',2); hold on
plot(tvec, MBR_model, 'r--','LineWidth',2)
xlabel('Lag time (s)')
ylabel('Mean Back Relaxation')
legend('Data','OU model')


figure;
plot(tcommon, meanMBR,'k','LineWidth',2);
hold on;
plot(tcommon,model_MBR);
xlabel('Lag time (s)')
ylabel('Mean Back Relaxation (MBR)')
title('Mean MBR across all trajectories')
grid on

%%

% Assuming:
% - TAMSD: Matrix of size [nFrames × nParticles]
% - del_t: Time step per frame
% - directory: Path to save plots

time = (1:size(TAMSD, 1))' * del_t;
valid_TAMSD = TAMSD;

mean_TAMSD = mean(valid_TAMSD, 2, 'omitnan');
std_TAMSD = std(valid_TAMSD, 0, 2, 'omitnan');

% Ensure no negative values for log scale
valid_mask = (mean_TAMSD - std_TAMSD) > 0;
time_valid = time;
mean_valid = mean_TAMSD;
std_valid = std_TAMSD;

figure;
hold on;


fill([time_valid; flipud(time_valid)], ...
    [mean_valid - std_valid; flipud(mean_valid + std_valid)], ...
    [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.75);

loglog(time_valid, mean_valid, 'b-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time, \tau (ms)');
ylabel('TAMSD');
title('Mean TAMSD ± 1 STD');
grid on;

% Save plot
savefig(fullfile(directory, 'Mean_TAMSD_All_Particles.fig'));
print(fullfile(directory, 'Mean_TAMSD_All_Particles.png'), '-dpng', '-r300');
close;



%% Time Averaged MSD TAMSD



% Number of particles and frames
[num_frames,num_particles] = size(particle_presence);

% Initialize an array to hold the rate of frames missed
missed_rate = zeros(num_particles, 1);

% Initialize an array to hold the average distance from the centroid
avg_distance = zeros(num_particles, 1);


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
imwrite(detected_movie_channel_1(:,:,1), gifFile, 'LoopCount', Inf, 'DelayTime', 0.1);

% Append subsequent frames
for i = 2:size(detected_movie_channel_1, 3)
    imwrite(detected_movie_channel_1(:,:,i), gifFile, 'WriteMode', 'append', 'DelayTime', 0.1);
end

%% SAVING KEY VARIABLES

% List of stuck particles
%stuck_list = find(total_lifetime_normalized > mobile_threshold_low & total_lifetime_normalized < mobile_threshold_high);
cd(directory);
% Save results
results_file_name = 'all_information';
save(results_file_name, 'kymograph_stitched_intensity', 'kymograph_stitched', 'total_lifetime_normalized','TAMSD','MSD','D','slope_vector','new_stats','Gp_mean','Gpp_mean','omega_unique','tcommon', 'meanMBR','trajectories');
cd(olddir);
function [xc, yc, R] = circfit(x, y)
% CIRCUMCENTER calculates the center (xc, yc) and radius R of a circle
% that best fits a set of 2D points (x, y) using least squares.

A = [-2*x, -2*y, ones(length(x), 1)];
B = -(x.^2 + y.^2);
X = A \ B;
xc = X(1);
yc = X(2);
R = sqrt(X(3) + xc^2 + yc^2);
end

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


% function cost = enhanced_cost_matrix(stats1, stats2, max_speed, alpha)
%     % Position cost with exponential decay
%     pos1 = vertcat(stats1.Centroid);
%     pos2 = vertcat(stats2.Centroid);
%     n1 = size(pos1,1);
%     n2 = size(pos2,1);
%     dist_matrix = pdist2(pos1, pos2);
%
%     % Dynamic distance threshold (density-adaptive)
%     median_dist = median(pdist(pos1));
%     dist_threshold = min(3*max_speed, 5*median_dist);
%
%     % Exponential distance penalty
%     dist_cost = 1 - exp(-(dist_matrix.^2)/(0.25*max_speed^2));
%     dist_cost(dist_matrix > dist_threshold) = Inf;
%
%     % Intensity similarity (mandatory)
%     int1 = [stats1.Intensity]';
%     int2 = [stats2.Intensity]';
%     int_cost = abs(int1 - int2')./max([int1; int2]);
%
%     % Shape features (safe implementation)
%     shape_cost = zeros(n1,n2);
%     if isfield(stats1, 'Solidity') && isfield(stats1, 'Eccentricity') && ...
%        isfield(stats2, 'Solidity') && isfield(stats2, 'Eccentricity')
%
%         % Ensure we have column vectors
%         sol1 = [stats1.Solidity]';
%         ecc1 = [stats1.Eccentricity]';
%         sol2 = [stats2.Solidity]';
%         ecc2 = [stats2.Eccentricity]';
%
%         % Safe shape computation
%         if ~isempty(sol1) && ~isempty(ecc1) && ~isempty(sol2) && ~isempty(ecc2)
%             shape1 = sol1 .* ecc1;
%             shape2 = sol2 .* ecc2;
%             shape_cost = abs(shape1 - shape2');
%         end
%     end
%
%     % Dynamic feature weighting
%     proximity_weight = 1./(1 + exp(-10*(dist_matrix/max_speed - 1.5)));
%     feature_weight = 0.8*(1 - proximity_weight);
%
%     % Density-based isolation penalty
%     isolation_penalty = zeros(n1,n2);
%     if n1 > 1 && n2 > 1
%         density1 = 1./(1 + pdist2(pos1, pos1));
%         density2 = 1./(1 + pdist2(pos2, pos2));
%         isolation_penalty = max(0, 2 - mean(density1,2) - mean(density2,2)');
%     end
%
%     % Final cost with all safeguards
%     cost = alpha*proximity_weight.*dist_cost + ...
%           (1-alpha)* feature_weight.*(0.6*int_cost + 0.4*shape_cost);
%     cost = cost .* (1 + isolation_penalty);
%
%     % Ensure proper matrix size
%     cost = cost(1:n1, 1:n2);
% end

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
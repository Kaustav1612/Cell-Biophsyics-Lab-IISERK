clc; close all;
clearvars -except displacement_table displacement_table trajectories D MSD;
olddir = pwd;
% --- 1. Setup Parameters ---

num_frames = 300; % Hardcoded as requested
num_p2_tracks = size(trajectories, 2);
K_nearest = 20; % Average distance to the 20 closest pixels
channel_2_dir= 'F:\Uday_data\15_02_24(ML7_live_img1)\ml7_10um\Cell 5\Actin_1';
[all_stats,labeled_img_binary] = actin_positions(channel_2_dir,displacement_table,olddir);
image_size = size(labeled_img_binary); 
% Initialize Result Matrix: [ParticleID, Frame] -> Distance
% We use a matrix for speed, then summarize later
dist_matrix = nan(num_p2_tracks, num_frames);

fprintf('Processing %d frames for %d trajectories...\n', num_frames, num_p2_tracks);

% --- 2. Main Loop Over Frames ---
for f = 1:num_frames
    
    % --- Step A: Select correct Particle 1 stats based on your logic ---
    % Note: The order of IF statements matters. >200 must be checked before >100.
    if f > 200
        current_stats = all_stats{1, 3};
    elseif f > 100
        current_stats = all_stats{1, 2};
    else
        current_stats = all_stats{1, 1};
    end
    
    % If no stats for this frame or empty structure, skip
    if isempty(current_stats)
        continue;
    end
    
    % --- Step B: Build Particle 1 Point Cloud ---
    % Collect ALL pixel indices from ALL objects in this frame
    if isstruct(current_stats)
        all_p1_indices = vertcat(current_stats.PixelIdxList);
    else
        continue;
    end
    
    if isempty(all_p1_indices)
        continue;
    end
    
    % Convert linear indices to [X, Y]
    % Note: ind2sub returns [row, col], which maps to [y, x]
    [rows, cols] = ind2sub(image_size, all_p1_indices);
    p1_point_cloud = [cols, rows]; % [X, Y]
    
    % --- Step C: Get Particle 2 Positions for this Frame ---
    % We find all Particle 2s that exist in frame 'f'
    
    p2_batch_xy = [];
    p2_batch_ids = [];
    
    for i = 1:num_p2_tracks
        pos_data = trajectories(i).positions;
        
        % Check if trajectory extends to this frame
        if f <= size(pos_data, 1)
            xy = pos_data(f, :);
            % Check if valid (not NaN)
            if ~any(isnan(xy))
                p2_batch_xy = [p2_batch_xy; xy]; %#ok<AGROW>
                p2_batch_ids = [p2_batch_ids; i]; %#ok<AGROW>
            end
        end
    end
    
    if isempty(p2_batch_xy)
        continue;
    end
    
    % --- Step D: Calculate Distances ---
    % Find K nearest pixels in Particle 1 cloud for every Particle 2
    [~, dists] = knnsearch(p1_point_cloud, p2_batch_xy, 'K', K_nearest);
    
    % Average the K distances
    avg_dists = mean(dists, 2);
    
    % Store in matrix
    for j = 1:length(p2_batch_ids)
        p2_id = p2_batch_ids(j);
        dist_matrix(p2_id, f) = avg_dists(j);
    end
    
    % Progress update
    if mod(f, 50) == 0
        fprintf('Frame %d processed.\n', f);
    end
end

% --- 3. Compile Statistics ---
results_table = table();
results_table.Particle2_ID = (1:num_p2_tracks)';
results_table.Avg_Min_Distance = nan(num_p2_tracks, 1);
results_table.Std_Dev = nan(num_p2_tracks, 1);
results_table.Frames_Tracked = zeros(num_p2_tracks, 1);

for i = 1:num_p2_tracks
    % Extract valid measurements (non-NaN)
    d_list = dist_matrix(i, :);
    d_list = d_list(~isnan(d_list)); 
    
    if ~isempty(d_list)
        results_table.Avg_Min_Distance(i) = mean(d_list);
        results_table.Std_Dev(i) = std(d_list);
        results_table.Frames_Tracked(i) = length(d_list);
    else
        results_table.Frames_Tracked(i) = 0;
    end
end

% --- 4. Plotting ---
disp('Calculation Complete.');
% disp(results_table);

% Ensure MSD and Distance vectors are same size before plotting
% Only plot for particles that have both valid MSD and Distance data
valid_indices = ~isnan(results_table.Avg_Min_Distance) & (results_table.Frames_Tracked > 0);
if length(MSD) == length(valid_indices)
    % Only keep indices where we have MSD data (assuming MSD is a vector matching trajectories)
    valid_indices = valid_indices & ~isnan(MSD(:)); 
    
    y_data = MSD(valid_indices);
    x_data = results_table.Avg_Min_Distance(valid_indices);
    
    figure('Color', 'w'); 
    scatter(x_data .* 0.065, y_data .* 0.065, 50, 'filled', 'MarkerFaceAlpha', 0.6);
    title('Correlation: MSD vs Distance to Structure');
    xlabel('Avg Distance to Structure (\mum)');
    ylabel('MSD (\mum^2)');
    grid on;
    
    % Optional: Add trendline
    h = lsline;
    set(h, 'LineWidth', 2, 'Color', 'r');
else
    warning('MSD array size (%d) does not match number of trajectories (%d). Cannot plot scatter.', length(MSD), num_p2_tracks);
end
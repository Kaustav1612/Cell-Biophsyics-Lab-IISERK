clear all;
close all;
clc;

%% Parameters (adjusted for better visibility)
width = 400;       % Image width (px)
height = 400;      % Image height (px)
nframes = 600;    
n_objects = randi([20,30]);
noise_level = 0.2; 
object_radius = 2; 
movement_speed = 0.2; 


%% Initialize objects
%% Initialize objects with Brownian motion parameters
positions = [randi([object_radius+5, width-object_radius-5], n_objects, 1), ...
             randi([object_radius+5, height-object_radius-5], n_objects, 1)];

% Brownian motion parameters
diffusion_coefficient = 0.4;  % Controls how fast particles diffuse (pixels^2/frame)
dt = 1;                      % Time step (1 frame)
step_std = sqrt(2*diffusion_coefficient*dt); % Standard deviation of step size
angles = 2*pi*rand(n_objects,1);
intensities = 0.7 + 0.3*rand(n_objects,1); % Brighter objects (70-100%)
%% Generate frames with Brownian motion
for f = 1:nframes
    frame = zeros(height, width);
    
    % Brownian motion update
    if f > 1
        % Random displacements (normally distributed)
        displacements = step_std * randn(n_objects, 2);
        new_positions = positions + displacements;
        
        % Boundary handling (reflective)
        % X boundaries
        out_left = new_positions(:,1) < object_radius+5;
        out_right = new_positions(:,1) > width-object_radius-5;
        new_positions(out_left,1) = 2*(object_radius+5) - new_positions(out_left,1);
        new_positions(out_right,1) = 2*(width-object_radius-5) - new_positions(out_right,1);
        
        % Y boundaries
        out_top = new_positions(:,2) < object_radius+5;
        out_bottom = new_positions(:,2) > height-object_radius-5;
        new_positions(out_top,2) = 2*(object_radius+5) - new_positions(out_top,2);
        new_positions(out_bottom,2) = 2*(height-object_radius-5) - new_positions(out_bottom,2);
        
        positions = new_positions;
    end
    
    centroid_history(:,:,f) = positions; % Store positions
    
    % Draw objects (same as before)
    [X,Y] = meshgrid(1:width, 1:height);
    for o = 1:n_objects
        x0 = positions(o,1);
        y0 = positions(o,2);
        gaussian = intensities(o) * exp(-((X-x0).^2 + (Y-y0).^2)/(2*object_radius^2));
        frame = frame + gaussian;
    end
    
    % Add noise and store (same as before)
    frame = imnoise(frame, 'poisson');
    frame = frame + noise_level*max(frame(:))*randn(size(frame));
    movie_stack(:,:,f) = uint16(2^12 * frame/max(frame(:))); 
end

image_stack=movie_stack;
hold on;
redrawing = true;
while redrawing
  imshow(image_stack (:,:,1),[]);  % Display the black and white image
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
  
%% GET ALL CENTROIDS OR READ ALL THE CENTROIDS

 
ThreshA = 15;         % Reduced minimum area to detect smaller puncta
ThreshB = 1000;      % Increased max area for variability
threshIntensity = 5; % Lowered intensity threshold to detect fainter puncta
[all_stats,detected_movie,labeled_movie_norm] = ObjectDetector(image_stack,mask,nframes,ThreshA,ThreshB,threshIntensity);

    
%% Visualization (verify all objects)
figure;
% Max projection (all objects should appear)
max_proj = max(movie_stack, [], 3);
imshow(max_proj, []);
title(['Max Projection - ', num2str(n_objects), ' objects']);
hold on;
plot(centroid_history(:,1,1), centroid_history(:,2,1), 'r+'); % Initial positions

% Animate with object counting
figure;
for f = 1:nframes
    imshow(movie_stack(:,:,f), []);
    hold on;
    plot(centroid_history(:,1,f), centroid_history(:,2,f), 'ro', 'MarkerSize', 4);
    title(['Frame ', num2str(f), ' | Objects: ', num2str(n_objects)]);
    hold off;
    pause(0.03); % Adjust speed
end


%% Choose an ROI
image_stack_processed = detected_movie;
run = 1;
nslot=nframes;

% szshl=zeros(nframes/nslot);
current_frame = detected_movie(:,:,1);

x_start = 1;
y_start = 1;
x_end = size(detected_movie, 2);
y_end = size(detected_movie, 1);
    
sizest=size(all_stats{1,1});
np=sizest(1);
line_scan_length = 6;
  
%% GET ALL THE OBJECTS AND GENERATE KYMOGRAPH AND STICH THEM TOGETHER
% Initialize parameters
num_particles_detected = zeros(1, nframes);
max_val = 0;
max_index = 1;
for i = 1:nframes
    num_particles_detected(i) = size(all_stats{i}, 1);
    if length(all_stats{i}) > max_val
        max_val = length(all_stats{i});
        max_index = i;
    end    
end

max_speed_per_frame = 3;    % pixels/frame (now consistent with 5px max_step)
cost_cutoff = 3;              % Maximum allowed matching cost
min_traj_length = 10;         % Minimum frames for valid trajectory
max_gap = 5;                  % Maximum allowed gap (frames)
max_step_distance = 7;        % Maximum distance between successive points
alpha = 1;                  % Position vs. features weight

%% Modified main tracking loop
centroids = cell(nframes, 1);
for frame = 1:nframes
    if ~isempty(all_stats{frame})
        centroids{frame} = vertcat(all_stats{frame}.Centroid);
    else
        centroids{frame} = []; 
    end
end

matched_pairs = cell(nframes-1, 1);
for t = 1:nframes-1
    pts_t = centroids{t};
    pts_t1 = centroids{t+1};
    if isempty(pts_t) || isempty(pts_t1)
        continue;
    end
    
    cost_matrix = enhanced_cost_matrix(all_stats{t}, all_stats{t+1}, max_speed_per_frame, alpha);
    [assignments, ~] = munkres(cost_matrix);
    
    % Enhanced validation
    for k = 1:length(assignments)
        if assignments(k) > 0
            % Calculate actual speed (px/frame)
            distance = norm(pts_t(k,:) - pts_t1(assignments(k),:));
            frame_diff = 1; % Since we're doing frame-to-frame
            actual_speed = distance/frame_diff;
            
            % Strict validation
            if cost_matrix(k, assignments(k)) > cost_cutoff || ...
               actual_speed > max_speed_per_frame
                assignments(k) = 0; % Reject
            end
        end
    end
    matched_pairs{t} = assignments;
end

%% Enhanced trajectory building

tracking_details = cell(max_val, nframes);
for t = 1:nframes-1
    if isempty(matched_pairs{t}), continue; end
    
    for k = 1:length(matched_pairs{t})
        if matched_pairs{t}(k) > 0
            tracking_details{k,t} = centroids{t}(k,:);
            tracking_details{k,t+1} = centroids{t+1}(matched_pairs{t}(k),:);
        end
    end
end

trajectories = {}; % Initialize as cell array
for particle = 1:max_val
    frames = find(~cellfun(@isempty, tracking_details(particle,:)));
    if isempty(frames), continue; end
    
    segments = [];
    current_segment = [];
    positions=[];
    current_segment{end+1} = frames(1);
    current_positions = tracking_details{particle, frames(1)};
    positions(end+1,:)=tracking_details{particle, frames(1)};
    for i = 2:length(frames)
        prev_pos = tracking_details{particle, frames(i-1)};
        curr_pos = tracking_details{particle, frames(i)};
        distance = norm(curr_pos - prev_pos);
        frame_diff = frames(i) - frames(i-1);
        
        if (frame_diff <= max_gap) && ...
           (distance <= max_step_distance)
            % Valid continuation
            current_segment{end+1} = frames(i);
            positions(end+1,:) = curr_pos;
        else
            % Terminate current segment if valid
            numerical_segment = [current_segment{:}];
            segments{end+1} = struct(...
                    'frames', numerical_segment, ...
                    'positions', positions, ...
                    'length', length(current_segment), ...
                    'max_step', max(sqrt(sum(diff(positions(:)).^2,2))), ...
                    'max_gap', max(diff(numerical_segment)-1));
            
            % Start new segment
            current_segment=[];
            positions=[];
            current_segment{end+1} = frames(i);
            positions(end+1,:)= curr_pos;
        end
    end
    
        segments{end+1} = struct(...
            'frames', [current_segment{:}], ...
            'positions', positions, ...
            'length', length(current_segment),...
            'max_gap', max(diff([current_segment{:}])-1));
 

    for seg = 1:length(segments)
        if ~isempty(segments{seg}) && isfield(segments{seg}, 'positions')
            trajectories{end+1} = struct(...
                'segment_number',seg,...
                'frames',[segments{seg}.frames],...
                'particle_id', particle, ...
                'positions', [segments{seg}.positions], ...
                'length', [segments{seg}.length], ...
                'max_gap', [segments{seg}.max_gap]);
        end
    end
end

for i = 1:length(trajectories)
    trajectories{i}.trajectory_id = i;
end

trajectories = trajectories(cellfun(@(x) x.length >= min_traj_length, trajectories));
%% Parameters for merging trajectories


    % Parameters
    max_merge_distance = 6;     % Maximum spatial distance for merging
    min_merged_length = 150;    % Minimum length after merging to keep

    % Function to calculate distance between trajectory endpoints
    endpoint_distance = @(t1,t2) min([
        pdist2(t1.positions(end,:),t2.positions(1,:)),   % end of t1 to start of t2
        pdist2(t1.positions(1,:),t2.positions(end,:))     % start of t1 to end of t2
    ]);

    % Initialize merge tracking
    merge_iteration = 0;
    max_iterations = 10;         % Safety limit to prevent infinite loops

    % Main merging loop
    while merge_iteration < max_iterations
        merge_iteration = merge_iteration + 1;
        fprintf('\n=== Merge Iteration %d ===\n', merge_iteration);
        
        % Initialize variables for this iteration
        n_trajectories = length(trajectories);
        merged_indices = false(1, n_trajectories);
        new_trajectories = {};
        merge_count = 0;
        
        % Sort trajectories by length (longest first) to prioritize merging longer trajectories
        [~, sort_idx] = sort(cellfun(@(x) x.length, trajectories), 'descend');
        sorted_trajectories = trajectories(sort_idx);
        
        % Try to merge all possible trajectory pairs
        for i = 1:n_trajectories
            current_idx = sort_idx(i);
            if merged_indices(current_idx), continue; end
            
            best_merge = [];
            best_j = 0;
            best_dist = 30;
            
            % Find the best candidate to merge with current trajectory
            for j = 1:n_trajectories
                candidate_idx = sort_idx(j);
                if current_idx == candidate_idx || merged_indices(1,candidate_idx), continue; end
                
                % Calculate distance and gap between trajectories
                dist = endpoint_distance(trajectories{current_idx}, trajectories{candidate_idx});
                  
                
                % Check merge conditions
                if trajectories{current_idx}.particle_id==trajectories{candidate_idx}.particle_id...
                        && dist <= max_merge_distance
                    if dist < best_dist 
                        best_dist = dist;
                        best_j = candidate_idx;
                    end
                end
            end
            
            % Perform the merge if a candidate was found
            if best_j > 0
                
                % Determine merge direction
                if trajectories{best_j}.segment_number > trajectories{current_idx}.segment_number
                    % t1 then t2 (normal case)
                    new_pos = [trajectories{current_idx}.positions; trajectories{best_j}.positions];
                    new_frames = [trajectories{current_idx}.frames, trajectories{best_j}.frames];
                elseif trajectories{current_idx}.segment_number > trajectories{best_j}.segment_number
                    % t2 then t1 (reverse case)
                    new_pos = [trajectories{best_j}.positions; trajectories{current_idx}.positions];
                    new_frames = [trajectories{best_j}.frames, trajectories{current_idx}.frames];
                end
 
                % Verify merged trajectory consistency
                if  length(new_frames) == size(new_pos,1)...
                  && size(new_pos,1)< nframes

                    seg = min(trajectories{best_j}.segment_number,trajectories{current_idx}.segment_number);
                    % Create merged trajectory
                    merged_traj = struct(...
                        'segment_number',seg,... 
                        'particle_id',trajectories{current_idx}.particle_id, ...
                        'frames', new_frames, ...
                        'positions', new_pos, ...
                        'length', length(new_frames), ...
                        'merged_from', [current_idx, best_j]);
                    
                    fprintf('Merged Trajectories %d and %d (dist=%.2f) with (particle id = %d)\n',...
                        current_idx, best_j, best_dist,trajectories{current_idx}.particle_id);
                    
                    % Add to new trajectories
                    new_trajectories{end+1} = merged_traj;
                    merged_indices(1,candidate_idx) = true;
                    merge_count = merge_count + 1;
                else
                    warning('Inconsistent merge between trajectories %d and %d', current_idx, best_j);
                end
            else
                % Keep unmerged trajectory
                new_trajectories{end+1} = trajectories{current_idx};
                merged_indices(1,current_idx) = true;
            end
        end
        
        % Only break if no merges occurred after first iteration
        if merge_count == 0 && merge_iteration > 1
            break;
        end
        
        % Update trajectories for next iteration
        trajectories = new_trajectories;
    end

    % Final filtering by minimum length
    % trajectories = trajectories(cellfun(@(x) x.length >= min_merged_length, trajectories));

    % Re-number trajectories
    for i = 1:length(trajectories)
        trajectories{i}.trajectory_id = i;
    end

    % Diagnostic output
    fprintf('\n=== Final Trajectory Statistics ===\n');
    fprintf('Total trajectories after merging: %d\n', length(trajectories));
    if ~isempty(trajectories)
        lengths = cellfun(@(x) x.length, trajectories);
        
        fprintf('Length stats: Min=%d, Mean=%.1f, Max=%d\n', ...
            min(lengths), mean(lengths), max(lengths));

    end

%% Enhanced Visualization

trajectories = trajectories(cellfun(@(x) x.length >= min_merged_length, trajectories));

figure;
imshow(image_stack(:,:,1), []);
hold on;
title('Trajectories Color-Coded by Length');
xlabel('X position (px)');
ylabel('Y position (px)');

% Create legend handles
legend_handles = [];
legend_labels = {};

% Define color scheme
color_scheme = struct(...
    'long', [0 0.447 0.741], ...    % Blue
    'medium', [0.466 0.674 0.188], ... % Green
    'short', [0.85 0.325 0.098]);      % Red

for i = 1:length(trajectories)
    traj = trajectories{i};
    
    % Determine trajectory category
    if traj.length >= 100
        category = 'long';
        lw = 1.5;
    elseif traj.length >= 50
        category = 'medium';
        lw = 1;
    else
        category = 'short';
        lw = 0.5;
    end
    
    % Plot trajectory
    h = plot(traj.positions(:,1), traj.positions(:,2), ...
        'Color', color_scheme.(category), 'LineWidth', lw);
    
    % Plot start/end markers
    plot(traj.positions(1,1), traj.positions(1,2), 'o', ...
        'MarkerFaceColor', color_scheme.(category), 'MarkerEdgeColor', 'k');
    plot(traj.positions(end,1), traj.positions(end,2), 's', ...
        'MarkerFaceColor', color_scheme.(category), 'MarkerEdgeColor', 'k');
    
end

% Add legend and grid
legend(legend_handles, legend_labels, 'Location', 'best');
grid on;

%% 


[kymograph_position,kymograph_intensity,particle_kymograph,dist_from_centroid] = ...
    generate_kymograph(labeled_movie_norm,binary_stack,x_coords,y_coords,new_stats,nframes,minimum_particle_size,...
    max_particles,num_particles_detected,current_frame,particle_presence,X_Centroids,Y_Centroids,frame_centroid);

%% SET THRESHOLD to get at time at a region and get "lifetime" along the kymograph
clear y kymograph_stitched kymograph_stitched_intensity;

[kymograph_stitched,kymograph_stitched_intensity,total_lifetime_normalized] = ...
    get_lifetime(num_particles_detected,kymograph_position,kymograph_intensity,particle_presence,nframes,max_particles);

%% Define thresholds for stuck and mobile particles

mobile_threshold_low = 8;
mobile_threshold_high = 50;
stuck_threshold_low = 0.8 * nframes;
stuck_threshold_high = nframes;

[num_stuck_particles,num_mobile_particles,stuck_image_full,mobile_image_full,stuck_image,mobile_image] = find_stuck(stuck_threshold_low,stuck_threshold_high,mobile_threshold_low,mobile_threshold_high,labeled_movie,num_particles_detected,total_lifetime_normalized,max_index);
 
%% PLOTS 
directory = "D:\Kaustav\CODES SET\TIRF MOVIE\CAVEOLAE & TNFR1 (TANMOY & JIBITESH)\tnfr1 data\Dish 3\HYPO_004";

figure(2); 
imshow(labeled_movie(:,:,max_index), []); 
colormap jet; 
colorbar; 
title("Labelled endosomes");
filename=sprintf('Labeled Objects.fig');
% savefig(fullfile(directory,filename));
% close;

figure(3);
imshow(kymograph_stitched);
colormap jet;
title("DETECTION KYMOGRAPH X");
image_data = getframe(gcf).cdata;
filename=sprintf('Detection Kymograph.fig');
% savefig(fullfile(directory,filename));
% close;

% Plot size of detected particles
figure(4);
plot(total_lifetime_normalized);
title('Lifetime Variation');
hold on;
close;

% Display mobile particles image
% figure(5);
% imshow(stuck_image_full, []);
% title('Stuck Particles');
% figure(6);
% imshow(mobile_image_full,[]);
% title('Mobile Particles');

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
% savefig(fullfile(directory,filename));
% close;

% Plot the trajectories
% figure(9);
% hold on;
% colors = lines(max_particles); % Generate distinguishable colors for each particle
% for j =1 :max_particles
%    plot(X_Centroids(j, :),frames, 'Color', colors(j, :), 'DisplayName', sprintf('Particle %d', j));
% end
% title('X-Trajectory of the Particle');
% hold off;
% close;
% 
% figure(10);
% hold on;
% colors = lines(max_particles); % Generate distinguishable colors for each particle
% for j =1 :max_particles
%     plot(X_Centroids(j, :), Y_Centroids(j, :),'Color', colors(j, :), 'DisplayName', sprintf('Particle %d', j));
% end
% title('Trajectory of the Particle');
% hold off;
% filename=sprintf('Trajectory_D1C1.fig');
% savefig(fullfile(directory,filename));
% close;



% centroid_intensity = ones(size(max_particles,nframes)).*average_intensity;
%  for j =1 :max_particles
%      for k =1:nframes
%          
%      end
%  end
% class = NaN(max_particles,1);
% for j =1 :max_particles
%     counter = 0;
%     missed_frame=0;
%     index=0;
%     for i  = 1 :nframes
%         if ~isnan(X_Centroids(j,i)) && ~isnan(Y_Centroids(j,i))     
%             counter=counter+1;
%         else 
%             if counter >= 30
%                 index= i;
%                 missed_frame = missed_frame+1;
%             elseif counter>=30 && missed_frame > 4
%                 counter=0;
%                 missed_frame=0;
%             else 
%                 counter=0;
%                  missed_frame=0;
%             end
%         end
% 
%     end
%     if counter >= (0.167*nframes) && counter < (0.83*nframes)
%         if index < (0.5*nframes)
%             class(j)= 1;
%         elseif index > (0.5*nframes) 
%             class(j) = 2;
%         end
%     elseif counter > (0.9167*nframes)
%         class(j)= 3;
%     else 
%         class(j)=4;
%     end
% end

% Class 1  = Exocytosis
% Class 2 = Endocytosis
% Class 3 = Purely Stuck in membrane (Postional Jitter)
% Class 4 = Kiss & Run

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
% savefig(fullfile(directory,filename));
% hold off;
% close;



 TAMSD=NaN(size(nframes,max_particles));
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
             TAMSD(n,j)=(sum_X+sum_Y)*(4225)/(N-n+1);
        elseif  sum_X+sum_Y == 0
             TAMSD(n,j)=NaN;
        end        
     end 
 end

del_t=1000; 
slope_vector =[];
a = log10(frames*del_t);
%  figure(12);
 for j = 1 : max_particles
    p = polyfit(a(1:length(find(TAMSD(:,j)>0))),log10(TAMSD(find(TAMSD(:,j)>0),j)),1);
    slope = p(1);
    intercept = p(2);
    if min(TAMSD(find(TAMSD(:,j)>0),j)) >= 0 
    x_fit = linspace(0,length(find(TAMSD(:,j)>0))*del_t, 100);
    y_fit =  slope*x_fit + intercept;
        if slope ~=0 && length(find(TAMSD(:,j)>0)) > 30
            if ~isnan(slope)
                slope_vector = [slope_vector;slope];    
                figure(j);
                loglog(a(1:length(find(TAMSD(:,j)>0))),TAMSD(find(TAMSD(:,j)>0),j));
                title("");
                xlabel('log(t) (ms)');
                ylabel('log(avg(MSD(t))) ');
                filename=sprintf('Particle_%d.fig',j);
%                 savefig(fullfile(directory,filename));
%                 close;
            end   
           
        end
    end
 end

MSD = NaN(max_particles,1);
D = NaN(max_particles,1);

for j =1 : max_particles
     if length(find(TAMSD(:,j)>0)) > 30
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
     D(j,1) = (MSD (j,1)*4225)/(4*N*del_t);
     end
end


%% Time Averaged MSD TAMSD



% Number of particles and frames
[num_frames,num_particles] = size(particle_presence);

% Initialize an array to hold the rate of frames missed
missed_rate = zeros(num_particles, 1);

% Initialize an array to hold the average distance from the centroid
avg_distance = zeros(num_particles, 1);


for i = 1:nframes
    % Calculate the number of missed frames for the current particle
    num_frames_miss = sum(particle_presence(i, :) == 0);
    
    % Calculate the rate of frames missed
    missed_rate(i) = num_frames_miss;
    
    % Calculate the average distance from the centroid across all frames
    avg_distance(i) = mean(dist_from_centroid(i, :));
end

% Number of frames missed vs distance from the centroid
figure(13);
plot(avg_distance, missed_rate, 'o');
xlabel('Distance from Centroid');
ylabel('No. of Frames Missed');
title('Total Number of Frames Missed vs. Distance from Centroid');
grid on;
filename=sprintf('Frame Missed Plot_D1C1.fig');
% savefig(fullfile(directory,filename));
% close;


gifFile = 'TNFR1_DISH3_hypo004.gif';

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

% List of stuck particles
%stuck_list = find(total_lifetime_normalized > mobile_threshold_low & total_lifetime_normalized < mobile_threshold_high);

% Save results
results_file_name = 'TNFR1_DISH3_hypo004';
save(results_file_name, 'kymograph_stitched_intensity', 'kymograph_stitched', 'total_lifetime_normalized','labeled_image','class','MSD','D','slope_vector','new_stats');

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
    cost = alpha*dist_cost + (1-alpha)*(0.6*int_cost + 0.4*area_cost);
end



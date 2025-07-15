% This code helps you analyze a movie of Tf. Any one of the can be used as
% a referrence see  variable 'frametostart' on line 34
% objects are detected. centring on centroids, 13 pixel linescan in
% performed and kymograph built for each such puncta. (For Chaning linescan width change variable 
% 'line_scan_length' on line 138)
% Min_Area_Thresh is the threshold of pixels in objects for considering them as
% puncta.

% This code is modularised with all the function are written as separate
% MATLAB function scripts please keep all the MATLAB Scripts in a single directory



clear all;
close all;
clc;

centroids = [];
connectivity = 1;
connection_type = 'TNFR1';
cell_number = '1';

% Give your correct data directory
olddir=cd('I:\Jibitesh\TIRF\041223\Dish3\Exported\hypo_004\set');

% PUT 'N' IF YOU DO NOT HAVE THE STATS DATA IN THE DIRECTORY
% PUT 'Y' IF YOU HAVE THE STAS DATA IN THE DIRECTORY
GETALLCEN='N';
clc
close all

Min_Area_Thresh=15;
Max_Area_Thresh=500; 
int_thresh=20;

frametostart=1;
A=dir('*.tif');
nframes=200; 
% nframes=30;
for i=1:nframes
    filename=A(i).name;
    image_stack (:,:,i)=imread(filename);
end
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
  
    [B, ~] = bwboundaries(mask, 'noholes');
    boundary_points = B{1};  % Extract the boundary points
    x_coords_mask = boundary_points(:, 2);  % x-coordinates
    y_coords_mask = boundary_points(:, 1);  % y-coordinates

    % Fit a circle to the boundary points
    [xc, yc, R] = circfit(x_coords_mask, y_coords_mask);
    frame_centroid=[xc,yc];
  cd(olddir);  
%% GET ALL CENTROIDS OR READ ALL THE CENTROIDS
sfname=[connection_type '_ALLSTATS' cell_number 'till framr' num2str(nframes, '%3d')]; % change this for every ROI
best_sensitivity = 0; 
if GETALLCEN=='N' 
for i =1 : nframes 
    current_frame = image_stack (:,:,i);
 
    
    % Set the threshold calibration which is the number of frames after
    % which sensitivity for adaptthresh is recalibrated based on the number of focused clustera
    
    threshold_calibration = 100;
    
      if rem(i,threshold_calibration)==1        
        Min_Area_Thresh = 15;
%         Max_Area_Thresh = 70; 
        correlation_threshold = 0.86;
        final_image =  zeros(size(current_frame));
        threshold = [];
        for k = 6 : 12    
            sigma = k;
            line_scan_length = k-2;
            [final_image,gauss_blur,threshold] = ObjectDetector_un(current_frame,correlation_threshold,line_scan_length,mask,sigma,final_image,threshold,Min_Area_Thresh,Max_Area_Thresh);
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
        imshow(binary_stack(:,:,i));
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
        % Display labeled endosomes 
        figure(1);
        imshow(labeled_image_norm,[],colormap=jet); 
        title("Labelled punctas");
        close;
        
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
        labeled_movie(:,:,i)=labeled_img;
        labeled_movie_norm(:,:,i)=labeled_image_norm;
        detected_movie(:,:,i)= labeled_image_norm;
end
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
x_coords = int16([all_stats{1,max_index}.XC]);
y_coords = int16([all_stats{1,max_index}.YC]);
Intparticle = [all_stats{1,max_index}.Intensity];

 
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
total_max_movement = 10;
max_movement_allowed_per_frame = 1;
[new_stats] = update_all_stats(all_stats,index_max,nframes,max_particles,max_movement_allowed_per_frame,total_max_movement);
[Centroids_3D,frames,X_Centroids,Y_Centroids,particle_presence] = follow_track(new_stats, nframes,max_particles); 
minimum_particle_size=6;
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
savefig(fullfile(directory,filename));
close;

% Plot the trajectories
figure(9);
hold on;
colors = lines(max_particles); % Generate distinguishable colors for each particle
for j =1 :max_particles
   plot(X_Centroids(j, :),frames, 'Color', colors(j, :), 'DisplayName', sprintf('Particle %d', j));
end
title('X-Trajectory of the Particle');
hold off;
close;

figure(10);
hold on;
colors = lines(max_particles); % Generate distinguishable colors for each particle
for j =1 :max_particles
    plot(X_Centroids(j, :), Y_Centroids(j, :),'Color', colors(j, :), 'DisplayName', sprintf('Particle %d', j));
end
title('Trajectory of the Particle');
hold off;
filename=sprintf('Trajectory_D1C1.fig');
savefig(fullfile(directory,filename));
close;



% centroid_intensity = ones(size(max_particles,nframes)).*average_intensity;
%  for j =1 :max_particles
%      for k =1:nframes
%          
%      end
%  end
class = NaN(max_particles,1);
for j =1 :max_particles
    counter = 0;
    missed_frame=0;
    index=0;
    for i  = 1 :nframes
        if ~isnan(X_Centroids(j,i)) && ~isnan(Y_Centroids(j,i))     
            counter=counter+1;
        else 
            if counter >= 30
                index= i;
                missed_frame = missed_frame+1;
            elseif counter>=30 && missed_frame > 4
                counter=0;
                missed_frame=0;
            else 
                counter=0;
                 missed_frame=0;
            end
        end

    end
    if counter >= (0.167*nframes) && counter < (0.83*nframes)
        if index < (0.5*nframes)
            class(j)= 1;
        elseif index > (0.5*nframes) 
            class(j) = 2;
        end
    elseif counter > (0.9167*nframes)
        class(j)= 3;
    else 
        class(j)=4;
    end
end

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
savefig(fullfile(directory,filename));
hold off;
close;



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
                savefig(fullfile(directory,filename));
                close;
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
savefig(fullfile(directory,filename));
close;


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

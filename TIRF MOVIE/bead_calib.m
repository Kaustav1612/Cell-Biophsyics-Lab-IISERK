%% Bead movement calibration

num_channels = 2;
olddir = pwd;

all_stats = cell(1, num_channels);

for i = 1:num_channels
    % Select directory and load image
    current_frame_dir = uigetdir('I:\Uday_data', sprintf('Select Channel %d Folder', i));
    cd(current_frame_dir);
    A = dir('*.tif');
    filename = A(1).name; % assuming only one image per channel folder
    current_frame = imread(filename);

    figure;
    imshow(current_frame, []);
    title(sprintf('Draw ROI for Channel %d', i));
    redrawing = true;
    while redrawing
        mask_handle = drawfreehand('Color', 'r', 'LineWidth', 2);
        choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
        redrawing = strcmp(choice, "Yes");
        if redrawing
            delete(mask_handle);
        end
    end
    mask = mask_handle.createMask;
    close;
    cd(olddir);

    % --- Object detection ---
    Min_Area_Thresh = 15;
    Max_Area_Thresh = 200;
    threshIntensity = 10;
    correlation_threshold = 0.85;
    final_image = zeros(size(current_frame));

    for k = 6:10
        sigma = k;
        line_scan_length = k - 2;
        [final_image, gauss_blur, threshold] = ObjectDetector_un( ...
            current_frame, correlation_threshold, line_scan_length, ...
            mask, sigma, final_image, threshIntensity, ...
            Min_Area_Thresh, Max_Area_Thresh);
    end
    figure();
    imshow(final_image,[]);

    % Optional: dilation if used
    dilated_image = imdilate(final_image > 0, strel('disk', 1));

    % Connected component analysis
    connected_components = bwconncomp(dilated_image);
    stats = regionprops(connected_components, 'Area', 'Eccentricity', ...
                        'PixelIdxList', 'Centroid', 'EquivDiameter');

    % Compute intensity and store centroid info
    for u = 1:numel(stats)
        pixel_indices = stats(u).PixelIdxList;
        intensity = mean(current_frame(pixel_indices));
        stats(u).Intensity = intensity;
        stats(u).XC = stats(u).Centroid(1);
        stats(u).YC = stats(u).Centroid(2);
    end

    all_stats{i} = stats; % store for this channel
end

stats_ch1 = all_stats{1};
stats_ch2 = all_stats{2};

% Extract centroids
centroids1 = [[stats_ch1.XC]', [stats_ch1.YC]'];
centroids2 = [[stats_ch2.XC]', [stats_ch2.YC]'];

% Match objects by nearest centroid (min Euclidean distance)
[idx, dist] = knnsearch(centroids2, centroids1);

% Compute displacement
deltaX = centroids2(idx, 1) - centroids1(:, 1);
deltaY = centroids2(idx, 2) - centroids1(:, 2);

% Store results
displacement_table = table((1:numel(deltaX))', centroids1(:,1), centroids1(:,2), ...
                           centroids2(idx,1), centroids2(idx,2), ...
                           deltaX, deltaY, dist, ...
                           'VariableNames', {'ObjectID', 'XC_Ch1', 'YC_Ch1', ...
                           'XC_Ch2', 'YC_Ch2', 'DeltaX', 'DeltaY', 'Distance'});

% Display the results
disp(displacement_table);
cd('K:\new_codes\mobility_analysis\Uday\CytoD(1micM)')
displacment_correction = 'displacement_table';
save(displacment_correction, 'displacement_table');
cd(olddir);

figure;
quiver(centroids1(:,1), centroids1(:,2), deltaX, deltaY, 0, 'r');
hold on;
scatter(centroids1(:,1), centroids1(:,2), 'bo');
title('Displacement Vectors (Channel 1 → Channel 2)');
xlabel('X'); ylabel('Y');
axis equal;
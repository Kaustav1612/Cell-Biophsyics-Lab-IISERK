

function [all_stats,labeled_img_binary] = actin_positions(channel_2_dir,displacement_table,olddir)

channel_2_drift_X = displacement_table.XC_Ch1-displacement_table.XC_Ch2;
channel_2_drift_Y = displacement_table.YC_Ch1-displacement_table.YC_Ch2;

mean_channel_drift_X = mean(channel_2_drift_X);
mean_channel_drift_Y = mean(channel_2_drift_Y);

cd(channel_2_dir);
A=dir('*.tif');
nframes=size(A,1);
for i=1:nframes
    filename=A(i).name;
    image_stack_channel_2_uncorrected(:,:,i)=imread(filename);
end
cd(olddir);

num_frames = nframes;
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

image_stack_channel_2 = corrected_channel_2;

hold on;
redrawing = true;
while redrawing
    imshow(corrected_channel_2(:,:,1),[]);  % Display the black and white image
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

k=0;
for i = 1 :100 : nframes
    current_frame = image_stack_channel_2(:,:,i);

    % Filtering image using gaussian blur
    gauss_blur=imgaussfilt(current_frame,8,"FilterSize",[13 13]);
    filtered_image = current_frame-gauss_blur;

    % Normalisation of the filtered image
    min_value = min(filtered_image(:));  % Get minimum value of all elements
    max_value = max(filtered_image(:));  % Get maximum value of all elements
    normalized_image = double(filtered_image - min_value) / double(max_value - min_value);

    %Binarising using threshold
    threshold = 0.03;
    binary_image = imbinarize(normalized_image, threshold);

    % Erosion and Dilation
    eroded_image = imerode(binary_image, strel('disk',1));
    dilated_image=imdilate(eroded_image,strel('disk',1));
    BW = dilated_image;
    BW_1=dilated_image;

    BW(~mask) = 0;
    image = BW;
    masked_normalized_image = normalized_image.*mask;

    dilated_image(~mask) = 0;

    [componentIndices,stats_label,~] = detect(image,15,masked_normalized_image);
    skel_image = bwmorph(image, 'skel', Inf);


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

    labeled_image = zeros(size(image));
    labeled_img_binary =  zeros(size(image));
    labeled_image_norm = zeros(size(image));
    for j = 1 : numel(stats)
        Min_Area_Thresh = 20;
        Max_Area_Thresh = 5000;
        threshIntensity = 200;
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
        final_stats(u).PixelIdxList = pixel_indices;
    end
    k = k+1;
    % Save stats for current frame in an all_stats array
    all_stats{k} = final_stats;
    labeled_movie_norm_channel_2(:,:,k)=masked_normalized_image;
    detected_movie_channel_2(:,:,k) = image;
    skel_detected_movie_channel_2(:,:,k)= skel_image;

    figure();
    imshow(labeled_img_binary,[]);

end

end
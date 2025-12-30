function [Allstats,binary_stack,norm_stack] = ObjectDetector(I,BW_Bound,nframes,ThreshA,ThreshB,threshIntensity)
    Allstats = cell(nframes,1);
    
for i = 1:nframes
    I_frame = I(:,:,i);
    labeled_image = zeros(size(I_frame));
    labeled_image_norm = zeros(size(I_frame));
    % Apply Laplacian of Gaussian (LoG) filtering
    sigma=4;
    LoG = imgaussfilt(I_frame, sigma) - imgaussfilt(I_frame, sigma * 1.5);
    
    % Normalize LoG response and threshold using Otsu’s method
    normLoG = mat2gray(LoG);
    level = graythresh(normLoG);
    BW = imbinarize(normLoG, level);
    
    % Morphological operations
    BW = imfill(BW, 'holes');
    BW = imdilate(BW, strel('disk',1));
    BW_final = BW & BW_Bound;
    
    % Watershed segmentation to separate closely packed puncta
    D = -bwdist(~BW_final);
    D(~BW_final) = -Inf;
    L = watershed(D);
    BW_final(L == 0) = 0;
    
    % Gaussian Mixture Model (GMM) for intensity-based segmentation
    %pixels = double(I_frame(BW_final));
    %gm = fitgmdist(pixels, 2, 'RegularizationValue', 0.1);
    %BW_final = BW_final & (I_frame > gm.mu(1));
    
    % Detect connected components
    cc = bwconncomp(BW_final);
    stats = regionprops(cc,'Area','Centroid','PixelIdxList','Eccentricity','Solidity');

   
for u = 1:numel(stats)
    % Calculate and add 'Intensity'
    pixel_indices = stats(u).PixelIdxList;
    intensity = mean(I_frame(pixel_indices));
    stats(u).Intensity = intensity;

    % Extract and add 'XC' and 'YC' from the existing 'Centroid'
    % 'Centroid' is already an Nx2 array where N is the number of objects, and Centroid(u, :) is [x,y] for object u
    centroid = stats(u).Centroid;
    stats(u).XC = centroid(1); % X-coordinate of centroid
    stats(u).YC = centroid(2); % Y-coordinate of centroid

end


    idx = find(([stats.Area]>ThreshA) & ([stats.Area]<ThreshB) & ([stats.Intensity]>threshIntensity));
    statsFiltered = stats(idx);
    
    for j = 1:length(stats)
        currentArea = stats(j).Area;
        % Check if the area is within the specified thresholds
        if currentArea >= ThreshA && currentArea<ThreshB
            if stats(j).Intensity>threshIntensity
            % Access object pixels using PixelIdxList
                object_pixels = stats(j).PixelIdxList;
                labeled_image(object_pixels) = j;  % Label with object index
                labeled_image_norm(object_pixels) = I_frame(object_pixels);
            end
        end
    end
    
    Allstats{i} = statsFiltered;
    labeled_stack(:,:,i) = labeled_image;
    binary_stack(:,:,i)=labeled_image>0;
    norm_stack (:,:,i) = labeled_image_norm;
    
end

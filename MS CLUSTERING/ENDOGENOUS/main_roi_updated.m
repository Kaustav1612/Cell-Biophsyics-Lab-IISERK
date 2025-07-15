clear all;
clc;
close all;

%% Enlist the list of workspace saved variables

global n_objects  confocal_image cluster_centers  all_properties image  roi_wise_properties class_centers classify_features gfp_image;

%% Define the structures to be saved

all_properties=struct('Cluster',[],'Object',[],'Haralick_Features',[],'Haralick_Features_Scaled',[]);
roi_wise_properties =  struct('Cluster',[]);

img_path = "F:\Jibitesh\STED\14.12.23\Control\1\Image 1 - STAR RED STED.tiff";
confocal_path =  "F:\Jibitesh\STED\14.12.23\Control\1\Image 1 - STAR RED.tiff";
image = imread(img_path);
confocal_image = imread(confocal_path);

start_path = 'D:\Kaustav\CODES SET\MS CLUSTERING ALGO\Space_Shape_Intensity\Analysis\Results_May2025';
directory_overlayrois = uigetdir(start_path,'Overlay ROI directory');
directory_raw = uigetdir(start_path,'Raw ROI directory');
directory_label = uigetdir(start_path,'Labeled ROI directory');
directory_scatter = uigetdir(start_path,'Scatter plot directory');
directory_filtered_roi=uigetdir(start_path,'Filtered ROI directory');
%% Key parameters for object detection

threshold_distance = 60;
thresh = 3;

%% Masking and Automatic ROI generator

hold on;
redrawing = true;
userinput=[];
while redrawing
  imshow(image,[0 50]); 
  mask_handle = drawfreehand('Color','r','LineWidth',2);
  choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
  redrawing = strcmp(choice, "Yes");
  if redrawing
    delete(mask_handle);
  end
end

mask = mask_handle.createMask;
image(~mask) = 0;

incomplete_drawing = true;
while incomplete_drawing
    imshow(image, [0 50]);
    title("Fit ROI");
    hold on;
    redrawing = true;
    while redrawing
      h = drawfreehand('Color','r','LineWidth',2);
      choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
      redrawing = strcmp(choice, "Yes");
      if redrawing
        delete(h);
      end
    end

    BW = createMask(h);
    
    % Put the size of each ROI in pixels assuming a square ROI
    roi_size = 334;

    rois = {};
    normalized_rois={};
    number=0;
    all_rois={};
    all_normalized_rois={};
    all_confocal_rois={};
    image(~BW)=0;
    figure();
    title("Overlaid ROIs");
    imshow(image, [0 50]);
    [rows, cols] = size(BW);
    for i = 1:roi_size:rows
        for j = 1:roi_size:cols
            roi_y1 = i;
            roi_y2 = min(i + roi_size - 1, rows);
            roi_x1 = j;
            roi_x2 = min(j + roi_size - 1, cols);
    
            roi_mask = BW(roi_y1:roi_y2, roi_x1:roi_x2);
            if all(roi_mask(:))
                number =number+1;
                rectanglePoints = [roi_x1 roi_y1 (roi_size) (roi_size)];
                rectangle('Position', rectanglePoints, 'LineWidth', 1, 'EdgeColor', 'yellow');
                text(rectanglePoints(1)+(roi_size/3), rectanglePoints(2)+(roi_size/3), num2str(number), 'Color', 'yellow', 'FontSize', 8);
                norm_roi = image(roi_y1:roi_y2, roi_x1:roi_x2);
                confocal_roi=confocal_image(roi_y1:roi_y2, roi_x1:roi_x2);
                roi =  norm_roi > thresh;
                all_rois{end+1} = roi;
                all_normalized_rois{end+1} = norm_roi;
                all_confocal_rois{end+1} = confocal_roi;
            end
        end
    end
 choice_2 = questdlg("Are the ROIs covering the cells or you would like to redraw the mask?", "Freehand Mask", "Yes", "No", "Yes");
      
      
      incomplete_drawing = strcmp(choice_2, "Yes");
      if incomplete_drawing
        delete(h);
        close all;
      end
      
      if ~incomplete_drawing
        filename=sprintf('Endogenous.fig');
        savefig(fullfile(directory_overlayrois,filename))
      end
      
end


%% Choosing which ROI to be analysed 

% Step 1: Ask the user for the number of inputs
prompt1 = {'How many ROIs you want to analyse?'};
dlgTitle1 = 'No. of ROIs';
numLines1 = 1;
defaultAnswer1 = {'0'};

answer1 = inputdlg(prompt1, dlgTitle1, numLines1, defaultAnswer1);

% Check if the user clicked the OK button or canceled the dialog
if ~isempty(answer1)
    numInputs = str2double(answer1{1});
    
    if isnan(numInputs) || numInputs < 0
        error('Invalid number of inputs.');
    end
    
    if numInputs == 0
        rois=all_rois;
        for k = 1 : length(rois)
            userinput=[userinput,k]; 
            figure();
            imshow(all_normalized_rois{k},[]);
            filename=sprintf('STEDROI_%d.fig',k);
            savefig(fullfile(directory_raw,filename))
            close;
            figure();
            imshow(all_confocal_rois{k},[]);
            filename=sprintf('ConfocalROI_%d.fig',k);
            savefig(fullfile(directory_raw,filename))
            close;
        end
    else
        % Step 2: Ask the user to enter the inputs as a comma-separated list
        prompt2 = {sprintf('Enter %d ROIs separated by commas:', numInputs)};
        dlgTitle2 = '1,2,3,';
        numLines2 = 1;
        defaultAnswer2 = {'1,2,3,'};

        answer2 = inputdlg(prompt2, dlgTitle2, numLines2, defaultAnswer2);

        % Check if the user clicked the OK button or canceled the dialog
        if ~isempty(answer2)
            userInputsStr = answer2{1};

            % Split the input string by commas and convert to numbers
            userInputs = str2double(strsplit(userInputsStr, ','));
            userinput=userInputs;

            if length(userInputs) ~= numInputs
                error('Less number of ROIs given as input');
            end
            for k = 1 : length(userInputs)
                rois{end+1}=all_rois{userInputs(k)};
            end
            % Display the inputs
            disp('Analysing the following ROIs:');
            disp(userInputs);
        else
            disp('User canceled the input dialog for ROI.');
        end
    end
else
    disp('User canceled the input dialog for number of inputs.');
end%% Iterate over each of the selected ROIs to be analysed
   
for m =1 :size(all_rois,2)
    roi_norm = all_normalized_rois{m};
    confocal_roi = all_confocal_rois{m};
    [row_l, col_l] = size(all_rois{1,m});
    roi = all_rois{1,m};
    
    x_axis_pixel_count = zeros(1, row_l);
    y_axis_pixel_count = zeros(1, col_l);
    bright_pixel_count =  size(find(all_rois{1,m}),1);
    for i = 1: row_l
        count = 0;
        for j = 1: col_l
            if all_rois{1,m}(i,j) == 1
                count = count + 1;
            end
        end
        x_axis_pixel_count(1,i) = count;
    end

    for j = 1: col_l
        count = 0;
        for i = 1: row_l
            if all_rois{1,m}(i,j) == 1
                count = count + 1;
            end
        end
        y_axis_pixel_count(1,j) = count;
    end
    
    
    
    %% Analysis for STED and Confocal Image with varying alpha
  for alpha = 1 : 5
         fprintf('Current alpha value = %d \n',alpha);

        % Analysis for STED ROI
        
        [filtered_img_binary_sted,filtered_img_norm_sted]= filtration(roi,alpha,bright_pixel_count,threshold_distance,roi_norm);
        figure();
        imshow(filtered_img_binary_sted,[]);
        title(sprintf('Binary Image STED(alpha=%d)',alpha));
        close;

        eroded_image_sted = imerode(filtered_img_binary_sted, strel('disk',1));
        dilated_image_sted=imdilate(eroded_image_sted,strel('disk',1));
        BW_sted = dilated_image_sted;
        
        eroded_image_sted_norm = imerode(filtered_img_norm_sted, strel('disk',1));
        dilated_image_sted_norm=imdilate(eroded_image_sted_norm,strel('disk',1));
        BW_sted_norm = dilated_image_sted_norm;

        figure();
        imshow(BW_sted_norm,[]);
        title(sprintf('Filtered Transfected STED ROI (alpha=%d)',alpha));
        filename=sprintf('STED_Transfected_Clusters,alpha=%d,roi=%d.fig',alpha,m);
        savefig(fullfile(directory_filtered_roi,filename))
        close;
        

   
    conf_yes = 0;
   [number_obj, componentIndices, stats,componentAreas,stats_label,roi_LabeledImage,centroids] = detect(BW_sted,6,BW_sted_norm,conf_yes);
   if number_obj > 0
       figure();
       imshow(roi_LabeledImage,[0 number_obj],'Colormap',jet);
       colorbar;
       title(sprintf('Initial Labeled Image (alpha=%d)',alpha));
       filename=sprintf('NC_STED_ROI_%d,alpha=%d.fig',m,alpha);
       savefig(fullfile(directory_label,filename));
       close;

       params.bandwidth = 2; 
       params.epsilon = params.bandwidth*0.01; % Convergence Checker Parameter for shifting of the cluster center until it reaches the attraction basin
       params.learning_parameter = 0.1; % Rate of Learning
       params.merge_distance = 0.05*params.bandwidth;
       [cluster_centers,cluster_assignments,initial_object_positions,final_object_positions]=cluster_object_space(BW_sted,2,filtered_img_norm_sted,params);
    n_objects = size(cluster_centers,1);
   
    
   roi_clusteredLabeledImage = zeros(size(BW_sted));
   for i = 1 : size(initial_object_positions,1)
       roi_clusteredLabeledImage(initial_object_positions(i,1),initial_object_positions(i,2)) = cluster_assignments(i);
   end
   initial_label_sted = [];
   final_label_sted= [];
   for k = 1 : size(cluster_centers,1)
       initial_label_sted = [initial_label_sted;roi_LabeledImage(round(cluster_centers(k,1)),round(cluster_centers(k,2)))];
       final_label_sted = [final_label_sted;roi_clusteredLabeledImage(round(cluster_centers(k,1)),round(cluster_centers(k,2)))];
   end
    final_label_sted = nonzeros(final_label_sted);
    % Get unique cluster assignments (excluding zeros if present)
    unique_values = nonzeros(unique(initial_label_sted)); % Remove zero if exists

    if ~isempty(unique_values)
        % Create bin edges for histcounts
        bin_edges = [unique_values; max(unique_values)+1]'; % Transpose to row vector

        % Count occurrences of each cluster
        counts_sted = histcounts(initial_label_sted, bin_edges);

        % Create figure - only plot clusters with counts > 0
        figure();

        % Filter out zero counts (though histcounts shouldn't produce them)
        valid_counts = counts_sted(counts_sted > 0);
        valid_labels = unique_values(counts_sted > 0);
        figure();
        histogram(valid_counts);
        filename=sprintf('Cluster_Dist.fig',m,alpha);
        savefig(fullfile(directory_label,filename));
        close;
    else
        warning('No valid clusters found (only zeros or empty array)');
    end

   [object_properties, cluster_properties, roi_cluster_properties,distances]=enlist_properties_updated(cluster_centers,initial_object_positions,roi_clusteredLabeledImage,roi_LabeledImage,filtered_img_norm_sted,stats_label,componentIndices);
    all_properties.Cluster=[all_properties.Cluster,cluster_properties];
    roi_wise_properties.Cluster = [roi_wise_properties.Cluster,roi_cluster_properties];
    all_properties.Object=[all_properties.Object,object_properties];
    
    figure();
    imshow(roi_clusteredLabeledImage,[0 size(cluster_centers,1)],'Colormap',jet);
    title(sprintf('Clustered Labeled Image (alpha=%d)',alpha));
    colorbar;
    filename=sprintf('C_STED_ROI_%d,alpha=%d.fig',m,alpha);
    savefig(fullfile(directory_label,filename));
    close;
  
    
    figure();
    gscatter(final_object_positions(:,1),final_object_positions(:,2),cluster_assignments);
    hold on;
    plot(cluster_centers(:,1),cluster_centers(:,2),'kx','MarkerSize', 3, 'LineWidth', 3);
    title(sprintf('Clustered Labeled Plot (alpha=%d)',alpha));
    filename=sprintf('C_STED_Plot_ROI_%d,alpha=%d.fig',m,alpha);
    savefig(fullfile(directory_scatter,filename));
    close;
   

    figure();
    gscatter(initial_object_positions(:,1),initial_object_positions(:,2),cluster_assignments);
    hold on;
    title(sprintf('Initial Labeled Plot (alpha=%d)',alpha));
    filename=sprintf('NC_STED_Plot_ROI_%d,alpha=%d.fig',m,alpha);
    savefig(fullfile(directory_scatter,filename)); 
    close;
   end 

   % For Confocal Analysis

   [filtered_img_binary_conf,filtered_img_norm_conf] = filtration(confocal_roi,alpha,bright_pixel_count,threshold_distance,roi_norm);
   figure();
   imshow(filtered_img_binary_conf,[]);
   title(sprintf('Binary Image Confocal (alpha=%d)',alpha));
   close;   
    
    eroded_image_conf = imerode(filtered_img_binary_conf, strel('disk',1));
    dilated_image_conf=imdilate(eroded_image_conf,strel('disk',1));
    BW_conf = dilated_image_conf;
    
    eroded_image_conf_norm = imerode(filtered_img_norm_conf, strel('disk',1));
    dilated_image_conf_norm=imdilate(eroded_image_conf_norm,strel('disk',1));
    BW_conf_norm = dilated_image_conf_norm;

        figure();
        imshow(BW_conf_norm,[]);
        title(sprintf('Filtered Transfected Confocal ROI (alpha=%d)',alpha));
        filename=sprintf('Confocal_Transfected_Clusters,alpha=%d.fig',alpha);
        savefig(fullfile(directory_filtered_roi,filename))
        close;
    
     conf_yes = 1;
   [number_obj, componentIndices, stats,componentAreas,stats_label,roi_LabeledImage,centroids] = detect(BW_conf,6,BW_conf_norm,conf_yes);
   if number_obj > 0
   figure();
   imshow(roi_LabeledImage,[0 number_obj],'Colormap',jet);
   colorbar;
   title(sprintf('Initial Labeled Image (alpha=%d)',alpha));
   filename=sprintf('NC_Confocal_ROI_%d,alpha=%d.fig',m,alpha);
   savefig(fullfile(directory_label,filename));
   close;
   end
   
   
%    [class_centers, class_assignments, classify_features,all_haralick_features,all_haralick_features_scaled] = cluster_object_shape(BW, 1.5, 2, filtered_img_norm, roi_clusteredLabeledImage, cluster_assignments, n_objects);
%     n_classes = size(class_centers, 1);
%     
% % Printing the cluster centers 
% for i = 1:n_classes
%     class_center = class_centers(i, :);
%     fprintf('Cluster %d: [', i);
%     fprintf('%f ', class_center);  % Print each feature value with a space
%     fprintf(']\n');
% end
% 
% 
% % Create the labeled image 
% classLabeledImage = zeros(size((all_rois{1,m})));
% for i = 1:size(class_assignments, 1)
%     classLabeledImage(roi_clusteredLabeledImage == i) = class_assignments(i);
% end
% 
% 
%     all_properties.Haralick_Features=[all_properties.Haralick_Features;all_haralick_features];
%     all_properties.Haralick_Features_Scaled=[all_properties.Haralick_Features_Scaled;all_haralick_features_scaled];
% 
%    
%     figure();
%     imshow(classLabeledImage,[0 size(class_centers,1)],'Colormap',jet);
%     close;
%    
    
   
  end
 end



save(['Properties_NonTransCntrl_Cell_01','.mat'],'all_properties','roi_wise_properties');

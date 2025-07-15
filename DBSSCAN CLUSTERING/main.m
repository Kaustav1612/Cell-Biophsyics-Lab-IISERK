clc;
close all;
clear all;

%% Enlist the list of workspace saved variables

global n_objects cluster_centers  all_properties image  roi_wise_properties gfp_image confocal_image valid_counts_all;

%% Define the structures to be saved

all_properties=struct('Cluster',[],'Object',[],'Haralick_Features',[],'Haralick_Features_Scaled',[]);
roi_wise_properties =  struct('Cluster',[]);

%% Path to the tiff file

img_path = "H:\STED\Cont\image8\Image 10 - STAR RED STED.tiff";
gfp_path = "H:\STED\Cont\image8\Image 10 - eGFP.tiff";
confocal_path =  "H:\STED\Cont\image8\Image 10 - STAR RED.tiff";
image = imread(img_path);
gfp_image = imread(gfp_path);
confocal_image = imread(confocal_path);
% start_path = 'D:\Kaustav\CODES SET\MS CLUSTERING ALGO\Space_Shape_Intensity\Analysis\Results_May2025';
% directory_overlayrois = uigetdir(start_path,'Overlay ROI directory');
% directory_raw = uigetdir(start_path,'Raw ROI directory');
% directory_label = uigetdir(start_path,'Labeled ROI directory');
% directory_scatter = uigetdir(start_path,'Scatter plot directory');
% directory_filtered_roi=uigetdir(start_path,'Filtered ROI directory');
%% Key parameters for object detection

threshold_distance = 60;
thresh = 3;

%% Masking and Automatic ROI generator

hold on;
redrawing = true;
userinput=[];
while redrawing
  imshow(image,[]);
  title("Select and Crop out the Cell");
   
  rect = getrect;
  choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
  if strcmp(choice, "No")
      redrawing=0;
  else
      redrawing=1;
  end
  
  choice_ROI = questdlg('Do you want to keep this ROI or redraw?', ...
                          'Automatic Masked', ...
                          'Automatic Masked', 'Manual', 'Manual');
                      
        % Manual Rectangular ROI

        if strcmp(choice_ROI, 'Manual')
            close all;
            image = imcrop(image, rect);
            gfp_image = imcrop(gfp_image, rect);
            confocal_image = imcrop(confocal_image,rect);
            [roi_coords,sideLength] = draw_multiple_square_rois(image, gfp_image);
            close all;
            % Initialize containers for ROIs and normalized ROIs
            all_rois = {};
            all_normalized_rois = {};
            all_gfp_rois={};
            all_confocal_rois={};

            for i = 1:length(roi_coords)
                % Get the square ROI coordinates
                rectanglePoints = roi_coords{i};

                % Extract top-left and bottom-right corners
                min_row = round(rectanglePoints(1, 2)); % Top-left Y (row)
                max_row = round(rectanglePoints(3, 2)); % Bottom-right Y (row)
                min_col = round(rectanglePoints(1, 1)); % Top-left X (column)
                max_col = round(rectanglePoints(2, 1)); % Bottom-right X (column)

                % Validate indices to stay within image bounds
                min_row = max(1, min_row);
                max_row = min(size(image, 1), max_row);
                min_col = max(1, min_col);
                max_col = min(size(image, 2), max_col);


                % Extract corresponding region in the image
                norm_roi = image(min_row:max_row, min_col:max_col);
                figure();
                imshow(norm_roi,[]);
                roi =  norm_roi > thresh;
                close;
                gfp_roi = gfp_image(min_row:max_row, min_col:max_col);
                confocal_roi = confocal_image(min_row:max_row, min_col:max_col);
                % Store the ROI and normalized ROI
                all_rois{end+1} = roi; 
                all_normalized_rois{end+1} = norm_roi;
                all_gfp_rois{end+1} = gfp_roi;
                all_confocal_rois{end+1} = confocal_roi;
                
                gfp_roi = gfp_roi> thresh;
                confocal_roi = confocal_roi>thresh;
                
            end

                  for k = 1 : length(all_rois)
                          figure();
                          imshow(all_normalized_rois{k},[]);
                          filename=sprintf('STEDROI_%d.fig',k);
%                           savefig(fullfile(directory_raw,filename))
%                           close;
                          figure();
                          imshow(all_confocal_rois{k},[]);
                          filename=sprintf('ConfocalROI_%d.fig',k);
%                           savefig(fullfile(directory_raw,filename))
%                           close;
                  end

        % Automatic  Rectangular ROI
        elseif strcmp(choice_ROI, 'Automatic Masked')
        gfp_image = imread(gfp_path);
        image = imcrop(image, rect);
        gfp_image = imcrop(gfp_image, rect);
        incomplete_drawing = true;
        while incomplete_drawing
            imshow(image, []);
            hold on;
            redrawing = true;
            while redrawing
              % Allow the user to draw a freehand boundary
              h = drawfreehand('Color','r','LineWidth',2);

              % Option to redraw or exit
              choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");

              % Update flag based on user choice
              redrawing = strcmp(choice, "Yes");

              % Clear the previous mask for redraw
              if redrawing
                delete(h);
              end
            end

            % Create a binary mask from the freehand boundary
            BW = createMask(h);

            % Size of the rectangular ROI
            roi_size = 334;

            % Initialize a cell array to store ROIs
            rois = {};
            normalized_rois={};
            number=0;
            all_rois={};
            all_normalized_rois={};
            % Loop through the mask and extract ROIs
            [rows, cols] = size(BW);
            for i = 1:roi_size:rows
                for j = 1:roi_size:cols
                    % Define the ROI boundaries
                    roi_y1 = i;
                    roi_y2 = min(i + roi_size - 1, rows);
                    roi_x1 = j;
                    roi_x2 = min(j + roi_size - 1, cols);

                    % Check if the entire ROI is within the mask
                    roi_mask = BW(roi_y1:roi_y2, roi_x1:roi_x2);
                    if all(roi_mask(:))
                        number =number+1;
                        rectanglePoints = [roi_x1 roi_y1 (roi_size) (roi_size)];
                        rectangle('Position', rectanglePoints, 'LineWidth', 1, 'EdgeColor', 'yellow');
                        text(rectanglePoints(1)+(roi_size/3), rectanglePoints(2)+(roi_size/3), num2str(number), 'Color', 'yellow', 'FontSize', 8);
                      
                        
                        norm_roi = image(roi_y1:roi_y2, roi_x1:roi_x2);
                        roi =  norm_roi > thresh;
                        all_rois{end+1} = roi;
                        all_normalized_rois{end+1} = norm_roi;
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
                filename=sprintf('RoiOverlay_Transfected_ZAFCntrl_Cell_%d.fig',05);
%                 savefig(fullfile(directory_overlayrois,filename))
%                 close;
              end

        end

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
                                imshow(all_rois{k},[]);
                                title('Roi');
                                directory='D:\Kaustav\CODES SET\MS CLUSTERING ALGO\Spatial Clustering\Analysis\Process_Images\RAW EXTRACTED ROIS\Transfected\ZAF\CELL_06';
                                filename=sprintf('ROI_%d.fig',k);
%                                 savefig(fullfile(directory,filename))
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
                        end
                    end
                end
        end
      
end



%% Choosing which ROI to be analysed 
valid_counts_all= [];
if strcmp(choice_ROI, 'Manual')   
    for m =1 :size(all_rois,2)
        roi_size = sideLength{m};
        roi_norm = all_normalized_rois{m};
        gfp_roi = all_gfp_rois{m};
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
        
      gfp_roi_bin=gfp_roi>0;
      for alpha = 1 : 3
             fprintf('Current alpha value = %d \n',alpha);

            % Analysis for STED ROI

            [filtered_img_binary_sted,filtered_img_norm_sted]= filtration(roi,alpha,bright_pixel_count,threshold_distance,roi_norm);
            figure();
            imshow(filtered_img_binary_sted,[]);
            title(sprintf('Binary Image STED(alpha=%d)',alpha));
            close;

            if isequal(size(filtered_img_binary_sted),size(gfp_roi_bin))
               boolean_out_sted = (filtered_img_binary_sted == 1) & (gfp_roi_bin == 1);
            end

            filtered_img_binary_sted = filtered_img_binary_sted & boolean_out_sted; % Proper binary output
            filtered_img_norm_sted = filtered_img_norm_sted .* boolean_out_sted; % Preserves original values

            eroded_image_sted = imerode(filtered_img_binary_sted, strel('disk',1));
            dilated_image_sted=imdilate(eroded_image_sted,strel('disk',1));
            BW_sted = dilated_image_sted;

            eroded_image_sted_norm = imerode(filtered_img_norm_sted, strel('disk',1));
            dilated_image_sted_norm=imdilate(eroded_image_sted_norm,strel('disk',1));
            BW_sted_norm = dilated_image_sted_norm;

            figure();
            imshow(BW_sted_norm,[]);
            title(sprintf('Filtered Transfected STED ROI (alpha=%d)',alpha));
            filename=sprintf('STED_Transfected_Clusters,alpha=%d.fig',alpha);
%             savefig(fullfile(directory_filtered_roi,filename))
%             close;



       conf_yes=0;
       [number_obj, componentIndices, stats,componentAreas,stats_label,roi_LabeledImage,centroids] = detect(BW_sted,6,BW_sted_norm,conf_yes);

        [features,~]=feature_extraction_space(roi,filtered_img_norm_sted);
if number_obj > 0
    % Extract centroid coordinates
    X = features; % Nx2 matrix of [x,y] coordinates
    
    % Set DBSCAN parameters
    epsilon = 5; % Use your existing threshold or adjust
    minpts = 1; % Start with this, can be adjusted
    
    % Perform DBSCAN clustering
    labels = dbscan(X, epsilon, minpts);
    
    % Visualize results
    figure;
    gscatter(X(:,1), X(:,2), labels);
    title(sprintf('DBSCAN Clustering (alpha=%d, ε=%d, minpts=%d)', alpha, epsilon, minpts));
    xlabel('X position');
    ylabel('Y position');
%     filename = sprintf('DBSCAN_STED_ROI_%d_alpha=%d.fig', m, alpha);
%     savefig(fullfile(directory_scatter, filename));
%     close;
    
    % Store cluster information
    n_clusters = max(labels); % -1 is noise
    cluster_sizes = arrayfun(@(x) sum(labels==x), unique(labels));
    
    % You can add this to your properties structure
    all_properties(m).DBSCAN_Clusters = n_clusters;
    all_properties(m).DBSCAN_ClusterSizes = cluster_sizes;
    all_properties(m).DBSCAN_Labels = labels;
    all_properties(m).DBSCAN_Centroids = X;
end
% After DBSCAN clustering
if number_obj > 0
    % Create a labeled image matrix
    cluster_labeled_image = zeros(size(BW_sted));
    for i = 1 : size(X,1)
           cluster_labeled_image(X(i,1),X(i,2)) = labels(i,1);
    end
    
    
    figure;
    imshow(roi_LabeledImage, [0 number_obj], 'Colormap', jet);
    colorbar;
    title(sprintf('Non-Clustered-labeled ROI (α=%d, %d clusters)', alpha, number_obj));
    
    % Display the cluster-labeled image
    figure;
    imshow(cluster_labeled_image, [0 max(labels)], 'Colormap', jet(max(labels)+1));
    colorbar;
    title(sprintf('Cluster-labeled ROI (α=%d, %d clusters)', alpha, max(labels)));
    
    % Save the figure
    filename = sprintf('ClusterLabeled_ROI%d_alpha%d.fig', m, alpha);
%     savefig(fullfile(directory_label, filename));
    
    % Also save as PNG
    filename_png = sprintf('ClusterLabeled_ROI%d_alpha%d.png', m, alpha);
%     saveas(gcf, fullfile(directory_label, filename_png));
%     close;
%     
    % Store in your properties structure
    all_properties(m).ClusterLabeledImage = cluster_labeled_image;
    end
   end
end
    
elseif strcmp(choice_ROI, 'Automatic Masked')
  for m =1 :size(all_rois,2)
    roi_size = sideLength{m};
    roi_norm = all_normalized_rois{m};
    gfp_roi = all_gfp_rois{m};
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
gfp_roi_bin=gfp_roi>0;
        for alpha = 1 : 3
             fprintf('Current alpha value = %d \n',alpha);

            % Analysis for STED ROI

            [filtered_img_binary_sted,filtered_img_norm_sted]= filtration(roi,alpha,bright_pixel_count,threshold_distance,roi_norm);
            figure();
            imshow(filtered_img_binary_sted,[]);
            title(sprintf('Binary Image STED(alpha=%d)',alpha));
            close;

            if isequal(size(filtered_img_binary_sted),size(gfp_roi_bin))
               boolean_out_sted = (filtered_img_binary_sted == 1) & (gfp_roi_bin == 1);
            end

            filtered_img_binary_sted = filtered_img_binary_sted & boolean_out_sted; % Proper binary output
            filtered_img_norm_sted = filtered_img_norm_sted .* boolean_out_sted; % Preserves original values

            eroded_image_sted = imerode(filtered_img_binary_sted, strel('disk',1));
            dilated_image_sted=imdilate(eroded_image_sted,strel('disk',1));
            BW_sted = dilated_image_sted;

            eroded_image_sted_norm = imerode(filtered_img_norm_sted, strel('disk',1));
            dilated_image_sted_norm=imdilate(eroded_image_sted_norm,strel('disk',1));
            BW_sted_norm = dilated_image_sted_norm;

            figure();
            imshow(BW_sted_norm,[]);
            title(sprintf('Filtered Transfected STED ROI (alpha=%d)',alpha));
            filename=sprintf('STED_Transfected_Clusters,alpha=%d.fig',alpha);
%             savefig(fullfile(directory_filtered_roi,filename))
%             close;



       conf_yes=0;
       [number_obj, componentIndices, stats,componentAreas,stats_label,roi_LabeledImage,centroids] = detect(BW_sted,6,BW_sted_norm,conf_yes);

        [features,~]=feature_extraction_space(roi,filtered_img_norm_sted);
if number_obj > 0
    % Extract centroid coordinates
    X = features; % Nx2 matrix of [x,y] coordinates
    
    % Set DBSCAN parameters
    epsilon = 5; % Use your existing threshold or adjust
    minpts = 1; % Start with this, can be adjusted
    
    % Perform DBSCAN clustering
    labels = dbscan(X, epsilon, minpts);
    
    % Visualize results
    figure;
    gscatter(X(:,1), X(:,2), labels);
    title(sprintf('DBSCAN Clustering (alpha=%d, ε=%d, minpts=%d)', alpha, epsilon, minpts));
    xlabel('X position');
    ylabel('Y position');
%     filename = sprintf('DBSCAN_STED_ROI_%d_alpha=%d.fig', m, alpha);
%     savefig(fullfile(directory_scatter, filename));
%     close;
    
    % Store cluster information
    n_clusters = max(labels); % -1 is noise
    cluster_sizes = arrayfun(@(x) sum(labels==x), unique(labels));
    
    % You can add this to your properties structure
    all_properties(m).DBSCAN_Clusters = n_clusters;
    all_properties(m).DBSCAN_ClusterSizes = cluster_sizes;
    all_properties(m).DBSCAN_Labels = labels;
    all_properties(m).DBSCAN_Centroids = X;
end
% After DBSCAN clustering
if number_obj > 0
    % Create a labeled image matrix
    cluster_labeled_image = zeros(size(BW_sted));
    for i = 1 : size(X,1)
           cluster_labeled_image(X(i,1),X(i,2)) = labels(i,1);
    end
    
    
    figure;
    imshow(roi_LabeledImage, [0 number_obj], 'Colormap', jet);
    colorbar;
    title(sprintf('Non-Clustered-labeled ROI (α=%d, %d clusters)', alpha, number_obj));
    
    % Display the cluster-labeled image
    figure;
    imshow(cluster_labeled_image, [0 max(labels)], 'Colormap', jet(max(labels)+1));
    colorbar;
    title(sprintf('Cluster-labeled ROI (α=%d, %d clusters)', alpha, max(labels)));
    
    % Save the figure
    filename = sprintf('ClusterLabeled_ROI%d_alpha%d.fig', m, alpha);
%     savefig(fullfile(directory_label, filename));
    
    % Also save as PNG
    filename_png = sprintf('ClusterLabeled_ROI%d_alpha%d.png', m, alpha);
%     saveas(gcf, fullfile(directory_label, filename_png));
%     close;
%     
    % Store in your properties structure
    all_properties(m).ClusterLabeledImage = cluster_labeled_image;

   end  
        end
  end
end


    
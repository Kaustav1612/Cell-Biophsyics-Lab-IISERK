clc;
close all;
clear all;

%% Path to the tiff file

img_path = "I:\Jibitesh\STED\11.12.2023\Control\transfected\image16\Image 18 - STAR RED STED.tiff";
gfp_path = "I:\Jibitesh\STED\11.12.2023\Control\transfected\image16\Image 18 - eGFP.tiff";
confocal_path =  "I:\Jibitesh\STED\11.12.2023\Control\transfected\image16\Image 18 - STAR RED.tiff";
directory = 'C:\Users\Bidisha Sinha\Downloads\drive-download-20250703T063204Z-1-001\sim_results\hill';
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  image = imread(img_path);
gfp_image = imread(gfp_path);
confocal_image = imread(confocal_path);

%% Key parameters for object detection

threshold_distance = 60;
thresh_sted = 3;
thresh_confocal = 15;

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
                roi =  norm_roi > thresh_sted;
                close;
                gfp_roi = gfp_image(min_row:max_row, min_col:max_col);
                confocal_roi = confocal_image(min_row:max_row, min_col:max_col);
                % Store the ROI and normalized ROI
                all_rois{end+1} = roi; 
                all_normalized_rois{end+1} = norm_roi;
                all_gfp_rois{end+1} = gfp_roi;
                confocal_roi = confocal_roi>thresh_confocal;
                all_confocal_rois{end+1} = confocal_roi;
                
                gfp_roi = gfp_roi> thresh_sted;

                
            end

            for k = 1 : length(all_rois)
                 figure();
                 imshow(all_normalized_rois{k},[]);

                 figure();
                 imshow(confocal_roi,[]);

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


        
            gfp_roi_bin=gfp_roi>0;
            alpha = 3;


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
            filename=sprintf('ROI=%d,STED_Transfected_Clusters,alpha=%d.fig',m,alpha);
            savefig(fullfile(directory, filename));
            close;


       % For Confocal Analysis
       filtered_img_norm_conf = zeros(size(confocal_roi));
       %[filtered_img_binary_conf,filtered_img_norm_conf] = filtration(confocal_roi>0,alpha,bright_pixel_count,threshold_distance,confocal_roi);
       filtered_img_binary_conf = confocal_roi>0;
       filtered_img_norm_conf_idx = find(confocal_roi>0);
       filtered_img_norm_conf(filtered_img_norm_conf_idx) = confocal_roi(filtered_img_norm_conf_idx);
       % figure();
       % imshow(confocal_roi>0,[]);
       % title(sprintf('Binary Image Confocal (alpha=%d)',alpha));



       if isequal(size(filtered_img_binary_sted),size(gfp_roi_bin))
          boolean_out_conf = (filtered_img_binary_conf == 1) & (gfp_roi_bin == 1);
       end


        filtered_img_binary_conf = filtered_img_binary_conf & boolean_out_conf; % Proper binary output
        filtered_img_norm_conf = filtered_img_norm_conf .* boolean_out_conf; % Preserves original values
 
        eroded_image_conf = imerode(filtered_img_binary_conf, strel('disk',1));
        dilated_image_conf=imdilate(eroded_image_conf,strel('disk',1));
        BW_conf = dilated_image_conf;

        eroded_image_conf_norm = imerode(filtered_img_norm_conf, strel('disk',1));
        dilated_image_conf_norm=imdilate(eroded_image_conf_norm,strel('disk',1));
        BW_conf_norm = dilated_image_conf_norm;

       figure();
       imshow(filtered_img_norm_conf,[]);
       title(sprintf('Filtered Transfected Confocal ROI (alpha=%d)',alpha));
       filename=sprintf('Confocal ROI = %d (alpha=%d).fig',m,alpha);
       savefig(fullfile(directory, filename));
       close;

       conf_yes=1;
       [number_obj, componentIndices, stats,componentAreas,stats_label,roi_LabeledImage,centroids,BW_conf] = detect(BW_conf,6,BW_conf_norm,conf_yes);
       if number_obj > 0
       figure();
       imshow(roi_LabeledImage,[0 number_obj]);
       colorbar;
       title(sprintf('Initial Labeled Image (alpha=%d)',alpha));

       results = monte_carlo(BW_conf,stats,m); 
       sprintf("Roi No %d completed",m);  

       end             
    end       
   end
  
%%     


function [results] = monte_carlo(BW_conf,stats,m)
%% 1. Initialize Parameters
num_simulations = 20;
antibody_molecules = linspace(1.0e-8, 1.0e-7, 3 );

image_size_pixels = size(BW_conf);
num_receptors = length(find(BW_conf>0));
pixel_idx_allowed = [];

for k = 1 :length(stats)
    pixel_idx_allowed = [pixel_idx_allowed;stats(k).PixelIdxList];
end

K_d = linspace(1.0e-8, 5.0e-8, 3) + 0.2e-7 * abs(randn(1,3)); % Dissociation constant in M (moles/L)

% Physical parameters
pixel_size_m = 2.0e-8; % Pixel size in meters (m)

% Convert pixel size to micrometers for volume calculation clarity
pixel_size_um = pixel_size_m * 1e6; % 1 m = 1e6 um
% pixel_size_um will be 0.02 um

image_area_um2 = image_size_pixels(1) * image_size_pixels(2) * (pixel_size_um)^2; % Area in um^2
cell_depth_um = 8; % Cell depth in micrometers (um)
volume_um3 = image_area_um2 * cell_depth_um * 1e-3; % Volume in um^3
volume_pL = volume_um3 * 1e-3; % Convert um^3 to pL (1 pL = 1000 um^3)
volume_ml = volume_pL * 1e-9; % Convert pL to mL (1 mL = 1e9 pL)

% Biological parameters
antibody_MW_kDa = 100; % Molecular Weight of antibody in kDa
% Mass of one antibody molecule in nanograms (1 Da = 1.660539e-24 g, 1 g = 1e9 ng)
mass_per_antibody_ug = antibody_MW_kDa * 1000 * 1.660539e-24 * 1e6;
% Simplified: mass_per_antibody_ng = antibody_MW_kDa * 1.660539e-12;
% mass_per_antibody_ng = 100 * 1.660539e-12 = 1.660539e-10 ng 
conc_antibody = (antibody_molecules .* mass_per_antibody_ug) ./ volume_ml;
directory = 'C:\Users\Bidisha Sinha\Downloads\drive-download-20250703T063204Z-1-001\sim_results\hill';

results = antibody_staining_single_cell(num_simulations, conc_antibody, num_receptors, image_size_pixels, K_d,pixel_idx_allowed,directory,BW_conf,m);
end
function [results] = antibody_staining_single_cell(num_simulations, antibody_concs, num_receptors, image_size, K_d, pixel_idx_allowed,directory,BW_conf,m)
    
    % 2. Output Storage
    results.intensities = zeros(length(antibody_concs), num_simulations, length(K_d));
    results.cv = zeros(length(antibody_concs), 1, length(K_d));
    results.snr = zeros(length(antibody_concs), 1, length(K_d));
    results.cluster_sizes = zeros(length(antibody_concs), num_simulations, length(K_d));
    results.cluster_counts = zeros(length(antibody_concs), num_simulations, length(K_d));
    results.bound_number = zeros(length(antibody_concs), num_simulations, length(K_d));
    results.cluster_radii = zeros(length(antibody_concs), num_simulations, length(K_d));

    % 3. Main Simulation Loop
    for K_d_idx = 1:length(K_d)
        current_K_d = K_d(K_d_idx);
        for conc_idx = 1:length(antibody_concs)
            current_conc = antibody_concs(conc_idx);
            conc_results = zeros(num_simulations, 1);
            cluster_sizes = zeros(num_simulations, 1);
            cluster_counts = zeros(num_simulations, 1);
            cluster_radius_all = zeros(num_simulations, 1);
            bound_number = zeros(num_simulations, 1);
            rand_sim = randi([1, num_simulations]);
           
            parfor sim = 1:num_simulations
                % A. Generate receptor coordinates
                receptor_map = generate_single_cell_receptors(image_size, num_receptors,pixel_idx_allowed);

                % B. Simulate staining
                [base_image, N_bound] = simulate_staining(receptor_map, current_conc, image_size, current_K_d);
                
                % C. Apply features and PSF
                synthetic_image = apply_advanced_features(base_image, current_conc, image_size);

                % D. Cluster Analysis
                non_zero_pixels = synthetic_image(synthetic_image > 0);
                if ~isempty(non_zero_pixels)
                    threshold = prctile(non_zero_pixels, 95) * (1 + 0.1*log10(current_conc)); % Log-scale adjustment
                else
                    threshold = 0;
                end
                imshow(synthetic_image,[]);
                binary_image = synthetic_image > threshold;
                eroded_image = imerode(binary_image, strel('disk',1));
                dilated_image = imdilate(eroded_image, strel('disk',1));
                binary_image = dilated_image;
                labeled_clusters = bwlabel(binary_image);
                
                props = regionprops(labeled_clusters, 'Area', 'EquivDiameter'); 
                cluster_areas = [props.Area];
                cluster_radius = [props.EquivDiameter];
                cluster_areas = cluster_areas(cluster_areas > 6);
                cluster_radius = cluster_radius(cluster_areas > 6) .* 15;

                conc_results(sim) = mean(synthetic_image(synthetic_image > 0));
                cluster_sizes(sim) = mean(cluster_areas);
                cluster_radius_all(sim) = mean(cluster_radius);
                cluster_counts(sim) = numel(cluster_areas);
                bound_number(sim) = N_bound;
                
                if conc_idx == 1 || conc_idx == length(antibody_concs) || mod(conc_idx,3) == 0
                    visualize_results(synthetic_image, current_conc, sim, BW_conf, rand_sim, current_K_d,directory,m);
                end
            end

            results.intensities(conc_idx, :, K_d_idx) = conc_results;
            results.cv(conc_idx, K_d_idx) = std(conc_results)/mean(conc_results);
            results.snr(conc_idx, K_d_idx) = mean(conc_results)/std(conc_results);
            results.cluster_sizes(conc_idx, :, K_d_idx) = cluster_sizes;
            results.cluster_counts(conc_idx, :, K_d_idx) = cluster_counts;
            results.bound_number(conc_idx, :, K_d_idx) = bound_number;
            results.cluster_radii(conc_idx, :, K_d_idx) = cluster_radius_all;
        end
    end
    analyze_results(results, antibody_concs, K_d,directory,m);
end


%% Generate random receptor locations within a cell
function receptor_map = generate_single_cell_receptors(image_size, num_receptors,pixel_idx_allowed)
    receptor_map = zeros(image_size);
    for i = 1:num_receptors
          receptor_map(pixel_idx_allowed) = receptor_map(pixel_idx_allowed) + 1;
    end
end

%% Simulate antibody staining
function [base_image, N_bound] = simulate_staining(receptor_map, current_conc, image_size, K_d)
    % For random binding
    % binding_prob = rand();
    % For Hill's binding conditions
    binding_prob = (current_conc)^2/((current_conc)^2+(K_d)^2);
    base_image = zeros(image_size);
    [y, x] = find(receptor_map > 0);
    N_bound = 0;
    for i = 1:length(y)
        if rand() < binding_prob
            chosenNumber = chooseNumberSimple(binding_prob);
            base_image(y(i), x(i)) = chosenNumber; 
            N_bound = N_bound + 1;
        end
    end
    non_specific = 0.01 * current_conc * (rand(image_size) < 0.05);
    base_image = base_image + non_specific;
end

%% Apply microscope effects, noise, and clustering
function synthetic_image = apply_advanced_features(base_image, current_conc, image_size)
    regional_variation = 1 + 0.3*imfilter(rand(image_size), fspecial('gaussian', 60, 10));
    synthetic_image = base_image .* regional_variation;
    non_specific = 0.05 * current_conc * rand(image_size);
    synthetic_image = synthetic_image + non_specific;
    cluster_prob = 0.1 * current_conc;
    cluster_mask = rand(image_size) < cluster_prob;
    synthetic_image(cluster_mask) = synthetic_image(cluster_mask) + 5*exprnd(1, sum(cluster_mask(:)), 1);
    psf = fspecial('gaussian', 15, 2);
    synthetic_image = imfilter(synthetic_image, psf, 'replicate');
    synthetic_image = synthetic_image + 0.1*randn(image_size);
end

%% Visualization helper
function visualize_results(image, conc, sim, binary_image, rand_sim, K_d,directory,m)
    if sim == rand_sim
        figure()
        imshow(image);
        filename = sprintf('ROI %d  Sim %.0f, Conc=(%.3f), K_d=(%.3f).fig',m, sim, (conc*1.0e6), (K_d*1.0e6));
        savefig(fullfile(directory, filename));
        close;
        
        figure('Position', [100 100 1200 400]);
        subplot(1,3,1);
        imagesc(image);
        title(sprintf('ROI %d conc=(%.3f), K_d=(%.3f)',m,(conc*1.0e6), (K_d*1.0e6))); colorbar; axis image;
        subplot(1,3,2);
        imshow(image,[]);
        title('Thresholded Image');
        subplot(1,3,3);
        histogram(image(image>0), 50); title('Intensity Histogram'); xlabel('Intensity'); ylabel('Count');
        set(gcf, 'Color', 'w');
        filename = sprintf('ROI %d Sim %.0f, Conc=(%.3f), K_d=(%.3f).fig',m,sim, (conc*1.0e6), (K_d*1.0e6));
        savefig(fullfile(directory, filename));
        close;
    end
end

%% Analysis & plotting
function analyze_results(results, concentrations, Kd_values,directory,m)
    concentrations  = concentrations.*1.0e6;
    Kd_values=Kd_values.*1.0e6;
    num_Kd = length(Kd_values);
    figure('Position', [100 100 1600 1200]);
    hold on;
    % Mean Concentration with varying K_d
    subplot(3,2,1); hold on;
    for k = 1:num_Kd
        plot(Kd_values(k), mean(results.bound_number(1,:,k)), 'o', 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration');
    ylabel('Mean number of antibody bound');
    title('Binding Number Curve');
    legend('show'); grid on;
   
    % SNR
    subplot(3,2,2); hold on;
    for k = 1:num_Kd
        plot(concentrations, results.snr(:,1,k), 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration');
    ylabel('SNR');
    title('Signal-to-Noise Ratio');
    legend('show'); grid on;
    
    % Coefficient of Variation (CV)
    subplot(3,2,3); hold on;
    for k = 1:num_Kd
        plot(concentrations, results.cv(:,1,k), 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration (ug/mL)');
    ylabel('CV');
    title('Coefficient of Variation');
    legend('show'); grid on;
   
    % Mean Cluster Size
    subplot(3,2,4); hold on;
    for k = 1:num_Kd
        mean_cluster_size = squeeze(mean(results.cluster_sizes(:,:,k), 2));
        plot(concentrations, mean_cluster_size, 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration (ug/mL)');
    ylabel('Mean Cluster Area (in pixels)');
    title('Cluster Area');
    legend('show'); grid on;
    
    % Cluster Count
    subplot(3,2,5); hold on;
    for k = 1:num_Kd
        mean_cluster_count = squeeze(mean(results.cluster_counts(:,:,k), 2));
        plot(concentrations, mean_cluster_count, 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration (ug/mL)');
    ylabel('Cluster Count');
    title('Number of Clusters');
    legend('show'); grid on;

    subplot(3,2,6); hold on;
    for k = 1:num_Kd
        mean_cluster_radii = squeeze(mean(results.cluster_radii(:,:,k), 2));
        plot(concentrations, mean_cluster_radii, 'LineWidth', 1, 'DisplayName', sprintf('K_d = %f', Kd_values(k)));
    end
    xlabel('Antibody Concentration (ug/mL)');
    ylabel('Cluster Radius');
    title('Cluster Radius');
    legend('show'); grid on;
    
    set(gcf, 'Color', 'w');
    filename = sprintf('Simulation_Results_All_Kd,ROI=%d.fig',m);
    savefig(fullfile(directory, filename));
    close;
end

function chosenNumber = chooseNumberSimple(bindingProb)
    if nargin < 1 || ~isscalar(bindingProb) || bindingProb < 0 || bindingProb > 1
        error('bindingProb must be a scalar between 0 and 1.');
    end
    r = rand();
    chosenNumber = r * (1 - bindingProb) + 1 * bindingProb;
end

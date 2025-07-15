clc;

% 1. Initialize Parameters
num_simulations = 300;
antibody_molecules = linspace(10.0e-8, 10.0e-7, 10);

num_receptors = 350;
image_size_pixels = [400, 400]; % Image dimensions in pixels 
K_d = linspace(1.0e-8, 5.0e-8, 3) + 0.2e-7 * abs(randn(1,3)); % Dissociation constant in M (moles/L)

% Physical parameters
pixel_size_m = 2.0e-8; % Pixel size in meters (m)

% Convert pixel size to micrometers for volume calculation clarity
pixel_size_um = pixel_size_m * 1e6; % 1 m = 1e6 um
% pixel_size_um will be 0.02 um

image_area_um2 = image_size_pixels(1) * image_size_pixels(2) * (pixel_size_um)^2; % Area in um^2
% image_area_um2 = 5000 * 3000 * (0.02)^2 = 15,000,000 * 0.0004 = 6000 um^2

cell_depth_um = 8; % Cell depth in micrometers (um)

volume_um3 = image_area_um2 * cell_depth_um * 1e-3; % Volume in um^3
% volume_um3 = 6000 um^2 * 6 um = 36000 um^3

volume_pL = volume_um3 * 1e-3; % Convert um^3 to pL (1 pL = 1000 um^3)
% volume_pL = 36000 * 1e-3 = 36 pL

volume_ml = volume_pL * 1e-9; % Convert pL to mL (1 mL = 1e9 pL)
% volume_ml = 36 * 1e-9 = 3.6e-8 mL

% Biological parameters
antibody_MW_kDa = 100; % Molecular Weight of antibody in kDa
% Mass of one antibody molecule in nanograms (1 Da = 1.660539e-24 g, 1 g = 1e9 ng)
mass_per_antibody_ug = antibody_MW_kDa * 1000 * 1.660539e-24 * 1e6;
% Simplified: mass_per_antibody_ng = antibody_MW_kDa * 1.660539e-12;
% mass_per_antibody_ng = 100 * 1.660539e-12 = 1.660539e-10 ng 
conc_antibody = (antibody_molecules .* mass_per_antibody_ug) ./ volume_ml;
directory = 'G:\Sim Data\rand_transfec';

results = antibody_staining_single_cell(num_simulations, conc_antibody, num_receptors, image_size_pixels, K_d,directory);

function [results] = antibody_staining_single_cell(num_simulations, antibody_concs, num_receptors, image_size, K_d,directory)
    
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
                receptor_map = generate_single_cell_receptors(image_size, num_receptors, 2);

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
                    visualize_results(synthetic_image, current_conc, sim, binary_image, rand_sim, current_K_d,directory);
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
    analyze_results(results, antibody_concs, K_d,directory);
end

%% Generate random receptor locations within a cell
function receptor_map = generate_single_cell_receptors(image_size, num_receptors, receptor_size_pixels)
    receptor_map = zeros(image_size);
    for i = 1:num_receptors
        center_x = round(rand() * image_size(2));
        center_y = round(rand() * image_size(1));
        half_size = floor(receptor_size_pixels / 2);
        x_coords = max(1, center_x - half_size):min(image_size(2), center_x + half_size);
        y_coords = max(1, center_y - half_size):min(image_size(1), center_y + half_size);
        for y = y_coords
            for x = x_coords
                receptor_map(y, x) = receptor_map(y, x) + 1;
            end
        end
    end
end

%% Simulate antibody staining
function [base_image, N_bound] = simulate_staining(receptor_map, current_conc, image_size, K_d)
    binding_prob = rand();
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
function visualize_results(image, conc, sim, binary_image, rand_sim, K_d,directory)
    if sim == rand_sim
        figure()
        imshow(image);
        filename = sprintf('OrgImg Sim %.0f, Conc=(%.3f), K_d=(%.3f).fig', sim, (conc*1.0e6), (K_d*1.0e6));
        savefig(fullfile(directory, filename));
        close;
        
        figure('Position', [100 100 1200 400]);
        subplot(1,3,1);
        imagesc(image);
        title(sprintf('Image conc=(%.3f), K_d=(%.3f)', (conc*1.0e6), (K_d*1.0e6))); colorbar; axis image;
        subplot(1,3,2);
        imshow(binary_image);
        title('Thresholded Image');
        subplot(1,3,3);
        histogram(image(image>0), 50); title('Intensity Histogram'); xlabel('Intensity'); ylabel('Count');
        set(gcf, 'Color', 'w');
        filename = sprintf('Sim %.0f, Conc=(%.3f), K_d=(%.3f).fig', sim, (conc*1.0e6), (K_d*1.0e6));
        savefig(fullfile(directory, filename));
        close;
    end
end

%% Analysis & plotting
function analyze_results(results, concentrations, Kd_values,directory)
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
    filename = 'Simulation_Results_All_Kd.fig';
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
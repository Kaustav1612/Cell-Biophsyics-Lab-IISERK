clc;
clearvars -except  filter_mobile
%%
[sRg_x,sRg_y,R_g_x,R_g_y,r_x,r_y] = rg_trajectory(all_trajectory);

filter_mobile = prctile(sRg_x,95);
clearvars -except  filter_mobile
%%

[sRg_x,sRg_y,R_g_x,R_g_y,r_x,r_y] = rg_trajectory(all_trajectory);
idx_mobile = find(sRg_x>filter_mobile);
idx_fixed = find(sRg_x<filter_mobile);

% 1. Extract the mobile trajectories
all_trajectory_mobile = all_trajectory(idx_mobile);

% 2. Extract the 9fixed trajectories
all_trajectory_fixed = all_trajectory(idx_fixed);

% Optional: Print a summary to the console
fprintf('Segregation Complete:\n');
fprintf(' - Mobile Clusters: %d\n', length(all_trajectory_mobile));
fprintf(' - Fixed Clusters:  %d\n', length(all_trajectory_fixed));

%%
clearvars -except all_trajectory filter_mobile all_trajectory_mobile all_trajectory_fixed idx_mobile idx_fixed
% Parameters
nframes = 300;
del_t = 0.2;
a = (1:nframes)' * del_t;  % Column vector of time points
old_dir= pwd;
directory='F:\Uday_data\mobility_analysis\Cyto_live_img_1\Control';
% Initialize arrays
TAMSD_mobile = NaN(nframes, size(all_trajectory_mobile,2));
TAMSD_fixed = NaN(nframes, size(all_trajectory_fixed,2));
slope_vector_mobile = [];
slope_vector_fixed = [];
intensity_fixed = NaN(nframes, size(all_trajectory_fixed,2));
intensity_mobile = NaN(nframes, size(all_trajectory_mobile,2));
area_fixed = NaN(nframes, size(all_trajectory_fixed,2));
area_mobile = NaN(nframes, size(all_trajectory_mobile,2));

% Calculate number of trajectories
num_traj_fixed = size(all_trajectory_fixed,2);
num_traj_mobile = size(all_trajectory_mobile,2);

% Initialize figures for plotting
figure('Position', [100, 100, 1200, 400]);

% Subplot for mobile particles
subplot(1,2,1);
hold on;
title('Mobile Particles TAMSD');
xlabel('Time Lag (s)');
ylabel('TAMSD (\mum^2)');
set(gca, 'XScale', 'log', 'YScale', 'log');

pixel = 0.065;
% Loop through each mobile particle
for j = 1:size(all_trajectory_mobile,2)
    X_Centroids_val = all_trajectory_mobile(j).positions(:,1)*0.065;
    Y_Centroids_val = all_trajectory_mobile(j).positions(:,2)*0.065;
    trajectory_length = all_trajectory_mobile(j).length;
    
    % Store intensity if field exists
    if isfield(all_trajectory_mobile(j), 'intensity')
        intensity_time = all_trajectory_mobile(j).intensity;
        intensity_mobile(1:length(intensity_time), j) = intensity_time';
    end
     if isfield(all_trajectory_mobile(j), 'area')
        area_time = all_trajectory_mobile(j).area;
        area_mobile(1:length(area_time), j) = area_time';
    end
    
    % Calculate TAMSD
    for n = 1:trajectory_length-1  % n is the time lag
        sum_MSD = 0;
        count = 0;
        
        % Average over all possible time origins
        for t = 1:(trajectory_length-n)
            dx = X_Centroids_val(t+n) - X_Centroids_val(t);
            dy = Y_Centroids_val(t+n) - Y_Centroids_val(t);
            sum_MSD = sum_MSD + dx^2 + dy^2;
            count = count + 1;
        end
        
        if count > 0 && sum_MSD ~= 0
            % Correct TAMSD calculation
            TAMSD_mobile(n,j) = sum_MSD / (count * del_t * n);  % n is the time lag in steps
        end
    end

    % Linear fit in log-log space
    valid_idx = find(TAMSD_mobile(:,j) > 0 & ~isnan(TAMSD_mobile(:,j)));
    if length(valid_idx) >= 5  % Increased minimum points for better fit
        valid_time = a(valid_idx);
        valid_TAMSD = TAMSD_mobile(valid_idx,j);
        
        % Fit only first 50% of data for better slope estimation
        max_fit_idx = min(valid_idx(round(length(valid_idx)*0.5)), valid_idx(end));
        fit_range = valid_idx(valid_idx <= max_fit_idx);
        
        if length(fit_range) >= 5
            fit_time = a(fit_range);
            fit_TAMSD = TAMSD_mobile(fit_range,j);
            
            p = polyfit(log(fit_time), log(fit_TAMSD), 1);
            slope = p(1);
            slope_vector_mobile = [slope_vector_mobile; slope];
            
            % Plot individual TAMSD curve
            plot(valid_time, valid_TAMSD, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
            
            % Plot fitted line
            log_fit = polyval(p, log(valid_time));
            plot(valid_time, exp(log_fit), 'r--', 'LineWidth', 1);
        end
    end
end

% Plot average TAMSD for mobile particles
avg_TAMSD_mobile = nanmean(TAMSD_mobile, 2);
valid_avg = find(~isnan(avg_TAMSD_mobile));
if ~isempty(valid_avg)
    plot(a(valid_avg), avg_TAMSD_mobile(valid_avg), 'b-', 'LineWidth', 2);
end

% Subplot for fixed particles
subplot(1,2,2);
hold on;
title('Fixed Particles TAMSD');
xlabel('Time Lag (s)');
ylabel('TAMSD (\mum^2)');
set(gca, 'XScale', 'log', 'YScale', 'log');


% Loop through each fixed particle
for j = 1:size(all_trajectory_fixed,2)
    X_Centroids_val = all_trajectory_fixed(j).positions(:,1)*pixel;
    Y_Centroids_val = all_trajectory_fixed(j).positions(:,2)*pixel;
    trajectory_length = all_trajectory_fixed(j).length;
    
    % Store intensity if field exists
    if isfield(all_trajectory_fixed(j), 'intensity')
        intensity_time = all_trajectory_fixed(j).intensity;
        intensity_fixed(1:length(intensity_time), j) = intensity_time';
    end
    % Store intensity if field exists
    if isfield(all_trajectory_fixed(j), 'area')
        area_time = all_trajectory_fixed(j).intensity;
        area_fixed(1:length(area_time), j) = area_time';
    end
    
    % Calculate TAMSD
    for n = 1:trajectory_length-1
        sum_MSD = 0;
        count = 0;
        
        for t = 1:(trajectory_length-n)
            dx = X_Centroids_val(t+n) - X_Centroids_val(t);
            dy = Y_Centroids_val(t+n) - Y_Centroids_val(t);
            sum_MSD = sum_MSD + dx^2 + dy^2;
            count = count + 1;
        end
        
        if count > 0 && sum_MSD ~= 0
            TAMSD_fixed(n,j) = sum_MSD / (count * del_t * n);  % n is the time lag in stepst;
        end
    end
    
    % Linear fit for fixed particles
    valid_idx = find(TAMSD_fixed(:,j) > 0 & ~isnan(TAMSD_fixed(:,j)));
    if length(valid_idx) >= 5
        valid_time = a(valid_idx);
        valid_TAMSD = TAMSD_fixed(valid_idx,j);
        
        p = polyfit(log(valid_time), log(valid_TAMSD), 1);
        slope = p(1);
        slope_vector_fixed = [slope_vector_fixed; slope];
        
        % Plot individual TAMSD
        plot(valid_time, valid_TAMSD, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    

end

% Plot average TAMSD for fixed particles
avg_TAMSD_fixed = nanmean(TAMSD_fixed, 2);
valid_avg = find(~isnan(avg_TAMSD_fixed));
if ~isempty(valid_avg)
    plot(a(valid_avg), avg_TAMSD_fixed(valid_avg), 'r-', 'LineWidth', 2);
end


%%

% --- User Parameters ---
window_size = 20;   % Number of time lag points in the rolling window
del_t = 0.2;        % Time step in seconds

        
   % --- Step 1: Characterize Photobleaching from Fixed Particles ---
% --- Calculate TAMSD_z from intensity first ---
[bleach_model, bleach_params] = characterize_photobleaching_from_fixed(all_trajectory_fixed, del_t);
    
[TAMSD_z, t_lags, intensity_corrected] = calculate_TAMSD_z(...
        all_trajectory_mobile, bleach_model, nframes, del_t);


% --- Setup ---
num_particles = size(TAMSD_mobile, 2);
rolling_results = cell(num_particles, 1);

% Initialize summary arrays
all_alpha_xy = [];
all_alpha_z = [];
all_D_xy = [];
all_D_z = [];
time_all = [];

fprintf('Processing %d particles with rolling window analysis...\n', num_particles);

% --- Main Processing Loop ---
for j = 1:num_particles
    % 1. Data Cleaning
    y_data = TAMSD_mobile(:, j);  % XY TAMSD
    z_data = TAMSD_z(:, j);       % Z TAMSD from intensity
    
    % Find valid indices (positive and non-NaN)
    valid_idx_y = find(y_data > 0 & ~isnan(y_data));
    valid_idx_z = find(z_data > 0 & ~isnan(z_data));
    
    % Skip if insufficient data
    if length(valid_idx_y) < window_size
        fprintf('Particle %d skipped: insufficient XY data (%d valid points)\n', j, length(valid_idx_y));
        continue;
    end
    
    % Extract valid XY data
    t_valid_y = t_lags(valid_idx_y);
    y_valid = y_data(valid_idx_y);
    
    % Extract valid Z data if available
    has_z_data = ~isempty(valid_idx_z) && length(valid_idx_z) >= window_size;
    if has_z_data
        t_valid_z = t_lags(valid_idx_z);
        z_valid = z_data(valid_idx_z);
    end
    
    % --- 2. Rolling Window for XY TAMSD ---
    num_windows_y = length(t_valid_y) - window_size + 1;
    alpha_xy = zeros(num_windows_y, 1);
    D_xy = zeros(num_windows_y, 1);
    t_center_y = zeros(num_windows_y, 1);
    
    for k = 1:num_windows_y
        win_start = k;
        win_end = k + window_size - 1;
        
        t_win = t_valid_y(win_start:win_end);
        y_win = y_valid(win_start:win_end);
        
        % Linear fit: log(TAMSD) = alpha*log(t) + log(4D)
        p = polyfit(log(t_win), log(y_win), 1);
        
        alpha_xy(k) = p(1);
        D_xy(k) = exp(p(2)) / 4;  % Extract D from intercept
        t_center_y(k) = mean(t_win);
    end
    
    % --- 3. Rolling Window for Z TAMSD ---
    if has_z_data
        num_windows_z = length(t_valid_z) - window_size + 1;
        alpha_z = zeros(num_windows_z, 1);
        D_z = zeros(num_windows_z, 1);
        t_center_z = zeros(num_windows_z, 1);
        
        for k = 1:num_windows_z
            win_start = k;
            win_end = k + window_size - 1;
            
            t_win = t_valid_z(win_start:win_end);
            z_win = z_valid(win_start:win_end);
            
            % Linear fit: log(TAMSD_z) = alpha*log(t) + log(2D)
            p_z = polyfit(log(t_win), log(z_win), 1);
            
            alpha_z(k) = p_z(1);
            D_z(k) = exp(p_z(2)) / 2;  % For 1D diffusion in z
            t_center_z(k) = mean(t_win);
        end
    else
        alpha_z = [];
        D_z = [];
        t_center_z = [];
    end
    
    % --- 4. Store Results ---
    rolling_results{j}.time_xy = t_center_y;
    rolling_results{j}.alpha_xy = alpha_xy;
    rolling_results{j}.D_xy = D_xy;
    rolling_results{j}.time_z = t_center_z;
    rolling_results{j}.alpha_z = alpha_z;
    rolling_results{j}.D_z = D_z;
    
    % Accumulate data for summary plots
    all_alpha_xy = [all_alpha_xy; alpha_xy];
    all_D_xy = [all_D_xy; D_xy];
    if has_z_data
        all_alpha_z = [all_alpha_z; alpha_z];
        all_D_z = [all_D_z; D_z];
    end
    time_all = [time_all; t_center_y];
end


% Alpha histogram
subplot(2, 3, 1);
histogram(all_alpha_xy, 30, 'FaceColor', 'b', 'FaceAlpha', 0.7, 'DisplayName', '\alpha_{xy}');
hold on;
if ~isempty(all_alpha_z)
    histogram(all_alpha_z, 30, 'FaceColor', 'r', 'FaceAlpha', 0.7, 'DisplayName', '\alpha_z');
end
xline(1, 'k--', 'LineWidth', 2);
xlabel('\alpha');
ylabel('Count');
title('Distribution of Diffusion Exponents');
legend;

% D histogram (log scale)
subplot(2, 3, 2);
histogram(log10(all_D_xy(all_D_xy > 0)), 30, 'FaceColor', 'b', 'FaceAlpha', 0.7, 'DisplayName', 'D_{xy}');
hold on;
if ~isempty(all_D_z)
    histogram(log10(all_D_z(all_D_z > 0)), 30, 'FaceColor', 'r', 'FaceAlpha', 0.7, 'DisplayName', 'D_z');
end
xlabel('log_{10}(D)');
ylabel('Count');
title('Distribution of Diffusion Coefficients');
legend;

% Alpha_xy vs Alpha_z scatter
subplot(2, 3, 3);
if ~isempty(all_alpha_z)
    % Match lengths for scatter plot
    min_len = min(length(all_alpha_xy), length(all_alpha_z));
    scatter(all_alpha_xy(1:min_len), all_alpha_z(1:min_len), 10, 'filled', ...
            'MarkerFaceAlpha', 0.3);
    hold on;
    plot([-1, 3], [-1, 3], 'k--');  
    xlabel('\alpha_{xy}');
    ylabel('\alpha_z');
    title('\alpha_{xy} vs \alpha_z');
    axis equal;
    grid on;
end


% Bin the data by time
time_bins = linspace(min(time_all), max(time_all), 20);
alpha_xy_mean = zeros(length(time_bins)-1, 1);
alpha_z_mean = zeros(length(time_bins)-1, 1);
alpha_xy_std = zeros(length(time_bins)-1, 1);
alpha_z_std = zeros(length(time_bins)-1, 1);

for i = 1:length(time_bins)-1
    bin_mask = time_all >= time_bins(i) & time_all < time_bins(i+1);
    alpha_xy_mean(i) = mean(all_alpha_xy(bin_mask));
    alpha_xy_std(i) = std(all_alpha_xy(bin_mask));
    if ~isempty(all_alpha_z)
        alpha_z_mean(i) = mean(all_alpha_z(bin_mask));
        alpha_z_std(i) = std(all_alpha_z(bin_mask));
    end
end
subplot(2, 3, 4:6);
time_centers = (time_bins(1:end-1) + time_bins(2:end)) / 2;
errorbar(time_centers, alpha_xy_mean, (alpha_xy_std/sqrt(size(all_trajectory_mobile,2))), 'b.-', 'LineWidth', 2, ...
         'MarkerSize', 10, 'DisplayName', '\alpha_{xy} (mean ± std)');
hold on;
if ~isempty(all_alpha_z)
    errorbar(time_centers, alpha_z_mean, (alpha_z_std/sqrt(size(all_trajectory_mobile,2))), 'r.-', 'LineWidth', 2, ...
             'MarkerSize', 10, 'DisplayName', '\alpha_z (mean ± std)');
end
yline(1, 'k--', 'LineWidth', 1);
xlabel('Time Lag \tau (s)');
ylabel('\alpha');
title('Time Evolution of Mean Diffusion Exponent');
legend;
grid on;

% --- Print Summary Statistics ---
fprintf('\n========== ROLLING WINDOW ANALYSIS SUMMARY ==========\n');
fprintf('Window size: %d points\n', window_size);
fprintf('Time step: %.2f s\n', del_t);
fprintf('Total particles analyzed: %d\n', sum(~cellfun(@isempty, rolling_results)));

fprintf('\n--- XY Diffusion ---\n');
fprintf('Mean α_xy = %.3f ± %.3f\n', mean(all_alpha_xy), (std(all_alpha_xy)/sqrt(size(all_trajectory_mobile,2))));
fprintf('Median α_xy = %.3f\n', median(all_alpha_xy));
fprintf('Mean D_xy = %.2e ± %.2e µm²/s\n', mean(all_D_xy), (std(all_D_xy)/sqrt(size(all_trajectory_mobile,2))));
fprintf('Median D_xy = %.2e µm²/s\n', median(all_D_xy));

if ~isempty(all_alpha_z)
    fprintf('\n--- Z Diffusion (from intensity) ---\n');
    fprintf('Mean α_z = %.3f ± %.3f\n', mean(all_alpha_z), (std(all_alpha_z)/sqrt(size(all_trajectory_mobile,2))));
    fprintf('Median α_z = %.3f\n', median(all_alpha_z));
    fprintf('Mean D_z = %.2e ± %.2e µm²/s\n', mean(all_D_z), (std(all_D_z)/sqrt(size(all_trajectory_mobile,2))));
    fprintf('Median D_z = %.2e µm²/s\n', median(all_D_z));
    
    % Test for anisotropy
    [h, p] = ttest2(all_alpha_xy, all_alpha_z);
    fprintf('\n--- Anisotropy Test ---\n');
    fprintf('α_xy vs α_z t-test: p = %.4f\n', p);
    if p < 0.05
        fprintf('Significant difference: Motion is ANISOTROPIC (p < 0.05)\n');
    else
        fprintf('No significant difference: Motion appears ISOTROPIC\n');
    end
end

% --- Optional: Export results ---
cd(directory);
save('rolling_window_results.mat', 'rolling_results', 'all_alpha_xy', ...
     'all_alpha_z', 'all_D_xy', 'all_D_z');
% Summary statistics
fprintf('\n=== Rolling Window Analysis Summary ===\n');
fprintf('Window size: %d points\n', window_size);
cd(old_dir);
valid_particles = find(~cellfun(@isempty, rolling_results));



% Display statistics
fprintf('\n=== TAMSD Analysis Results ===\n');
fprintf('Mobile Particles: %d trajectories\n', num_traj_mobile);
fprintf('  Mean slope: %.3f ± %.3f\n', mean(slope_vector_mobile), (std(slope_vector_mobile)/sqrt(size(all_trajectory_mobile,2))));
fprintf('  Slope range: [%.3f, %.3f]\n', min(slope_vector_mobile), max(slope_vector_mobile));
fprintf('\nFixed Particles: %d trajectories\n', num_traj_fixed);
fprintf('  Mean slope: %.3f ± %.3f\n', mean(slope_vector_fixed), (std(slope_vector_fixed)/sqrt(size(all_trajectory_fixed,2))));
fprintf('  Slope range: [%.3f, %.3f]\n', min(slope_vector_fixed), max(slope_vector_fixed));

% Create histogram of slopes
figure('Position', [100, 600, 800, 300]);
subplot(1,2,1);
histogram(slope_vector_mobile, 20);
xlabel('Slope (α)');
ylabel('Count');
title('Mobile Particles: Slope Distribution');
xline(1, 'r--', 'Normal Diffusion');

subplot(1,2,2);
histogram(slope_vector_fixed, 20);
xlabel('Slope (α)');
ylabel('Count');
title('Fixed Particles: Slope Distribution');
xline(0, 'r--', 'Immobile');
%%


% Calculate TAMSD_z directly from intensity fluctuations
function [TAMSD_z, t_lags, intensity_corrected] = calculate_TAMSD_z(all_trajectory_mobile, bleach_params, nframes, del_t)
    % Calculate TAMSD_z from intensity with photobleaching correction
    % Uses pre-characterized bleaching parameters from fixed particles
    %
    % Inputs:
    %   all_trajectory_mobile - structure array of mobile particles
    %   bleach_params - structure from characterize_photobleaching_from_fixed()
    %                  containing: model, params, fit_fun
    %   nframes - number of frames
    %   del_t - time step between frames
    %
    % Outputs:
    %   TAMSD_z - corrected axial TAMSD (intensity-based)
    %   t_lags - time lag vector
    %   intensity_corrected - cell array of photobleach-corrected intensities
    
    num_mobile = size(all_trajectory_mobile, 2);
    TAMSD_z = NaN(nframes, num_mobile);
    t_lags = (1:nframes)' * del_t;
    intensity_corrected = cell(num_mobile, 1);
    
    % Extract the bleaching correction function from fixed particle analysis
    if isfield(bleach_params, 'fit_fun')
        bleach_fun = bleach_params.fit_fun;
    else
        error('bleach_params must contain fit_fun field from characterize_photobleaching_from_fixed()');
    end
    
    fprintf('Processing %d mobile particles with photobleaching correction...\n', num_mobile);
    fprintf('Using bleaching model: %s\n', bleach_params.model);
    
    for j = 1:num_mobile
        % Extract intensity time series
        if ~isfield(all_trajectory_mobile(j), 'intensity')
            fprintf('Particle %d: No intensity data, skipping\n', j);
            continue;
        end
        
        I_raw = all_trajectory_mobile(j).intensity(:);
        trajectory_length = length(I_raw);
        
        if trajectory_length < 10
            fprintf('Particle %d: Trajectory too short (%d frames), skipping\n', j, trajectory_length);
            continue;
        end
        
        % --- Apply photobleaching correction ---
        t = (0:trajectory_length-1)' * del_t;
        
        % Calculate bleaching trend for this particle's time points
        I_bleach = bleach_fun(t);
        
        % Normalize bleaching curve to start at 1
        if I_bleach(1) > 0
            I_bleach_norm = I_bleach / I_bleach(1);
        else
            I_bleach_norm = ones(size(I_bleach));
        end
        
        % Method 1: Multiplicative correction (divide by bleaching trend)
        I_corrected_raw = I_raw ./ I_bleach_norm;
        
        % Rescale to preserve original intensity scale
        % Use first few points as reference (assuming particle starts near focus)
        n_ref = min(10, trajectory_length);
        I_ref_original = mean(I_raw(1:n_ref));
        I_ref_corrected = mean(I_corrected_raw(1:n_ref));
        
        if I_ref_corrected > 0
            I_corrected_raw = I_corrected_raw * (I_ref_original / I_ref_corrected);
        end
        
        % Method 2 (Alternative): Additive correction for comparison
        % I_trend = I_raw(1) * I_bleach_norm;
        % I_corrected_additive = I_raw - I_trend + mean(I_raw(1:min(20, trajectory_length)));
        
        % Optional: Light smoothing to reduce noise
        if trajectory_length > 20
            I_corrected_raw = smoothdata(I_corrected_raw, 'gaussian', 3);
        end
        
        % Ensure no negative values
        I_corrected_raw = max(I_corrected_raw, 1e-10);
        
        % Store corrected intensity
        intensity_corrected{j} = I_corrected_raw;
        
        % --- Calculate TAMSD_z from corrected intensity fluctuations ---
        % This measures axial mobility through intensity changes
        for n = 1:trajectory_length-1  % time lag
            sum_MSD_z = 0;
            count = 0;
            
            for tau = 1:(trajectory_length-n)
                % Direct intensity difference (measures total axial movement)
                dI = I_corrected_raw(tau+n) - I_corrected_raw(tau);
                
                % Alternative: Normalized intensity difference
                % dI_norm = dI / mean([I_corrected_raw(tau+n), I_corrected_raw(tau)]);
                
                sum_MSD_z = sum_MSD_z + dI^2;
                count = count + 1;
            end
            
            if count > 0
                TAMSD_z(n, j) = sum_MSD_z / count;
            end
        end
        
        % --- Quality check: Verify correction worked ---
        % Calculate trend in corrected vs raw intensity
        if trajectory_length > 30
            % Linear fit to check for remaining trend
            p_raw = polyfit(t, I_raw, 1);
            p_corr = polyfit(t, I_corrected_raw, 1);
            
            trend_reduction = (1 - abs(p_corr(1))/abs(p_raw(1))) * 100;
            
            if trend_reduction < 50
                fprintf('Particle %d: Warning - photobleaching correction may be incomplete (%.1f%% reduction)\n', ...
                        j, trend_reduction);
            end
        end
        
        % Progress update
        if mod(j, 20) == 0
            fprintf('Processed %d/%d particles\n', j, num_mobile);
        end
    end
    
    % --- Summary statistics ---
    n_processed = sum(~cellfun(@isempty, intensity_corrected));
    fprintf('\n=== Photobleaching Correction Summary ===\n');
    fprintf('Total mobile particles: %d\n', num_mobile);
    fprintf('Successfully processed: %d\n', n_processed);
    fprintf('Bleaching model used: %s\n', bleach_params.model);
    
    % Calculate average bleaching correction factor
    t_end = nframes * del_t;
    correction_at_end = bleach_fun(t_end) / bleach_fun(0);
    fprintf('Bleaching correction at %.1f s: %.1f%% of initial intensity\n', ...
            t_end, correction_at_end * 100);
    

    sgtitle('Photobleaching Correction Results');
end

function [bleach_model, bleach_params] = characterize_photobleaching_from_fixed(all_trajectory_fixed, del_t)
    
    num_fixed = size(all_trajectory_fixed, 2);
    all_intensities = [];
    all_times = [];
    
    fprintf('Characterizing photobleaching from %d fixed particles...\n', num_fixed);
    
    for j = 1:num_fixed
        if ~isfield(all_trajectory_fixed(j), 'intensity')
            continue;
        end
        
        I_raw = all_trajectory_fixed(j).intensity(:);
        trajectory_length = length(I_raw);
        
        % Only use particles with sufficient data
        if trajectory_length < 50
            continue;
        end
        
        % Normalize each particle's intensity to its initial value
        I_norm = I_raw / mean(I_raw(1:min(10, trajectory_length)));  % Normalize to first 10 frames
        
        % Time vector for this particle
        t = (0:trajectory_length-1)' * del_t;
        
        % Accumulate data for ensemble averaging
        all_intensities = [all_intensities; I_norm];
        all_times = [all_times; t];
    end
    
    % --- Fit ensemble photobleaching curve ---
    
    % Bin the data by time for robust fitting
    time_bins = 0:del_t:max(all_times);
    I_binned = zeros(size(time_bins));
    I_std = zeros(size(time_bins));
    
    for i = 1:length(time_bins)-1
        bin_mask = all_times >= time_bins(i) & all_times < time_bins(i+1);
        if sum(bin_mask) > 0
            I_binned(i) = mean(all_intensities(bin_mask));
            I_std(i) = std(all_intensities(bin_mask));
        end
    end
    
    % Remove empty bins
    valid_bins = I_binned > 0;
    time_bins = time_bins(valid_bins);
    I_binned = I_binned(valid_bins);
    I_std = I_std(valid_bins);
    
    % Try multiple photobleaching models and select the best
    models = {'single_exp', 'double_exp', 'exp_plus_linear'};
    best_model = '';
    best_r2 = -Inf;
    best_params = [];
    
    for m = 1:length(models)
        try
            switch models{m}
                case 'single_exp'
                    % I(t) = A*exp(-t/tau) + C
                    exp_fun = @(p, t) p(1) * exp(-t/p(2)) + p(3);
                    p0 = [0.8, 10, 0.2];  % A, tau, C
                    lb = [0, 0.1, 0];
                    ub = [1, 100, 1];
                    
                case 'double_exp'
                    % I(t) = A1*exp(-t/tau1) + A2*exp(-t/tau2) + C
                    exp_fun = @(p, t) p(1)*exp(-t/p(2)) + p(3)*exp(-t/p(4)) + p(5);
                    p0 = [0.5, 2, 0.3, 20, 0.2];
                    lb = [0, 0.1, 0, 1, 0];
                    ub = [1, 10, 1, 100, 1];
                    
                case 'exp_plus_linear'
                    % I(t) = A*exp(-t/tau) + B*t + C
                    exp_fun = @(p, t) p(1)*exp(-t/p(2)) + p(3)*t + p(4);
                    p0 = [0.8, 10, -0.001, 0.2];
                    lb = [0, 0.1, -0.01, 0];
                    ub = [1, 100, 0, 1];
            end
            
            % Fit the model
            p_fit = lsqcurvefit(exp_fun, p0, time_bins, I_binned, lb, ub);
            
            % Calculate R²
            I_pred = exp_fun(p_fit, time_bins);
            SS_res = sum((I_binned - I_pred).^2);
            SS_tot = sum((I_binned - mean(I_binned)).^2);
            R2 = 1 - SS_res/SS_tot;
            
            fprintf('Model: %s, R² = %.4f\n', models{m}, R2);
            
            if R2 > best_r2
                best_r2 = R2;
                best_model = models{m};
                best_params = p_fit;
            end
            
        catch
            fprintf('Model %s fitting failed\n', models{m});
        end
    end
    
    fprintf('\nBest model: %s (R² = %.4f)\n', best_model, best_r2);
    
    % Store the model
    bleach_model.model = best_model;
    bleach_model.params = best_params;
    bleach_model.time_bins = time_bins;
    bleach_model.I_binned = I_binned;
    bleach_model.I_std = I_std;
    bleach_params = best_params;
    

    

    
    % Return the fitted function handle
    bleach_model.fit_fun = @(t) exp_fun(best_params, t);
end

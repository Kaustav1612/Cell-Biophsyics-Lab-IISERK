function [Gp, Gpp, Gabs,Gp_mean,Gpp_mean,omega_unique] = mason_weitz_MSD(TAMSD,sizes,T,num_trajectories)
kB   = 1.380649e-23;
% Preallocate arrays to store all data
all_omega = [];
all_Gp = [];
all_Gpp = [];
all_Gabs = [];
all_G_r = [];

% Arrays to store data for individual particle plotting (optional)
omega_cell = {};
Gp_cell = {};
Gpp_cell = {};
Gabs_cell = {};
G_r_cell = {};

for particle = 1:num_trajectories
    dt_frame = 0.2;
    frames = size(find(TAMSD(:,particle)>0));
    lag_frames = 1:frames;                 % compute MSD up to 100-lag
    t = lag_frames * dt_frame;          % in seconds
    t = t';
    a = TAMSD(:,particle);
    % Create a logical index mask
    valid_indices = TAMSD(:,particle) > 0;
    
    % Use the logical mask to extract the corresponding values from 'a'
    MSD = a(valid_indices)*1e-12; % conversion into SI units um^2/s to m^2/s 
    rad = sizes(particle);
    
    if length(MSD) > 2
        % Choose frequencies as 1./t (reciprocal of lag times)
        omega = 1./t;  % rad/s approximately
        omega = omega(:);

        % Compute local slope alpha = d log(MSD)/d log(t)
        log_t = log(t);
        log_MSD = log(MSD);
        dlogMSD_dlogt = gradient(log_MSD)./gradient(log_t);
        alpha = dlogMSD_dlogt;

        % Interpolate MSD and alpha at t=1/omega
        MSD_interp = interp1(t', MSD, 1./omega, 'linear', 'extrap');
        alpha_interp = interp1(t', alpha, 1./omega, 'linear', 'extrap');

        % Compute |G*|
        Gamma_val = gamma(1 + alpha_interp); % Gamma function
        Gabs_particle = (kB*T) ./ (pi*rad*MSD_interp .* Gamma_val);
        
        Gp_particle = Gabs_particle .* cos(pi*alpha_interp/2);
        Gpp_particle = Gabs_particle .* sin(pi*alpha_interp/2);

        Gw = Gp_particle + i* Gpp_particle;
        G_r_particle = ifft(real(Gw));
        
        % Filter out any non-positive or NaN values
        valid_idx = (omega > 0) & (Gp_particle > 0) & (Gpp_particle > 0) & ~isnan(Gp_particle) & ~isnan(Gpp_particle);
        
        % Store filtered data
        omega_filtered = omega(valid_idx);
        Gp_filtered = Gp_particle(valid_idx);
        Gpp_filtered = Gpp_particle(valid_idx);
        G_r_filtered = G_r_particle(valid_idx);
        
        % Append to collective arrays
        all_omega = [all_omega; omega_filtered];
        all_Gp = [all_Gp; Gp_filtered];
        all_Gpp = [all_Gpp; Gpp_filtered];
        all_G_r = [all_G_r; G_r_filtered];
        
        % Also store in cell arrays for individual tracking
        omega_cell{particle} = omega_filtered;
        Gp_cell{particle} = Gp_filtered;
        Gpp_cell{particle} = Gpp_filtered;
        G_r_cell{particle} = G_r_filtered;
        
        % Store in the original cell arrays (if needed elsewhere)
        Gp{particle} = Gp_particle;
        Gpp{particle} = Gpp_particle;
        Gabs{particle} = Gabs_particle;
    end
    size(all_omega);
end

% Now create the collective plots
if ~isempty(all_omega)

    
    % Optionally add mean lines
    [omega_unique, ~, idx] = unique(all_omega);
    Gp_mean = accumarray(idx, all_Gp, [], @mean);
    Gpp_mean = accumarray(idx, all_Gpp, [], @mean);
    
    % plot(omega_unique, Gp_mean, 'b-', 'LineWidth', 3);
    % plot(omega_unique, Gpp_mean, 'r-', 'LineWidth', 3);
    

    % Plot 2: G' and G'' on log-log scale
    figure;
    hold on;
    scatter(all_omega, all_Gp, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    scatter(all_omega, all_Gpp, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.3);
    plot(omega_unique, Gp_mean, 'b-', 'LineWidth', 3);
    plot(omega_unique, Gpp_mean, 'r-', 'LineWidth', 3);
    
    xlabel('\omega (rad/s)');
    ylabel('G'' and G'''' (Pa)');
    legend('G''(\omega)', 'G''''(\omega)', 'Mean G''', 'Mean G''''');
    title('Complex modulus from MSD (Mason-Weitz) - All Particles (Log-Log)');
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');

    % Plot 3: G(r) on linear scale
    figure;
    hold on;
    scatter(all_omega, all_G_r, 20, 'g', 'filled', 'MarkerFaceAlpha', 0.3);
    
    G_r_mean = accumarray(idx, all_G_r, [], @mean);
    plot(omega_unique, G_r_mean, 'g-', 'LineWidth', 3);
    
    xlabel('\omega (rad/s)');
    ylabel('G(r) (Pa)');
    legend('G(r)', 'Mean G(r)');
    title('Visco-elastic modulus from MSD (Mason-Weitz) - All Particles');
    grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');

end
end
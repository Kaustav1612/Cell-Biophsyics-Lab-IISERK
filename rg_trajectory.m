function [sRg_x, sRg_y, R_g_x, R_g_y, r_x, r_y] = rg_trajectory_vectorized(trajectories)
    % Vectorized version for better performance
    
    num_traj = length(trajectories);
    
    % Initialize outputs
    R_g_x = NaN(num_traj, 1);
    R_g_y = NaN(num_traj, 1);
    r_x = NaN(num_traj, 1);
    r_y = NaN(num_traj, 1);
    
    for i = 1:num_traj
        pos = trajectories(i).positions;
        N = size(pos, 1);
        
        if N < 2
            continue;
        end
        
        % Mean position
        mean_pos = mean(pos, 1, 'omitnan');
        
        % Radius of gyration
        R_g_x(i) = sqrt(mean((pos(:,1) - mean_pos(1)).^2, 'omitnan'));
        R_g_y(i) = sqrt(mean((pos(:,2) - mean_pos(2)).^2, 'omitnan'));
        
        % Mean step size
        steps = abs(diff(pos, 1, 1));
        r_x(i) = mean(steps(:,1), 'omitnan');
        r_y(i) = mean(steps(:,2), 'omitnan');
    end
    
    % Scaled radius of gyration
    sRg_x = sqrt(pi/2) * (R_g_x ./ r_x);
    sRg_y = sqrt(pi/2) * (R_g_y ./ r_y);
end
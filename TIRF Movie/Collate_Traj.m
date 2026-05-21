% 1. Initialize structure with the TAMSD field
all_trajectory = struct(...
    'trajectory_id', {}, ...
    'particle_id', {}, ...
    'frames', {}, ...
    'positions', {}, ...
    'intensity', {}, ...
    'area', {}, ...
    'length', {}, ...
    'max_gap', {}, ...
    'valid_frames', {});           

%% 2. Loop through trajectories
for i = 1:length(trajectories)
    
    % Get the next available index in the master structure
    idx = length(all_trajectory) + 1;
    
    all_trajectory(idx).trajectory_id = idx; 
    all_trajectory(idx).particle_id   = trajectories(i).particle_id;
    all_trajectory(idx).frames        = trajectories(i).frames;
    all_trajectory(idx).positions     = trajectories(i).positions;
    all_trajectory(idx).intensity     = trajectories(i).intensity;
    all_trajectory(idx).area          = trajectories(i).area;
    all_trajectory(idx).max_gap       = trajectories(i).max_gap;
    all_trajectory(idx).valid_frames  = trajectories(i).valid_frames;
    
    % Calculate length
    all_trajectory(idx).length        = length(trajectories(i).frames);
    
end
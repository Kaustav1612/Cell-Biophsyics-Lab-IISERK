function [drift_corrected_traj] = drift_removal(trajs,nframes)
% trajs: cell array where trajs{i} is Nx2 (NaN for missing)
% compute reference positions (median per particle)
num_particles = size(trajs,2);
% Initialize 'ref' to hold 2D median coordinates
ref = zeros(num_particles, 2); 
drift_corrected_traj  = trajs;
% Loop to calculate the median position for each particle
for k = 1 : num_particles
    % nanmedian returns a 1x2 array [median_x, median_y].
    % 'ref' is now correctly sized to accept this.
    ref(k, :) = nanmedian(trajs(k).positions);
end
shifts = zeros(nframes,2);
for t=1:nframes
    diffs = [];
    for id=1:num_particles
        if t<=size(trajs(id).frames,2)
        if ~isnan(trajs(id).positions(t,1))
            diffs(end+1,:) = trajs(id).positions(t,:) - ref(id,:);
        end
        else
            continue
        end
    end
    if ~isempty(diffs)
        shifts(t,:) = median(diffs,1);
    else
        shifts(t,:) = 0;
    end
end
for t = 1 :nframes    
for i = 1:num_particles
    if ismember(trajs(id).frames,t)
        c(id).positions = trajs(id).positions - shifts(t,:);
    else
        continue
    end
end
end

end
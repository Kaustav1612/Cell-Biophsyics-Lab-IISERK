function [Centroids_3D,frames,X_Centroids,Y_Centroids,particle_presence,max_particles] = follow_track(new_stats, nframes,max_particles,index_max)
    % Initialize arrays to collect centroids
   

    % Initialize arrays to store all X and Y centroids
    X_Centroids = NaN(max_particles, nframes);
    Y_Centroids = NaN(max_particles, nframes);
    Centroids_3D = NaN(max_particles, nframes, 3); 
    
    particle_presence = NaN(nframes, max_particles);

    % Initialize arrays to store all X and Y centroids
    %X_Centroids = NaN(initial_max_particles, nframes);
    %Y_Centroids = NaN(initial_max_particles, nframes);
    
    
    
%   counter=0;
    % Populate the centroid arrays
    for i = 1:nframes
        for j = 1:max_particles
            current_frame = new_stats{i,1};
                if ~isnan(new_stats{i,1}(j,:).Centroid(1))         
                      X_Centroids(j, i) =current_frame(j).Centroid(1);
                      Y_Centroids(j, i) = current_frame(j).Centroid(2);
                      particle_presence(i,j)=1;
                else
                      particle_presence(i,j)=0;
                end

        end
    end

    % Calculate the norm of the centroids
    for i = 1:nframes
        for j = 1:max_particles
            Centroids_3D(j, i, 1) = X_Centroids(j, i); % X centroid
            Centroids_3D(j, i, 2) = Y_Centroids(j, i); % Y centroid
            Centroids_3D(j, i, 3) = i;                % Frame number
        end
    end

    % Generate frame numbers array
    frames = 1:nframes;
end
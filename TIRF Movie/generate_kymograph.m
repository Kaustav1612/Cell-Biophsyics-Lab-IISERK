function [kymograph_position,kymograph_intensity,particle_kymograph,dist_from_centroid] = ...
    generate_kymograph(image_stack,binary_stack,x_coords,y_coords,all_stats,nframes,minimum_particle_size,max_particles,num_particles_detected,current_frame,X_Centroids,Y_Centroids)
% GENERATE KYMOGRAPH ALONG X
kymograph_position = zeros(nframes, (2*minimum_particle_size)+1, max_particles);
kymograph_intensity = zeros(nframes, (2*minimum_particle_size)+1, max_particles);
particle_kymograph=zeros(nframes,max_particles);
dist_from_centroid=zeros(nframes,max_particles);
x=round(x_coords);
y=round(y_coords);
x_round = round(X_Centroids);
y_round = round(Y_Centroids);

    for j = 1:max_particles
        for i = 1:nframes
            actual_n_particles = num_particles_detected(i);
            stats_frame = all_stats{i};
            binary_frame = binary_stack(:,:,i); 
            intensity_frame = double(image_stack(:,:, i));           
           if actual_n_particles>=j && i < length(X_Centroids)
               if ~isnan(y_round(j,i)) &&  ~isnan(x_round(j,i)) 
                    kymograph_position(i, :, j) = binary_frame(y(j), x(j)-minimum_particle_size:x(j)+minimum_particle_size);
                    kymograph_intensity(i, :, j) = intensity_frame(y(j), x(j)-minimum_particle_size:x(j)+minimum_particle_size);

                    valid_indices_kymo = find(([stats_frame.XC] > x(j)-minimum_particle_size) & ([stats_frame.XC] < x(j)+minimum_particle_size) & ([stats_frame.YC] > y(j)-minimum_particle_size) & ([stats_frame.YC] < y(j)+minimum_particle_size));

                    if numel(valid_indices_kymo) > 0
                        is_particle_detected = 1;
                    else
                        is_particle_detected = 0;
                    end
                    %particle_kymograph{j} = binary_stack(y_coords(j), x_coords(j)-minimum_particle_size:x_coords(j)+minimum_particle_size, :);
                 
                    
               else
                   kymograph_intensity(i, :, j)=mean(current_frame(:));
               end
            elseif actual_n_particles<j 
                    kymograph_intensity(i, :, j)=mean(current_frame(:));
                    continue;
            end
            %mean_intensity_kymograph(i,j) = mean(image_stack(y_coords(j), x_coords(j)-minimum_particle_size:x_coords(j)+minimum_particle_size, :));       
        end       
    end
end

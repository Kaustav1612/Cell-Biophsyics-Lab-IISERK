function [kymograph_stitched,kymograph_stitched_intensity,total_lifetime_normalized] = get_lifetime(num_particles_detected,kymograph_position,kymograph_intensity,particle_presence,nframes,max_particles)

kymograph_stitched = [];
kymograph_stitched_intensity = [];
total_lifetime_normalized = zeros(max(num_particles_detected),1);
max_lifetime_detected=0;
threshold_for_non_detection=4;
particle_presence_smooth=zeros(nframes,max_particles);
    for j = 1:max(num_particles_detected)
                kymograph_stitched = [kymograph_stitched, kymograph_position(:,:,j)];
                kymograph_stitched_intensity = [kymograph_stitched_intensity, kymograph_intensity(:,:,j)];
    end
   for j = 1 : max(num_particles_detected)
       counter=0;
       lifetime=0;
        for i =1 :nframes
            if particle_presence(i,j)==1
                lifetime=lifetime+1;
            elseif particle_presence(i,j)==0 && counter <= threshold_for_non_detection
                counter = counter+1;
            end
        end
        total_lifetime_normalized(j)=lifetime;
   end
end


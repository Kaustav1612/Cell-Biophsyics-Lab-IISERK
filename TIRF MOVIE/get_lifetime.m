function [kymograph_stitched,kymograph_stitched_intensity,total_lifetime_normalized] = get_lifetime(num_particles_detected,kymograph_position,kymograph_intensity,nframes,max_particles)

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
   
end


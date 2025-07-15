function [num_stuck_particles,num_mobile_particles,stuck_image_full,mobile_image_full,stuck_image,mobile_image] = find_stuck(stuck_threshold_low,stuck_threshold_high,mobile_threshold_low,mobile_threshold_high,labeled_movie,num_particles_detected,total_lifetime_normalized,max_index)

% Initialize images for stuck and mobile particles
stuck_image = zeros(size(labeled_movie(:,:,max_index)));
mobile_image = zeros(size(labeled_movie(:,:,max_index)));
num_stuck_particles=0;
num_mobile_particles=0;
for i = 1:max(num_particles_detected)  
    if total_lifetime_normalized(i) >= stuck_threshold_low && total_lifetime_normalized(i) <= stuck_threshold_high
        stuck_image(labeled_movie(:,:,1) == i) = 1;
        num_stuck_particles=num_stuck_particles+1;
    end
    if total_lifetime_normalized(i) >= mobile_threshold_low && total_lifetime_normalized(i) <= mobile_threshold_high
        mobile_image(labeled_movie(:,:,1) == i) = 1;
        num_mobile_particles=num_mobile_particles+1;
    end
end

stuck_image(stuck_image>0)=1;
mobile_image(mobile_image>0)=1;
% Store stuck image for each cell
num_cells = 1;
mobile_image_full(:,:,num_cells) = mobile_image;
stuck_image_full(:,:,num_cells) = stuck_image;

 
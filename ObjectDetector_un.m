function [new_image, gauss_blur,threshold] = ObjectDetector_un(image, correlation_threshold, line_scan_length, mask,sigma,new_image,threshold,Min_Area_Thresh,Max_Area_Thresh)
    % Step 1: Convert image to double
    image = double(image);

    % Step 2: Apply Gaussian blur
    gauss_blur = imgaussfilt(image, 1, 'FilterSize', [1 1]);
    image(~mask) = 0;
    % Step 3: Initialize variables
    valid_objects = zeros(0, 2);
    correlation_values = zeros(0, 2);

    % Step 4: Find pixels above the mean intensity
    if ~isempty(threshold)
        linear_indices = find(image>=mean(threshold));
    else 
        linear_indices = find(image>=mean(image));
    end
    [rows, cols] = ind2sub(size(image), linear_indices);

    % Step 5: Loop through each selected region
  
    for k = 1:length(rows)
        % Check if masks are valid
        roi_mask1 = mask(rows(k), max(1, cols(k) - line_scan_length) : min(size(mask, 2), cols(k) + line_scan_length));
        roi_mask2 = mask(max(1, rows(k) - line_scan_length): min(size(mask, 1), rows(k) + line_scan_length), cols(k));

        if all(roi_mask1(:)) && all(roi_mask2(:))
         
             % Get the line scans
            x_line_scan_main = gauss_blur(rows(k), max(1, cols(k) - line_scan_length) : min(size(image, 2), cols(k) + line_scan_length));
            y_line_scan_main = gauss_blur(max(1, rows(k) - line_scan_length) : min(size(image, 1), rows(k) + line_scan_length), cols(k));
            % Create 1D Gaussian PSFs
            x_vals = linspace(-5, 5, length(x_line_scan_main));
            y_vals = linspace(-5, 5, length(x_line_scan_main));
            gaussian_psf_x = normpdf(x_vals, 0, sigma);
            gaussian_psf_y = normpdf(y_vals, 0, sigma);

            % Normalize the line scans and Gaussian PSFs
            x_line_scan = (x_line_scan_main - mean(x_line_scan_main)) / std(x_line_scan_main);
            y_line_scan = (y_line_scan_main - mean(y_line_scan_main)) / std(y_line_scan_main);
            gaussian_psf_x = (gaussian_psf_x - mean(gaussian_psf_x)) / std(gaussian_psf_x);
            gaussian_psf_y = (gaussian_psf_y - mean(gaussian_psf_y)) / std(gaussian_psf_y);

            % Calculate correlations
            correlation_x = corr(x_line_scan(:), gaussian_psf_x(:));
            correlation_y = corr(y_line_scan(:), gaussian_psf_y(:));

            % Check if both correlations exceed the threshold
            if correlation_x > correlation_threshold && correlation_y > correlation_threshold
                % Create a circular mask around the current point
                threshold = [threshold,min(x_line_scan_main)];
                threshold = [threshold,min(y_line_scan_main)];
                [radius_x, radius_y] = calculate_fwhm_and_radius(x_line_scan, y_line_scan);
                elliptical_mask = create_elliptical_mask(image, rows(k), cols(k), radius_x, radius_y);

                % Update valid objects and new_image
                valid_idx = find(elliptical_mask);
                valid_idx = valid_idx(image(valid_idx(:))>=prctile(image(valid_idx(:)),30));
                if length(valid_idx) >= Min_Area_Thresh && length(valid_idx) < Max_Area_Thresh
                    [rows_obj, cols_obj] = ind2sub(size(image), valid_idx);
                    new_pixels = [rows_obj, cols_obj];
                    new_objects = setdiff(new_pixels, valid_objects, 'rows');  

                    if ~isempty(new_objects)
                        new_image(sub2ind(size(image), new_objects(:,1), new_objects(:,2))) = image(sub2ind(size(image), new_objects(:,1), new_objects(:,2)));
                    end
                elseif length(valid_idx) >= Min_Area_Thresh && length(valid_idx) > Max_Area_Thresh
                    valid_idx = valid_idx(image(valid_idx(:))>=prctile(image(valid_idx(:)),50));
                    [rows_obj, cols_obj] = ind2sub(size(image), valid_idx);
                    new_pixels = [rows_obj, cols_obj];
                    new_objects = setdiff(new_pixels, valid_objects, 'rows');  

                    if ~isempty(new_objects)
                        new_image(sub2ind(size(image), new_objects(:,1), new_objects(:,2))) = image(sub2ind(size(image), new_objects(:,1), new_objects(:,2)));
                    end
                end
            end
            
        end
    end
    
    end


function [fwhm_x, fwhm_y, radius_x, radius_y, anisotropy] = calculate_fwhm_and_radius(linescan_x, linescan_y)

% Calculate FWHM for each linescan
fwhm_x = calculate_fwhm(linescan_x);
fwhm_y = calculate_fwhm(linescan_y);

% Calculate radii
sigma_x = fwhm_x / (2 * sqrt(2 * log(2)));
sigma_y = fwhm_y / (2 * sqrt(2 * log(2)));
radius_x = sigma_x;
radius_y = sigma_y;

% Calculate anisotropy
anisotropy = radius_x / radius_y;

end

function [fwhm] = calculate_fwhm(linescan)

% Find the maximum value and its index
[max_value, max_index] = max(linescan);

% Find the half maximum value
half_max = max_value / 2;

% Find the indices of the points closest to half maximum on both sides
left_index = find(linescan(1:max_index) <= half_max, 1, 'last');
right_index = find(linescan(max_index:end) <= half_max, 1, 'first') + max_index - 1;

% Calculate the full width at half maximum
fwhm = right_index - left_index + 1;
end


function elliptical_mask = create_elliptical_mask(image, rows_k, cols_k, radius_x, radius_y)

    [rows_ext, cols_ext] = size(image);

    % Create a meshgrid of x and y coordinates
    [X, Y] = meshgrid(1:cols_ext, 1:rows_ext);

    % Calculate the elliptical mask using the equation of an ellipse
    elliptical_mask = (X - cols_k).^2 / radius_x^2 + (Y - rows_k).^2 / radius_y^2 <= 1;

end
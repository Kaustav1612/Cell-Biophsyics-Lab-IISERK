

    function [features,intensity_map] = feature_extraction_space(img,norm_image)
        % Calculate features
        [rows, cols] = find(img);
        intensity_map  = zeros (size(img));
        intensity_map(rows, cols) = norm_image(rows, cols);
         
                % The following represent the intensity and spatial features of the objects
                    % detected using the original image
                features = [...
                    rows,...
                    cols...
                ];
                
end
function [new_stats] = update_all_stats(all_stats,index_max,nframes,max_particles,max_movement_allowed_per_frame,total_max_movement)

% Create the known structure with all elements set to NaN
knownStructure = all_stats{index_max,1};  % Get the structure

% Loop through all fields of the structure and assign NaN
fields = fieldnames(knownStructure);  % Get all field names of the structure
for k = 1:numel(knownStructure)
    for f = 1:numel(fields)
        fieldValue = knownStructure(k).(fields{f});
        
        % Assign NaN based on the field type
        if isnumeric(fieldValue)
            knownStructure(k).(fields{f}) = NaN(size(fieldValue));  % Set to NaN
        elseif iscell(fieldValue)
            knownStructure(k).(fields{f}) = {NaN};  % Set cell field to NaN
        % Add other conditions if you have more complex fields
        end
    end
end

% Create the cell array 
new_stats = cell(nframes,1);

% Fill the cell array with copies of the known structure
for i = 1:nframes
    new_stats{i} = knownStructure;  % Copy the structure with NaNs
end

 new_stats{index_max,1}=all_stats{index_max,1};   


    for i = index_max-1:-1:1
        current_frame = all_stats{i,1};
        prev_frame = all_stats{i+1,1};
        for k = 1 : size(prev_frame,1)
            for j = 1 : size(current_frame,1)
                    if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        for l = 1 : max_particles
                           if  pdist2([all_stats{index_max,1}(l,:).Centroid],[current_frame(j,:).Centroid]) < total_max_movement
                          frame = all_stats{i,1}(j,:);
                          new_stats{1,i}(l,:).Area=frame.Area;
                          new_stats{1,i}(l,:).Centroid=frame.Centroid;
                          new_stats{1,i}(l,:).Eccentricity=frame.Eccentricity;
                          new_stats{1,i}(l,:).PixelIdxList=frame.PixelIdxList;
                          new_stats{1,i}(l,:).Intensity=frame.Intensity;
                          new_stats{1,i}(l,:).EquivDiameter=frame.EquivDiameter;
                          new_stats{1,i}(l,:).XC=frame.XC;
                          new_stats{1,i}(l,:).YC=frame.YC;
                           end
                        end
                 
                    end
                 
                    end
            end
        end
    end
    
    for i  = index_max+1 : nframes
        current_frame = all_stats{i,1};
        prev_frame = all_stats{i-1,1};
        for k = 1 : size(prev_frame,1)
            for j = 1 : size(current_frame,1)
                    if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        for l = 1 : max_particles
                           if  pdist2([all_stats{index_max,1}(l,:).Centroid],[current_frame(j,:).Centroid]) < total_max_movement
                          frame = all_stats{i,1}(j,:);
                          new_stats{i,1}(l,:).Area=frame.Area;
                          new_stats{i,1}(l,:).Centroid=frame.Centroid;
                          new_stats{i,1}(l,:).Eccentricity=frame.Eccentricity;
                          new_stats{i,1}(l,:).PixelIdxList=frame.PixelIdxList;
                          new_stats{i,1}(l,:).Intensity=frame.Intensity;
                          new_stats{i,1}(l,:).EquivDiameter=frame.EquivDiameter;
                          new_stats{i,1}(l,:).XC=frame.XC;
                          new_stats{i,1}(l,:).YC=frame.YC;
                           end
                        end
                 
                    end
            end
        end
    end
end

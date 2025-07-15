function [new_stats] = update_all_stats(all_stats,index_max,nframes,max_particles,max_movement_allowed_per_frame,total_max_movement)

% Create the known structure with all elements set to NaN
knownStructure = all_stats{1, index_max};  % Get the structure

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

% Create the cell array with 1x200 elements
new_stats = cell(1, nframes);

% Fill the cell array with copies of the known structure
for i = 1:nframes
    new_stats{i} = knownStructure;  % Copy the structure with NaNs
end

%  current_frame = all_stats{1,1};
%  for j = 1 : size(current_frame,1)
%          frame = all_stats{1,1}(j,:);
%          new_stats{1,1}(j,:).Area=frame.Area;
%          new_stats{1,1}(j,:).Centroid=frame.Centroid;
%          new_stats{1,1}(j,:).Eccentricity=frame.Eccentricity;
%          new_stats{1,1}(j,:).PixelIdxList=frame.PixelIdxList;
%          new_stats{1,1}(j,:).Intensity=frame.Intensity;
%          new_stats{1,1}(j,:).EquivDiameter=frame.EquivDiameter;
%          new_stats{1,1}(j,:).XC=frame.XC;
%          new_stats{1,1}(j,:).YC=frame.YC;
%  end
 new_stats{1,index_max}=all_stats{1,index_max};   


    for i = index_max-1:-1:1
        current_frame = all_stats{1,i};
        prev_frame = all_stats{1,i+1};
        for k = 1 : size(prev_frame,1)
            for j = 1 : size(current_frame,1)
                    if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        for l = 1 : max_particles
                           if  pdist2([all_stats{1,index_max}(l,:).Centroid],[current_frame(j,:).Centroid]) < total_max_movement
                          frame = all_stats{1,i}(j,:);
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
        current_frame = all_stats{1,i};
        prev_frame = all_stats{1,i-1};
        for k = 1 : size(prev_frame,1)
            for j = 1 : size(current_frame,1)
                    if abs(current_frame(j,:).XC-prev_frame(k,:).XC) < max_movement_allowed_per_frame ...
                            && abs(current_frame(j,:).YC-prev_frame(k,:).YC) < max_movement_allowed_per_frame
                        for l = 1 : max_particles
                           if  pdist2([all_stats{1,index_max}(l,:).Centroid],[current_frame(j,:).Centroid]) < total_max_movement
                          frame = all_stats{1,i}(j,:);
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

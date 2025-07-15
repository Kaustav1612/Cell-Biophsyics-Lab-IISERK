function [new_stats] = update_all_stats(all_stats, params)
    % Initialize new_stats as a copy of all_stats
    new_stats = all_stats;
    start_frame = params.FrameStart;
    end_frame = params.FrameEnd;
    
    % Add MergeParents and Parent fields if missing
    for i = start_frame:end_frame
        if ~isempty(new_stats{i})
            for k = 1:length(new_stats{i})
                if ~isfield(new_stats{i}(k), 'MergeParents')
                    new_stats{i}(k).MergeParents = [];
                end
                if ~isfield(new_stats{i}(k), 'Parent')
                    new_stats{i}(k).Parent = [];
                end
                % Precompute equivalent diameter for motion threshold
                new_stats{i}(k).EquivDiameter = 2 * sqrt(new_stats{i}(k).Area / pi);
            end
        end
    end
    
    % Forward tracking 
    for i = start_frame:end_frame-1
        current_frame = new_stats{i};
        next_frame_raw = new_stats{i+1};
        

        % Create empty template structure with all fields
        all_fields = union(fieldnames(new_stats{2}(1)), fieldnames(new_stats{2}(1)));
        emptyStruct = struct();
        for f = 1:numel(all_fields)
            emptyStruct.(all_fields{f}) = NaN;
        end
        emptyStruct.Parent = [];
        emptyStruct.MergeParents = [];
        emptyStruct.Centroid_X = [];
        emptyStruct.Centroid_Y = [];
        
        % Initialize next_frame with proper structure
        next_frame = repmat(emptyStruct, size(next_frame_raw));
        for k = 1:length(next_frame_raw)
            next_frame(k) = copy_struct_fields(next_frame_raw(k), emptyStruct);
            % Ensure EquivDiameter exists in next frame
            if ~isfield(next_frame(k), 'EquivDiameter')
                next_frame(k).EquivDiameter = 2 * sqrt(next_frame(k).Area / pi);
            end
        end

          
         % --- MODIFIED MATRIX CALCULATION AND NEW SCORE MATRIX ---
        overlap_matrix = zeros(length(current_frame), length(next_frame));
        distance_matrix = zeros(length(current_frame), length(next_frame)); % Stores raw distances
        score_matrix = Inf(length(current_frame), length(next_frame)); % Lower score is better

        % Parameters for scoring
        overlap_weight = 0.5; % How much to weigh overlap
        distance_weight = 1; % How much to weigh distance
        

        for j = 1:length(current_frame)
            if is_valid_object(current_frame(j))
                centroid_j = current_frame(j).Centroid;
                if current_frame(j).Area > 35
                    max_allowed_dist = current_frame(j).EquivDiameter+current_frame(j).EquivDiameter;
                    min_overlap_pixels = 1; % Minimum required overlap
                else
                    max_allowed_dist = current_frame(j).EquivDiameter+current_frame(j).EquivDiameter/2;
                    min_overlap_pixels = 5; % Minimum required overlap
                end
                
                for k = 1:length(next_frame)
                    if is_valid_object(next_frame(k))
                        centroid_k = next_frame(k).Centroid;
                        actual_distance = norm(centroid_j - centroid_k);
                        
                        if actual_distance <= max_allowed_dist % Motion constraint
                            current_overlap = numel(intersect(...
                                current_frame(j).PixelIdxList, ...
                                next_frame(k).PixelIdxList));
                            
                            overlap_matrix(j,k) = current_overlap;
                            distance_matrix(j,k) = actual_distance;
                            
                            % Calculate score: Lower is better.
                            % Inverse overlap (add a small epsilon to avoid div by zero)
                            % Normalized distance (e.g., by max_allowed_dist to make it unitless)
                            if current_overlap >= min_overlap_pixels
                                score_matrix(j,k) = (overlap_weight * (1 / (current_overlap + eps))) + ...
                                                    (distance_weight * (actual_distance / max_allowed_dist));
                            end
                        else
                            % If movement too large, mark with infinite score
                            overlap_matrix(j,k) = 0;
                            distance_matrix(j,k) = Inf;
                            score_matrix(j,k) = Inf;
                        end
                    end
                end
            end
        end
        
        % --- Matching Logic using Score Matrix ---
        updated_next = repmat(emptyStruct, 0, 1);
        matched_indices_next = false(1, length(next_frame));
        matched_indices_current = false(1, length(current_frame));

        % 1. Prioritize Merge Detection
        for k = 1:length(next_frame)
            if is_valid_object(next_frame(k)) && ~matched_indices_next(k)
                % Find all potential parents with a finite (valid) score
                potential_parents_indices = find(isfinite(score_matrix(:, k)));

                if numel(potential_parents_indices) > 1
                    % This is a merge event
                    new_obj = emptyStruct;
                    new_obj = copy_struct_fields(next_frame(k), new_obj);
                    
                    % Sort parents by score (ascending, as lower score is better)
                    parent_scores = score_matrix(potential_parents_indices, k);
                    [~, sort_idx] = sort(parent_scores, 'ascend'); 
                    
                    % The 'Parent' field is often the one with the best score (lowest cost)
                    % MergeParents includes all significant parents, sorted by score
                    new_obj.Parent = potential_parents_indices(sort_idx(1));
                    new_obj.MergeParents = potential_parents_indices(sort_idx); 

                    updated_next = [updated_next; new_obj];
                    matched_indices_next(k) = true;
                    matched_indices_current(potential_parents_indices) = true; 
                end
            end
        end

        % 2. Handle 1-to-1 Matches
        % Iterate through each object in current_frame to find its best child
        for j = 1:length(current_frame)
            if is_valid_object(current_frame(j)) && ~matched_indices_current(j) 
                % Find potential children in next_frame that haven't been matched yet
                % and have a finite score (valid connection)
                valid_children_indices = find(isfinite(score_matrix(j, :)) & ~matched_indices_next);
                
                if ~isempty(valid_children_indices)
                    % Sort potential children by score (ascending, as lower score is better)
                    child_scores = score_matrix(j, valid_children_indices);
                    [~, sort_idx] = sort(child_scores, 'ascend');
                    
                    % The best match is the one with the lowest score
                    best_match_k = valid_children_indices(sort_idx(1));
                    
                    % This is a 1-to-1 match
                    new_obj = emptyStruct;
                    new_obj = copy_struct_fields(next_frame(best_match_k), new_obj);
                    new_obj.Parent = j;
                    new_obj.MergeParents = []; 
                    
                    updated_next = [updated_next; new_obj];
                    matched_indices_next(best_match_k) = true;
                    matched_indices_current(j) = true;          
                end
            end
        end
        
        % 3. Identify New Objects (Appearances)
        for k = 1:length(next_frame)
            if is_valid_object(next_frame(k)) && ~matched_indices_next(k)
                new_obj = emptyStruct;
                new_obj = copy_struct_fields(next_frame(k), new_obj);
                new_obj.Parent = [];        
                new_obj.MergeParents = [];  
                updated_next = [updated_next; new_obj];
                matched_indices_next(k) = true; 
            end
        end
        
        % Store updated frame
        new_stats{i+1} = updated_next;
    end
end

% Helper function to safely copy fields (same as before)
function dest = copy_struct_fields(src, dest)
    if isempty(src) || ~isstruct(src)
        return;
    end
    src_fields = fieldnames(src);
    for f = 1:numel(src_fields)
        try
            dest.(src_fields{f}) = src.(src_fields{f});
        catch
            % Skip if field can't be copied
        end
    end
end

% Validation function (updated to include Centroid for robustness)
function valid = is_valid_object(obj)
    required_fields = {'Area', 'PixelIdxList', 'Centroid'};
    valid = isstruct(obj) && all(isfield(obj, required_fields)) && ...
            ~isnan(obj.Area) && ~isempty(obj.PixelIdxList) && ~any(isnan(obj.Centroid));
end
function tracked_objects = track_object_lifetimes(new_stats, nframes,params)
    start_frame = params.FrameStart;
    end_frame = params.FrameEnd;
    
    object_id = 1;
    tracked_objects = struct('id', {}, 'start_frame', {}, 'end_frame', {}, ...
                             'ancestry', {}, 'lifetime', {}, 'merged', {},'XC',{},'YC',{});

    % Track if objects are already assigned
    object_assigned = cell(nframes, 1);
    for i = start_frame:end_frame
        object_assigned{i} = false(1, length(new_stats{i}));
    end

    for frame = nframes:-1:1
        frame_objects = new_stats{frame};

        for obj_idx = 1:length(frame_objects)
            if object_assigned{frame}(obj_idx)
                continue;  % Already assigned especially when merge occurs
            end

            ancestry = zeros(1, nframes);
            ancestry(frame) = obj_idx;
            start_frame = frame;
            end_frame = frame;
            merged_flag = false;

            current_frame = frame;
            current_idx = obj_idx;
            XC=[];
            YC=[];

            while current_frame > 1
                current_obj = new_stats{current_frame}(current_idx);
                XC  =   [XC,current_obj.Centroid(1)];
                YC  =   [YC,current_obj.Centroid(2)];
                % Check for merge
                if ~isempty(current_obj.MergeParents)
                    merged_flag = true;
                    for p = current_obj.MergeParents
                        object_assigned{current_frame - 1}(p) = true;
                    end
                    break;
                elseif ~isempty(current_obj.Parent)
                    parent_idx = current_obj.Parent;
                    current_frame = current_frame - 1;
                    ancestry(current_frame) = parent_idx;
                    object_assigned{current_frame}(parent_idx) = true;
                    start_frame = current_frame;
                    current_idx = parent_idx;

                else
                    break;
                end
            end

            for f = start_frame:end_frame
                if ancestry(f) > 0
                    object_assigned{f}(ancestry(f)) = true;
                end
            end

            % Save to tracked_objects
            tracked_objects(end+1).id = object_id;
            tracked_objects(end).start_frame = start_frame;
            tracked_objects(end).end_frame = end_frame;
            tracked_objects(end).ancestry = ancestry(ancestry > 0);
            tracked_objects(end).lifetime = end_frame - start_frame + 1;
            tracked_objects(end).merged = merged_flag;
            tracked_objects(end).XC = XC;
            tracked_objects(end).YC = YC;

            object_id = object_id + 1;
        end
    end
end

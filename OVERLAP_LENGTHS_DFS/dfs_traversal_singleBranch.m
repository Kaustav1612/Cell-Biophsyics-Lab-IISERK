function [all_lengths, lengths, visited] = dfs_traversal_singleBranch(image, startPixel, visited, all_lengths,branch_points,start_points)
    % Initialize variables
    lengths = 0;
    [rows, cols] = size(image);
    directionality = 0; % Current direction of movement
       % Number of direction changes allowed
    stack = startPixel; % Initialize stack with the starting pixel
    prev_directionality =0;
    if ~isempty(startPixel)
        visited(startPixel(1), startPixel(2)) = 1;

        while ~isempty(stack) || length(find(visited==0))>0
            if isempty(stack) && length(find(visited==0))>0
                [not_visited_rows,not_visited_cols] = find(visited==0);
                non_visited = [not_visited_rows(1),not_visited_cols(1)];
                stack = [stack;non_visited];
               
            else
            currentPixel = stack(end, :);
            stack(end, :) = []; % Pop the last element

            % Define possible neighbors (8-connectivity)
            neighbors = [
                currentPixel(1) , currentPixel(2);    % Down
                currentPixel(1) + 1, currentPixel(2) + 1; % Down-Right (Diagonal)
                currentPixel(1), currentPixel(2) + 1;    % Right
                currentPixel(1) - 1, currentPixel(2) - 1; % Up-Left (Diagonal)
                currentPixel(1) - 1, currentPixel(2);    % Up
                currentPixel(1) - 1, currentPixel(2) + 1; % Up-Right (Diagonal)
                currentPixel(1), currentPixel(2) - 1;    % Left
                currentPixel(1) + 1, currentPixel(2) - 1; % Down-Left (Diagonal)
            ];

            % Find valid neighbors
            validIndices = [];
            for k = 1:size(neighbors, 1)
                neighbor = neighbors(k, :);
                if neighbor(1) > 0 && neighbor(1) <= rows && ...
                   neighbor(2) > 0 && neighbor(2) <= cols && ...
                   image(neighbor(1), neighbor(2)) == 1 && ...
                   ~visited(neighbor(1), neighbor(2))
                    validIndices = [validIndices; k];
                end
            end

            % Handle branching or continuation
            if length(validIndices) >= 2
                for i = 1:length(validIndices)
                    neighbor = neighbors(validIndices(i), :);
                    stack = [stack; neighbor]; % Add to stack
                    visited(neighbor(1), neighbor(2)) = 1;
                    
                    % Update lengths and directionality
                    if directionality ~= 0
                        if validIndices(i) == directionality
                            prev_directionality=directionality;
                            lengths = lengths + 1;
                        elseif directionality ~= validIndices(i) && abs(directionality-validIndices(i))<=2
                            prev_directionality=directionality;
                            directionality = validIndices(i);
                            lengths = lengths + 1;
                            
                        elseif directionality ~= validIndices(i) && abs(directionality-validIndices(i))>2 
                            % New branch detected, initiate a recursive call
                         
%                                 prev_directionality
%                                 directionality
%                                 
%                                 validIndices
%                                 neighbor
%                                 currentPixel
                                lengths = lengths + 1;
                                startPixel = neighbor;
                                [all_lengths, ~, visited] = dfs_traversal_singleBranch(image, startPixel, visited, all_lengths,branch_points,start_points);   
                            end    
                        end

                end
            elseif length(validIndices) == 1
                
                % Continue in the same direction
                neighbor = neighbors(validIndices, :);
                prev_directionality=directionality;
                directionality = validIndices;
                stack = [stack; neighbor];
                visited(neighbor(1), neighbor(2)) = 1;
                lengths = lengths + 1;
                
            elseif isempty(validIndices) 
                directionality=prev_directionality;
                prev_directionality=0;
            end
            end
        end
    end
    if lengths > 5 
    % Append current branch length to all_lengths
    all_lengths = [all_lengths, lengths];
    end
end

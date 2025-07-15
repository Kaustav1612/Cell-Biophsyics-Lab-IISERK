        function [lengths,Branches] =  first_dfs_traversal(image,startPixel,visited)
            Branches=struct('startPixel',[],'visited',[],'directionality',[]);
            lengths = 0;
            % Initialize variables
            [rows, cols] = size(image);
            directionality = 0;
            stack=[];
            stack = [stack; startPixel];
            
            visited(startPixel(1), startPixel(2)) = 1;
            while ~isempty(stack)
                currentPixel = stack(end, :);
                stack(end, :) = [];  % Pop the last element

                % Check neighbors (8- possible connectivity)
                    neighbors = [
                                    currentPixel(1) + 1, currentPixel(2);    % Down
                                    currentPixel(1) + 1, currentPixel(2) + 1; % Down-Right (Diagonal)
                                    currentPixel(1), currentPixel(2) + 1;    % Right
                                    currentPixel(1) - 1, currentPixel(2) - 1  % Up-Left (Diagonal)
                                    currentPixel(1) - 1, currentPixel(2);    % Up
                                    currentPixel(1) - 1, currentPixel(2) + 1; % Up-Right (Diagonal)    
                                    currentPixel(1), currentPixel(2) - 1;    % Left
                                    currentPixel(1) + 1, currentPixel(2) - 1; % Down-Left (Diagonal)               
                                ];


    
                        indices=[];
                        for k = 1:size(neighbors, 1)
                            neighbor = neighbors(k, :);

                            % Check bounds and if the neighbor is within the object and not visited
                            if neighbor(1) > 0 && neighbor(1) <= rows && ...
                               neighbor(2) > 0 && neighbor(2) <= cols && ...
                               image(neighbor(1), neighbor(2)) == 1 && ...
                               ~visited(neighbor(1), neighbor(2))
                           
                                indices=[indices;k];
                   
                            end
                        end
                        indices
                   
                     if length(indices) > 1
                    
                                for i = 1:length(indices)
                                    defect = 0;
                                    counter=0;
                                    neighbor = neighbors(indices(i), :);
                                    stack = [stack; neighbor];
                                    visited(neighbor(1), neighbor(2)) = 1;
                                    if directionality ~=0
                                        if indices(i) == directionality
                                            lengths = lengths +1;
                                            counter=counter+1;
                                        end
                                         if directionality~=indices(i) && defect <= 5
                                            lengths = lengths +1;
                                            directionality=indices(i);
                                            defect = defect + 1;
                                            continue;
                                        elseif directionality~=indices(i) && defect > 5
                                            Branches =  branch_storage(neighbors(indices(i), :),visited,directionality,Branches);
                                            defect = 0;
                                         end
                                    end
                                end
                     
                             
                            elseif  length(indices) == 1   
                                directionality = indices;
                                stack = [stack; neighbor];

                                if neighbor(1)~=0 && neighbor(2)~=0
                                    visited(neighbor(1), neighbor(2)) = 1;
                                    lengths = lengths +1;
                                else
                                    continue;
                                end
                            end
            end  
        end

        
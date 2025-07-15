function [filtered_img_binary, filtered_img_norm] = filtration(image,alpha,bright_pixel_count,threshold_distance,roi_norm)
%  Making a random image with the same number of bright pixels as the initial thresholded image to simulate noise and use the degree to filter out the noise

    rand_img = zeros(size(image));
    i = randi([1 size(image,1)],1, bright_pixel_count);
    j = randi([1 size(image,2)],1, bright_pixel_count);

    for k = 1:bright_pixel_count
        rand_img(i(k),j(k)) = 1;
    end

    [row_l, col_l] = size(rand_img);
    x_axis_pixel_count = zeros(1, row_l);
    y_axis_pixel_count = zeros(1, col_l);

    for i = 1: row_l
        count = 0;
        for j = 1: col_l
            if rand_img(i,j) == 1
                count = count + 1;
            end
        end
        x_axis_pixel_count(1,i) = count;
    end

    for j = 1: col_l
        count = 0;
        for i = 1: row_l
            if rand_img(i,j) == 1
                count = count + 1;
            end
        end
        y_axis_pixel_count(1,j) = count;
    end
    figure();
    imshow(rand_img,[]);
    
    

    %%  Using the generated random image to make the graph and compute the average degree -> 
    % Then do the filtering with alpha = 4 and visualise the filtered ROI

    ramdomised_log_array = rand_img > 0;
    [row, col] = find(ramdomised_log_array>0);

    % Coverting the coordinates in nanometer length scale 
    x_rand = row .* 15;
    y_rand = col .* 15;

    % Now computing the pairwise distance heavy on the memory
    rand_distances = pdist([x_rand, y_rand]);
    adjacency_matrix_rand = squareform(rand_distances<threshold_distance);

    % Weighted graph can be made as made previously (if needed)
    % Here we make a unweighted and undirected graph
    graph_rand = graph(adjacency_matrix_rand);
    deg_rand = degree(graph_rand);
    avg_random_deg = uw_avg_deg_und(deg_rand);

    % Defining Filtering Threshold using some fraction degree of the noised image:
    bg_filter = floor(alpha * avg_random_deg);

    %% Before filtering computing the graph 

    [row, col] = find(image>0); 
    % Getting the coordinates in nanometer 
    x = row .* 15;
    y = col .* 15;

    % Computing pairwise distance 
    distances = pdist([x, y]);
    adjacency_matrix_unfil = squareform(distances);
    graph_unfiltered_roi = graph(adjacency_matrix_unfil);
    adjacency_matrix = squareform(distances <= threshold_distance);
    graph_dist_filtered_roi = graph(adjacency_matrix);

    figure();
    plot(graph_dist_filtered_roi);
    title('Graph dist filtered roi');
    
    figure();
    plot(graph_unfiltered_roi);
    title('Graph unfiltered roi');
    
    figure();
    plot(graph_rand);
    title('Graph random roi');
    
    
    n = length(x); % or n = length(Y);
    XYTable = NaN(n, 3); 
    %Calulate the 
    degree_of_graph = degree(graph_dist_filtered_roi);

    XYTable(:, 1) = x;
    XYTable(:, 2) = y;
    XYTable(:, 3) = degree_of_graph;
    T = array2table(XYTable,...
        'VariableNames',{'X','Y','Degree'});

    % Identify pixel pairs with degree below bg_filter then set the degree of  those pixel pairs to be 0
    % and the filtered pixel pairs are selected

    belowThreshold = T.Degree <= bg_filter;
    T.Degree(belowThreshold) = 0;
    FilteredT = T(T.Degree ~= 0, :);

    
    filtered_img_binary = zeros(size(image));
    filtered_img_norm = zeros(size(image));
    
    

    for k = 1:length(FilteredT.X)
        filtered_img_binary(FilteredT.X(k)/15,FilteredT.Y(k)/15) = 1;
        filtered_img_norm(FilteredT.X(k)/15,FilteredT.Y(k)/15)=roi_norm(FilteredT.X(k)/15,FilteredT.Y(k)/15);
    end

    [row, col] = find(filtered_img_binary); 
    % Getting the coordinates in nanometer 
    x = row .* 15;
    y = col .* 15;

    % Computing pairwise distance 
    distances_filt = pdist([x, y]);
    adjacency_matrix_fil = squareform(distances_filt);
    graph_filtered_roi = graph(adjacency_matrix_fil);

    
    figure();
    plot(graph_filtered_roi);
    title('Graph filtered roi');

    figure();
    imshow(filtered_img_norm,[])
    title('Filtered Image');
    close;
%     directory='D:\Kaustav\CODES SET\MS CLUSTERING ALGO\Space_Shape_Intensity\Analysis\Process_Images\FILTERED ROIS\Transfected\ZAF\CLUSTERED\CELL_06';
%     filename=sprintf('Filtered_ROI_%d.fig',m);
%     savefig(fullfile(directory,filename));
%     close;
       
    figure();
    imshow(filtered_img_binary,[]);
    title('Binary Image');
    close;
end


% Function for computing average unweighted degree measure for undirected
% and unweighted network graph
function [uw_avg_deg] = uw_avg_deg_und(vector)
    % Computes average degree of undirected unweighted network graph
    % Input: A row/col vector containing degree of all the nodes
    % Returns: Unweighted average degree of the graph 
    uw_avg_deg = mean(vector);
end

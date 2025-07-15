function [roi_coords,sideLength] = draw_multiple_square_rois(image, gfp_image)
    % Display the images side by side
    figure;
    imshowpair(image, gfp_image);
    title('Draw Perfect Square ROIs');
    hold on;

    % Initialize a cell array to store ROI coordinates
    roi_coords = {};
    sideLength = {};
    continue_drawing = true;

    while continue_drawing
       
        
        h = drawrectangle('InteractionsAllowed', 'all', 'Color', 'r');

        % Extract top-left corner and dimensions
        top_left = h.Position(1:2);
        width = h.Position(3);
        height = h.Position(4);

        % Ensure the shape is a square
        side_length = min(width, height);
        square_coords = [top_left; 
                         top_left + [side_length, 0];
                         top_left + [side_length, side_length];
                         top_left + [0, side_length];
                         top_left];  % Close the square

        % Display the square ROI
        square_plot = plot(square_coords(:, 1), square_coords(:, 2), 'r-', 'LineWidth', 2);

        % Confirm the ROI or redraw
        choice = questdlg('Do you want to keep this ROI or redraw?', ...
                          'Confirm ROI', ...
                          'Keep', 'Redraw', 'Keep');
        if strcmp(choice, 'Redraw')
            delete(h); % Remove the drawn rectangle
            delete(square_plot); % Remove the displayed red square
            continue;  % Restart the loop
        end

        % Store the square ROI coordinates
        roi_coords{end + 1} = square_coords;
        sideLength{end + 1} = side_length;

        % Ask if the user wants to draw another ROI
        choice = questdlg('Do you want to draw another ROI?', ...
                          'Draw More?', ...
                          'Yes', 'No', 'Yes');
        if strcmp(choice, 'No')
            continue_drawing = false;
        end
    end

    % Close the figure
    hold off;
    filename=sprintf('RoiOverlay_Transfected.fig');
%     savefig(fullfile(directory,filename));
%     close;


end

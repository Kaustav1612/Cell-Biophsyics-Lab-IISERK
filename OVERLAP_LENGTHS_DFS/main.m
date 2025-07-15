
clear all;
warning off;
clc;
%% Image Analysis 

% Declare global variables
global  image masked_normalized_image object_length ;


% Read the image
img = imread('D:\Kaustav\ActCtrl_001.tif');

% Filtering image using gaussian blur
gauss_blur=imgaussfilt(img,8,"FilterSize",[13 13]);
filtered_image = img-gauss_blur;

% Normalisation of the filtered image
min_value = min(filtered_image(:));  % Get minimum value of all elements
max_value = max(filtered_image(:));  % Get maximum value of all elements
normalized_image = double(filtered_image - min_value) / double(max_value - min_value);

%Binarising using threshold 
threshold = 0.05;    
binary_image = imbinarize(normalized_image, threshold);

% Erosion and Dilation
eroded_image = imerode(binary_image, strel('disk',1));
dilated_image=imdilate(eroded_image,strel('disk',1));
BW = dilated_image;
BW_1=dilated_image;

imshow(BW, []);  % Display the black and white image
hold on;

redrawing = true;
while redrawing
  % Allow the user to draw a freehand boundary
  mask_handle = drawfreehand('Color','r','LineWidth',2);
  
  % Option to redraw or exit
  choice = questdlg("Redraw mask?", "Freehand Mask", "Yes", "No", "Yes");
  
  % Update flag based on user choice
  redrawing = strcmp(choice, "Yes");
  
  % Clear the previous mask for redraw (if needed)
  if redrawing
    delete(mask_handle);
  end
end

% Create a binary mask from the freehand boundary
mask = mask_handle.createMask;
image = BW .* mask;
masked_normalized_image = normalized_image.*mask;
figure();
imshow(image,[]);
title('Masked Image');
close;

[n_actin_located, componentIndices, stats,componentAreas,stats_label,~] = detect(image,15,masked_normalized_image);
object_length=calculateObjectLength(image,componentIndices,stats_label);
figure();
object_length= object_length(object_length>2);
object_length=object_length';
histogram(object_length.*0.065,30);
figure();
boxplot(object_length.*0.065)
title('Length Distribution of Actin in micro-meters');
xlabel('Length (micro-m)');
ylabel('Number');


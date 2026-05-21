# Cell-Biophsyics-Lab-IISERK

These scripts are written during my masters in IISER K at Cell Biophysics Lab, primarily related to cell imaging object detection and object tracking.

These scripts are written for a particular kind of analysis work, however these scripts can be tweaked to be used for analysis of biological images.

For any assistance feel free to contact me 

Email: kaustavpal1612@gmail.com
X: https://x.com/KaustavPal1612
Bluesky: https://bsky.app/profile/kaustavpal1612.bsky.social
LinkedIn: [https://www.linkedin.com/in/kaustav-pal-b524a5223/](https://www.linkedin.com/in/kaustavpal1612/)

Keep all the constituent functions in an individual folder all in the same directory 
Folder 1 and Folder 2 : DBSCAN CLUSTERING and MS CLUSTERING
Run the main.m file for DBSCAN or main_roi_updated.m for MS CLUSTERING in MATLAB

Put the paths in respective variable 
Put the alpha and threshold distance for filtration
Put the bandwidth value for the clustering

Folder 3: OVERLAP_LENGTHS_DFS
This is a MATLAB program to calculate the lengths of detected objects, targetted to understand distribution of actin fibres or other elongated structures
Run the main.m file in MATLAB adjust the detection by altering the sensitivity of binarisation

Folder 4: TIRF Movie
This is a MATLAB program to detect and track objects in high resolution movies, these program provides the a framework to classify objects as fixed or mobile. Run the file main.m with the movie path with TIFF image sequence.
Do provide the line-scan range based on the size of the object of interest, as well the correlation threshold on how focused the objects are. Also another parameter threshold calibration, to account for defocussing issues or bleaching. (Run main.m)

However, this code can be tweaked to use custom detection based on your method, find binary_image_stack replace that with the binary_stack of the custom detection, then the program analyzes trajectories based on that detection. 

Moreover it provides the diffusion charactersistics of mobile objects, for that run Traj_Analysis.m, collate the trajectories in a single mat file to get statistics for that run Collate_Traj.m, we look at the rolling time window to characterise diffusivity, also we look at the z direction diffusivity through intensity fluctuations.

Folder 5: Tracking Particle Spin Disk 
This MATAB is also similiar to particle or object tracking however it was custom written for larger objects with faster rates of movements, done for a collaborative work for Morphogenesis Lab.

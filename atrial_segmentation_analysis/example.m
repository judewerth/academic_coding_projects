%% Example to show how we read in .nrrd files
clear
clc
close all

% read in mri scan
%[scan, scan_info] = nrrdread('NRRD_Files/Atrial Wall/Segmentation_010_AW.nrrd');
scan = nrrdread ('NRRD_Files/Atrial Wall/Segmentation_010_AW.nrrd');
%'Example_nrrds/AF20_LGE_MRI.nrrd'
% scan contains the data from the MRI scan in a 3D matrix
% we can view it here:
scan = scan .* 300;
slice = 25;
figure(1)
image(scan(:,:,slice))
colormap(gca,gray(256)) 
numel( find(scan(:,:,25) == 1)) 
%%
% the second output of the nrrdread function is a structure containing 
% metadata from the nrrd
scan_info
% we can use this to obtain information about the image, such as voxel
% dimension
scan_info.spacings

% we can also read in our segmentations after saving them as nrrd files
% read in binary mask (segmentation)
[mask, mask_info] = nrrdread('Example_nrrds/AF20_Endo.nrrd');
figure(2)
image(mask(:,:,slice))
colormap(gca,gray(2))
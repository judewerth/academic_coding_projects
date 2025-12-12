clear
clc

slice = 25;

% Use nrrd function to get vectors
scan = nrrdread ('NRRD_Files/Jude/Atrial Wall/Segmentation_010_AW.nrrd');
BP = nrrdread ('NRRD_Files/Jude/Blood Pool/Segmentation_010_BP.nrrd');
MRI = nrrdread('NRRD_Files/Jude/MRI Scan/Segmentation_010_MRI.nrrd');
MRI = MRI(:,:,slice); %????
scan = uint16(scan);
BP = uint16(BP);

[fibrosis , intensity_array] = FibrosisFunction(BP , scan, MRI);

scan = scan(:,:,slice) .* 100;


intensity_slice = intensity_array(:,:,slice);
intensity = intensity_slice > 1.2;
intensity = uint16(intensity) .* 200;

final_scan = scan + intensity;

figure(1)
%scan = scan(200:400 , 300:500 , :);
image(final_scan);
colormap(gca,gray(256)) 
% show original MRI file underneath
% use d variable to make one array

figure(2)
image(MRI)
colormap(gca,gray(256)) 
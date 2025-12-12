function [fibrosis , Intensity_Ratio] = ...
    FibrosisFunction(Blood_Pool , Atrial_Wall , MRI , threshold , participants)

fibrosis = zeros(1 , participants);

for i = 1:participants

    % find number of voxyls 
    BP_count = sum(Blood_Pool{i}(:));
    AW_count = sum(Atrial_Wall{i}(:));

    % find MRI values of marked voxyls
    BP_values = Blood_Pool{i} .* MRI;
    AW_values = Atrial_Wall{i} .* MRI;
    
    % find average MRI value of marked voxyls
    BP_mean = (sum(BP_values(:)) / BP_count);
    
    % find ratio of Atrial Wall Values / Blood Pool Mean
    Intensity_Ratio = AW_values ./ BP_mean;
    
    % find number of voxyls over the set threshold
    threshold_pixels = numel( find(Intensity_Ratio >= threshold) );
    % threh = sum(Intensity_Ratio >= threshold)
    % find the percentage of voxyls over the treshold (in Atrial Wall)
    fibrosis(i) = threshold_pixels / AW_count;

end








clear
clc

% Datapath Information
% folder names
global_folder = "nrrd files/";
names = ["Jude/" , "Ben/" , "Eric/"];
folders = ["Left Atrium/" ; "Blood Pool/" ; "Atrial Wall/"];
% file names
file_numbers = ["001" "002" "003" "004" "005"  "007" "008" "009" "010"...
                "011" "012"  "014" "016" "017" "018" "019" "021" "022"...
                "025" "026" "027" "028" "045"];
% missing "006" and "013"
suffix = ["_LA.nrrd" ; "_BP.nrrd" ; "_AW.nrrd"];

% Parameter Information
r = length(file_numbers); % number of files (usually length(file_numbers))
dimensions = [.625 , .625 , 2.5]; % size of each segmentation voxel (mm^3)
threshold = 1.2; % intensity threshold to quantify a voxel as fibrotic
threshold_vector = .98 : .01 : 1.6;

% Size Abreviations
x = length(folders); % number of different types of segmentations (3)
y = length(names); % number of participants (3)

% Zero Vectors
data_array = cell(x , y); % type x name
volume = zeros(r , y); % file x name
fibrosis = zeros(r , y, length(threshold_vector)); % file x name
dice = zeros(r , (y*(y-1)/2)); % file x name pairs


for i = 1:r % loop which iterates through files

    for k = 1:length(names) % loop which iterates through names (to generate datapaths)

        initial = char(names(k));
        initial = initial(1);

        datapath(:,k) = strcat ( ...
        global_folder , names(k) , folders(:) , file_numbers(i) , initial , suffix(:) );
        % 3 x 3 (array of every datapath for a file) oragnized by type x name
        %   Jude Ben Eric
        % LA
        % BP
        % AW
        % could use advice to make this more efficient
        % converting between char and str , having "zero" vector

    end

    for n = 1:y % loops through participants
        for m = 1:x % loops through segmentation types

            data_array{m , n} = nrrdread( datapath(m , n) );
            % creates data array for names and segmentation types
            % (same organization as datapath array)
            data_array{m , n} = uint16( data_array{m , n} );
            % converts data to universal unit of uint16

        end
    end

    % MRI data
    MRI = nrrdread(strcat(...
        global_folder , "MRI Scans/" , file_numbers(i) , "_MRI.nrrd")); 

    % Volume
    volume(i , :) = VolumeFunction(data_array(1 , :) , dimensions , y);
    % volume = [Jude Ben Eric]

    for h = 1:length(threshold_vector)
        % Fibrosis Rate
        fibrosis(i , : , h) = FibrosisFunction(data_array(2 , :) , data_array(3 , :) , MRI , threshold_vector(h) , y);
        % fibrosis = [Jude Ben Eric]
        if sum(fibrosis(i , : , h)) <= .01
            break
        end
    end

    % Dice Coeffefient
    dice(i , :) = DiceFunction(data_array(1 , :) , y);
    % dice = [ Ben & Eric , Jude & Eric , Jude & Ben ]

end

quality_scores = [4 4 4 4 4  3 2 4 2 2  3 3 3 2 3 4 4 3 3 3 3 3 2];

%   file number : 01 02 03 04 05  07 08 09 10 11 12  14 16 17 18 19 21 22 25 26 27 28 45
% quality score : 4  4  4  4  4   3  2  4  2  2  3   3  3  2  3  4  4  3  3  3  3  3  2

data_stats = zeros(x , 2 ,3); 
% 3D vector :           mean | standard deviation
%             2 quality
%             3 quality
%             4 quality
          
% volume (difference) -- fibrosis (difference) -- dice coefficient


for i = 1:3

    index = find(quality_scores == i+1); 
    % finds the files with the specific quality scores (2 , 3 , 4)

    quality_data = [ mean( abs( diff( volume(index , :)' ) ) ) % mean difference for volume
                     mean( abs( diff( fibrosis(index, :)' ) ) ) % mean difference for fibrosis
                     mean( dice(index , :) , 2)' ];  % mean of dice coefficient

    data_stats(i,1,:) = mean( quality_data , 2 );
    data_stats(i,2,:) = std( quality_data , 0 , 2);
    
    fig = figure(i);
    graphs = tiledlayout( ceil( length(index)/2 ) , 2 , 'TileSpacing','tight' , 'Padding','tight' );
    labels = file_numbers(index);
    tit = title( graphs , strcat("Fibrosis Graphs for " , string(i+1) , " Quality Imaages"));

    for j = 1:length(index)

        clear("fibrosis_values")
        fibrosis_values(:,:) = fibrosis(index(j) , : , :);
        fibrosis_values = fibrosis_values';

        nexttile
        plot(threshold_vector , fibrosis_values )

        xlabel("Threshold Value")
        ylabel("Fibrosis")
        title(labels(j))


    end

    leg = legend( 'Jude' , 'Ben' , 'Eric');
    leg.Orientation = 'horizontal';
    leg.Layout.Tile = 'north';

end


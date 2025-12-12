function volume = VolumeFunction(Left_Atrium , dimensions , participants)
% Just input the name of the folder or the data path of the folder (use '')

volume = zeros(1, participants);

for i = 1:participants

    % amount of voxyls marked in file
    voxel = sum(Left_Atrium{i}(:));
    
    % find volume
    volume(i) = voxel * prod(dimensions) * .001; 
    % multiply number of voxyls by the size each voxyl for total size of file
    % going from mm^3 --> cm^3 (* by .001)

end



% find the dimensions of each voxyl
% a = LA_info.spacedirections;
% b = find(a == ',');
% dimensions = [str2double(a(2:b(1)-1)) , ...
%         str2double(a(b(3)+1:b(4)-1)) , ...
%         str2double(a(b(6)+1:end-1))];
% I wanted it to work so it would find the dimension regardless of the
% lengths of the dimensions (which is why it's weirld complicated)
% Seg3D does size = (x,0,0)(0,y,0)(0,0,z) 
% Coreview does size = x y z



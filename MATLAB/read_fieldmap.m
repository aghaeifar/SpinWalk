function [fieldmap, mask, fov] = read_fieldmap(filename)


fileID = fopen(filename);
dim = fread(fileID, 3, 'int32=>int32'); %this is the matrix size
fov = fread(fileID, 3, 'single=>single');% FOV in [m]
fieldmap = fread(fileID, prod(dim), 'single=>single'); % dB0 map
mask     = fread(fileID, prod(dim), 'uint8=>uint8'  ); % mask of vessels
fclose(fileID);

fieldmap = reshape(fieldmap, dim');
mask = reshape(mask, dim');

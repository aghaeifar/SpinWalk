function [fieldmap, mask, fov] = read_fieldmap(filename)

fileID = fopen(filename);
fseek(fileID, 0, 'eof');
filesize = ftell(fileID);
fclose(fileID);


fileID = fopen(filename);
dim = fread(fileID, 3, 'int32=>int32'); %this is the matrix size
fov = fread(fileID, 3, 'single=>single');% FOV in [m]

fieldmap = [];
mask     = [];
header_size = 6 * 4;
if filesize == (prod(dim) * 5 + header_size)
    fieldmap = fread(fileID, prod(dim), 'single=>single'); % dB0 map
    mask     = fread(fileID, prod(dim), 'uint8=>uint8'  ); % mask of vessels
    fieldmap = reshape(fieldmap, dim');
    mask     = reshape(mask, dim');
elseif filesize == (prod(dim) + header_size)
    mask     = fread(fileID, prod(dim), 'uint8=>uint8'  ); 
    mask     = reshape(mask, dim');
end
fclose(fileID);



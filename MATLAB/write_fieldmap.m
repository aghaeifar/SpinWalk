function write_fieldmap(fieldmap, mask, fov, filename)

fileID = fopen(filename, 'w');

fwrite(fileID, size(mask), 'int32'); %this is the matrix size
fwrite(fileID, single(fov), 'single');% FOV in [m]
if isempty(fieldmap) == false
    fwrite(fileID, single(fieldmap), 'single'); % dB0 map
end
fwrite(fileID, uint8(mask), 'uint8'  ); % mask of vessels
fclose(fileID);


function write_fieldmap(fieldmap, mask, fov, filename)

fileID = fopen(filename, 'w');
fwrite(fileID, size(fieldmap), 'int32'); %this is the matrix size
fwrite(fileID, single(fov), 'single');% FOV in [m]
fwrite(fileID, single(fieldmap), 'single'); % dB0 map
fwrite(fileID, uint8(mask), 'uint8'  ); % mask of vessels
fclose(fileID);


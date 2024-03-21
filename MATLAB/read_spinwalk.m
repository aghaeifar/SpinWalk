function [m_xyz, dims, hdr_extra] = read_spinwalk(filename)


fileID = fopen(filename);
fseek(fileID, 0, 'eof');
filesize = ftell(fileID);
fclose(fileID);


fileID      = fopen(filename);
hdr_size    = fread(fileID, 1, 'int32=>int32');
dims        = fread(fileID, 4, 'int32=>int32'); % 3 * n_echo * n_spins * n_sample_length_scales
hdr_extra   = fread(fileID, (hdr_size - 4*numel(dims))/8, 'double=>double');

% guess data type
datasize = filesize - hdr_size - 4;
if datasize == prod(dims)
    m_xyz   = fread(fileID, prod(dims), 'uint8=>uint8');
elseif datasize == 4*prod(dims)
    m_xyz   = fread(fileID, prod(dims), 'single=>single');
elseif datasize == 8*prod(dims)
    m_xyz   = fread(fileID, prod(dims), 'double=>double');
else
    error('Data type not recognized!')
end

fclose(fileID);
m_xyz       = reshape(m_xyz, dims(:)');
m_xyz       = reshape(m_xyz, [dims(1), dims(2), dims(3), dims(4)]);



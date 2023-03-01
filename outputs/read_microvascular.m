function [m_xyz, dims, hdr_extra] = read_microvascular(filename)


fileID = fopen(filename);
hdr_size  = fread(fileID, 1, 'int32=>int32');
dims      = fread(fileID, 4, 'int32=>int32'); % 3 * n_spins * n_sample_length_scales * device_count
hdr_extra = fread(fileID, hdr_size/4-numel(dims), 'single=>single');
m_xyz = fread(fileID, prod(dims), 'single=>single');
fclose(fileID);
m_xyz = reshape(m_xyz, dims(:)');
m_xyz = reshape(m_xyz, [dims(1), dims(2)*dims(3), dims(4)]);

clear
fname = 'results.dat';

fileID = fopen(fname);
n_spins         = fread(fileID, 1, 'int32=>int32');
n_fieldmaps     = fread(fileID, 1, 'int32=>int32');
n_sample_length = fread(fileID, 1, 'int32=>int32');
device_count    = fread(fileID, 1, 'int32=>int32');
m_xyz = fread(fileID, 3 * n_spins * n_fieldmaps * n_sample_length * device_count, 'single=>single');
fclose(fileID);
m_xyz = reshape(m_xyz, [3, n_spins, n_fieldmaps, n_sample_length, device_count]);
m_xyz = permute(m_xyz, [1 2 5 3 4]);
m_xyz = reshape(m_xyz, [3, n_spins*device_count, n_fieldmaps, n_sample_length]);

figure
spins_xy = squeeze(complex(sum(m_xyz(1,:,:,:), 2), sum(m_xyz(2,:,:,:), 2) ));
signal_magnitude = abs(spins_xy);
relative_signal = 100 * (1 - signal_magnitude(1,:)./ signal_magnitude(2,:));

disp(relative_signal')

vessel_radius = [0.5 0.8 1 1.5 2 4 5.5 8 10 12 16 20 32 45 64 75 90 100 110 128 140 180 256 400 512 620 800 1024];

semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');

%% 


clear
clc
rad_ref_um = 53.367;
fname{1} = './results_gre_fieldmap1.dat';
fname{2} = './results_gre_fieldmap2.dat';

fname{1} = './results_ssfp_fieldmap1.dat';
fname{2} = './results_ssfp_fieldmap2.dat';

spins_xy = [];
for i=1:numel(fname)
    fileID = fopen(fname{i});
    hdr_size  = fread(fileID, 1, 'int32=>int32');
    dims      = fread(fileID, 4, 'int32=>int32'); % 3 * n_spins * n_sample_length_scales * device_count
    if dims(3) ~= hdr_size/4-numel(dims)
        warning('Why header info is confusing here?')
    end
    scales    = fread(fileID, hdr_size/4-numel(dims), 'single=>single');
    m_xyz = fread(fileID, prod(dims), 'single=>single');
    fclose(fileID);
    m_xyz = reshape(m_xyz, dims(:)');
    m_xyz = permute(m_xyz, [1 2 4 3]);
    m_xyz = reshape(m_xyz, [dims(1), dims(2)*dims(4), dims(3)]);
    spins_xy = cat(4, spins_xy, m_xyz);
end

spins_xy = squeeze(complex(sum(spins_xy(1,:,:,:), 2), sum(spins_xy(2,:,:,:), 2) ));

%figure
signal_magnitude = abs(spins_xy);
relative_signal = 100 * (1 - signal_magnitude(:,1)./ signal_magnitude(:,2));

disp(relative_signal')

vessel_radius = rad_ref_um * scales;

semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');
hold on
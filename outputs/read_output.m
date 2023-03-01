clear
clc
rad_ref_um = 53.367;
fname{1} = './gre_m1_0.dat';
fname{2} = './gre_m1_1.dat';

% fname{1} = './results_ssfp_M1_0_fieldmap1.dat';
% fname{2} = './results_ssfp_M1_1_fieldmap2.dat';

% fname{1} = './results_se_M1_0_fieldmap1.dat';
% fname{2} = './results_se_M1_1_fieldmap2.dat';

spins_xy = [];
for i=1:numel(fname)
    [m_xyz, dims, scales] = read_microvascular(fname{i});
    if dims(4) ~= numel(scales)
        warning('Why header info is confusing here?')
    end
    spins_xy = cat(4, spins_xy, m_xyz);
end

spins_xy = squeeze(complex(sum(spins_xy(1,:,:,:), 2), sum(spins_xy(2,:,:,:), 2) ));
signal_magnitude = abs(spins_xy);
relative_signal  = 100 * (1 - signal_magnitude(:,1)./ signal_magnitude(:,2));
vessel_radius    = rad_ref_um * scales;
semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');
hold on

%% Read Boston data
clear
clc
folder = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/';
names = {'gre' , 'se', 'ssfp'};
n_fnames = 10;

for n=1:numel(names)

    spins_xy = [];    
    for i=0:n_fnames-1
        filename = fullfile(folder, [names{n} '_boston_m1_' num2str(i) '.dat']);
        [m_xyz, dims, scales] = read_microvascular(filename);
        spins_xy = cat(3, spins_xy, m_xyz);
    end
    
    spins_xy = squeeze(complex(sum(spins_xy(1,:,:), 2), sum(spins_xy(2,:,:), 2) ));
    signal_magnitude = abs(spins_xy);
    
    subplot(2,2,n);
    plot(signal_magnitude(1:5), 'b'); hold on
    plot(signal_magnitude(6:end), 'r');
    title(names{n});
end
clc
rad_ref_um = 53.367;
fname{1} = './gre_m1_0.dat';
fname{2} = './gre_m1_1.dat';

% fname{1} = './ssfp_m1_0.dat';
% fname{2} = './ssfp_m1_1.dat';

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

%% Read Boston data orientation
clc
folder = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs';
f      = dir(fullfile(folder, '*boston*.dat'));
names = {'gre', 'se', 'ssfp'};
rot   = linspace(0,90,6);
time  = [0, 1, 2, 3, 4];

spins_xy = [];
for seq=names
    spins_xy_r = [];
    for r=rot
        spins_xy_t = [];
        for t=time
            filename = fullfile(folder, [seq{1} '_boston_m1_' num2str(t) '_rot' num2str(r, '%02d') '.dat']);
            [m_xyz, dims, ~] = read_microvascular(filename);
            spins_xy_t = cat(3, spins_xy_t, m_xyz);            
        end
        spins_xy_r = cat(4, spins_xy_r, spins_xy_t);
    end
    spins_xy = cat(5, spins_xy, spins_xy_r);
end

spins_xy = squeeze(complex(sum(spins_xy(1,:,:,:,:), 2), sum(spins_xy(2,:,:,:,:), 2) ));
signal_magnitude = abs(spins_xy);

figure(1); 
for seq=1:numel(names)
    subplot(211); hold on
    signal = squeeze(signal_magnitude(:,1,seq));
    relative_signal  = 100 * (signal ./ signal(1) - 1);
    plot(time, relative_signal);
end
legend('GRE', 'SE', 'SSFP')

for seq=1:numel(names)
    subplot(212); hold on
    signal = squeeze(signal_magnitude(1,:,seq));
    relative_signal  = 100 * (signal ./ signal(1) - 1);
    plot(rot, relative_signal);
end
legend('GRE', 'SE', 'SSFP')




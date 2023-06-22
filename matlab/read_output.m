clc
rad_ref_um = 53.367;
fname{1} = './gre_m1_0.dat';
fname{2} = './gre_m1_1.dat';

% fname{1} = './ssfp_m1_0.dat';
% fname{2} = './ssfp_m1_1.dat';

fname{1} = './se_m1_0.dat';
fname{2} = './se_m1_1.dat';

fname{1} = './grase_m1_0.dat';
fname{2} = './grase_m1_1.dat';

spins_xy = [];
for i=1:numel(fname)
    [m_xyz, dims, scales] = read_microvascular(fname{i});
    if dims(4) ~= numel(scales)
        warning('Why header info is confusing here?')
    end
    spins_xy = cat(5, spins_xy, m_xyz);
end

spins_xy = squeeze(complex(sum(spins_xy(1,:,:,:,:), 3), sum(spins_xy(2,:,:,:,:), 3) ));
signal_magnitude = abs(spins_xy);
relative_signal  = 100 * (1 - signal_magnitude(:,:,1)./ signal_magnitude(:,:,2));

vessel_radius    = rad_ref_um * scales;
h = semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');
set(h, {'color'}, [num2cell(distinguishable_colors(size(relative_signal,1)), 2)]);


%% GRASE
fname{1} = './grase_m1_0.dat';
fname{2} = './grase_m1_1.dat';

TR = 0.2;
EcoSpc = 0.005;
vessel_radius    = rad_ref_um * scales;

spins_xy = [];
for i=1:numel(fname)
    [m_xyz, dims, scales] = read_microvascular(fname{i});
    if dims(4) ~= numel(scales)
        warning('Why header info is confusing here?')
    end
    spins_xy = cat(5, spins_xy, m_xyz);
end

spins_xy = squeeze(complex(sum(spins_xy(1,:,:,:,:), 3), sum(spins_xy(2,:,:,:,:), 3) ));
signal_magnitude = abs(spins_xy);
relative_signal  = 100 * (1 - signal_magnitude(:,:,1)./ signal_magnitude(:,:,2));
relative_signal  = signal_magnitude(:,:,2) - signal_magnitude(:,:,1);

imagesc(relative_signal');
ind = floor(linspace(1, TR/EcoSpc, 8));
set(gca,'XTick', ind, 'XTickLabel', EcoSpc * ind * 1e3, 'color','none');  
xlabel('Echo Time (ms)');
ind = floor(linspace(1, numel(scales), 6));
set(gca,'YTick', ind, 'YTickLabel', round(vessel_radius(ind)), 'color','none');  box off;
ylabel('Vessl radius (um)');
title('BOLD Signal')
colormap(inferno); 
% h = semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');
% set(h, {'color'}, [num2cell(distinguishable_colors(size(relative_signal,1)), 2)]);
% legend('Echo1 (GRE)', 'Echo2 (GRE)', ' Echo3 (SE)', 'Echo4 (GRE)', 'Echo5 (GRE)', 'Echo6 (SE)', 'Echo7 (GRE)')


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




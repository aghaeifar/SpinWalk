clear
clc

%% 
rad_ref_um = 53.367;
fname{1,1} = './gre_m1_0.dat';
fname{1,2} = './gre_m1_1.dat';
fname{2,1} = './se_m1_0.dat';
fname{2,2} = './se_m1_1.dat';
fname{3,1} = './ssfp_m1_0.dat';
fname{3,2} = './ssfp_m1_1.dat';
fname{4,1} = './grase_m1_0.dat';
fname{4,2} = './grase_m1_1.dat';

spins_xy = cell(size(fname,1), 1);
signal_magnitude = cell(size(fname,1), 1);
relative_signal = cell(size(fname,1), 2);
for i=1:size(fname,1)
    spins_xy{i} = [];
    for j=1:size(fname, 2)
        [m_xyz, dims, scales] = read_microvascular(fname{i,j});
        if dims(4) ~= numel(scales)
            warning('Why header info is confusing here?')
        end
        spins_xy{i} = cat(5, spins_xy{i}, m_xyz);
    end

    temp = sum(spins_xy{i}, 3);
    temp = complex(temp(1,:,:,:,:), temp(2,:,:,:,:));
    signal_magnitude{i} = abs(temp);
    relative_signal{i,1}  = 100 * (1 - signal_magnitude{i}(:,:,:,:,1)./ signal_magnitude{i}(:,:,:,:,2));
    relative_signal{i,2}  = signal_magnitude{i}(:,:,:,:,2) - signal_magnitude{i}(:,:,:,:,1);

    subplot(2,2,i);
    vessel_radius    = rad_ref_um * scales;
    if i ~= 4
        semilogx(vessel_radius, squeeze(relative_signal{i,1})); xlabel('Vessel radius (um)'); ylabel('Relative Signal %');
    else
        imagesc(squeeze(relative_signal{i,2})');
        ind = floor(linspace(1, size(relative_signal{i,2},2), 8));
        set(gca,'XTick', ind, 'XTickLabel', 0.005 * ind * 1e3, 'color','none'); xlabel('Echo Time (ms)');
        ind = floor(linspace(1, numel(vessel_radius), 5));
        set(gca,'YTick', ind, 'YTickLabel', round(vessel_radius(ind)), 'color','none');  box off; ylabel('Vessel radius (um)');
    end
end



%%


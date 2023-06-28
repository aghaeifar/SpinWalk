clc
rad_ref_um = 53.367;
fname{1} = '../../outputs/gre_m1_0.dat';
fname{2} = '../../outputs/gre_m1_1.dat';

fname{1} = '../../outputs/se_m1_0.dat';
fname{2} = '../../outputs/se_m1_1.dat';

fname{1} = '../../outputs/ssfp_m1_0.dat';
fname{2} = '../../outputs/ssfp_m1_1.dat';


dim_echo = 2;
dim_spin = 3;
dim_vessel_size = 4;
spins_xy = [];
for i=1:numel(fname)
    [m_xyz, dims, scales] = read_spinwalk(fname{i});
    if dims(4) ~= numel(scales)
        warning('Why header info is confusing here?')
    end
    spins_xy = cat(ndims(m_xyz)+1, spins_xy, m_xyz);
end

spins_xy = complex(sum(spins_xy(1,:,:,:,:), dim_spin), sum(spins_xy(2,:,:,:,:), dim_spin) );
signal_magnitude = abs(spins_xy);
relative_signal  = 100 * (1 - signal_magnitude(:,:,:,:,1)./ signal_magnitude(:,:,:,:,2));
relative_signal  = squeeze(relative_signal);

vessel_radius    = rad_ref_um * scales;
h = semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('BOLD Signal %');


%% GRASE
fname{1} = '../../outputs/grase_m1_0.dat';
fname{2} = '../../outputs/grase_m1_1.dat';

TR = 0.2;
EcoSpc = 0.005;
vessel_radius    = rad_ref_um * scales;

spins_xy = [];
for i=1:numel(fname)
    [m_xyz, dims, scales] = read_spinwalk(fname{i});
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





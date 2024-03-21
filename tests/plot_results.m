clc
clear

file_m1{1}{1} = '../../outputs/gre_m1_fieldmap_0.dat';
file_m1{1}{2} = '../../outputs/gre_m1_fieldmap_1.dat';
file_m1{2}{1} = '../../outputs/se_m1_fieldmap_0.dat';
file_m1{2}{2} = '../../outputs/se_m1_fieldmap_1.dat';
file_m1{3}{1} = '../../outputs/ssfp_m1_fieldmap_0.dat';
file_m1{3}{2} = '../../outputs/ssfp_m1_fieldmap_1.dat';

file_T{1}{1} = '../../outputs/gre_T_fieldmap_0.dat';
file_T{1}{2} = '../../outputs/gre_T_fieldmap_1.dat';
file_T{2}{1} = '../../outputs/se_T_fieldmap_0.dat';
file_T{2}{2} = '../../outputs/se_T_fieldmap_1.dat';
file_T{3}{1} = '../../outputs/ssfp_T_fieldmap_0.dat';
file_T{3}{2} = '../../outputs/ssfp_T_fieldmap_1.dat';

dim_xyz  = 1;
dim_echo = 2;
dim_spin = 3;
dim_vessel_size = 4;
rad_ref_um = 53.367;

tissue_type = 0; % 0 = extra-vascular, 1 = intra-vascular, [0,1] = combined

signal_magnitude = cell(numel(file_m1), numel(file_m1{1}));
for seq = 1:numel(file_m1)   
    for i=1:numel(file_m1{seq})
        [m1, dims, scales] = read_spinwalk(file_m1{seq}{i});
        T = read_spinwalk(file_T{seq}{i});        
        if dims(4) ~= numel(scales)
            warning('Header info is confusing here?')
        end
        
        m1_t = zeros(numel(scales), 1);
        for s=1:numel(scales)
            ind = ismember(T(1,end,:,s), tissue_type);
            m1_f = m1(:,:,ind(:),s);
            m1_t(s) = abs(complex(sum(m1_f(1,:)), sum(m1_f(2,:))));
        end
        signal_magnitude{seq, i} = m1_t;
    end
end

clf
vessel_radius  = rad_ref_um * scales;
for seq = 1:numel(file_m1)
    relative_signal  = 100 * (1 - signal_magnitude{seq, 1} ./ signal_magnitude{seq, 2});    
    h = semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('BOLD Signal %'); 
    hold on;
    ylim([0, 7])
end
legend('GRE', 'SE', 'SSFP')


%% grase sequence

% file_m1{4}{1} = '../../outputs/grase_m1_0.dat';
% file_m1{4}{2} = '../../outputs/grase_m1_1.dat';
TR = 0.2;
EcoSpc = 0.005;
%         subplot(2,1,2)
%         relative_signal  = signal_magnitude{seq}(:,:,:,:,2) - signal_magnitude{seq}(:,:,:,:,1);
%         relative_signal  = 100 * (1 - signal_magnitude{seq}(:,:,:,:,1)./ signal_magnitude{seq}(:,:,:,:,2));
%         relative_signal  = squeeze(relative_signal);
%         imagesc(relative_signal);
%         axis image
%         set(gca,'YDir','normal') 
% 
%         ind = floor(linspace(1, numel(scales), 6));
%         set(gca,'XTick', ind, 'XTickLabel', round(vessel_radius(ind)), 'color','none');  box off;
%         xlabel('Vessl radius (um)');
% 
%         ind = floor(linspace(1, TR/EcoSpc, 8));
%         set(gca,'YTick', ind, 'YTickLabel', EcoSpc * ind * 1e3, 'color','none');  
%         ylabel('Echo Time (ms)');
%         
% %         title('BOLD Signal')
%         colormap(inferno); 
%         colorbar

%% stimulated echo
 clc
% clear
file_m1{1} = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/ste_m1_0.dat';
file_m1{2} = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/ste_m1_1.dat';

dim_echo = 2;
dim_spin = 3;
dim_vessel_size = 4;
rad_ref_um = 53.367;


spins_xyz = [];
for i=1:numel(file_m1)
    [m_xyz, dims, scales] = read_spinwalk(file_m1{i});
    spins_xyz = cat(ndims(m_xyz)+1, spins_xyz, m_xyz);
end
vessel_radius    = rad_ref_um * scales;

spins_xy = complex(sum(spins_xyz(1,:,:,:,:), dim_spin), sum(spins_xyz(2,:,:,:,:), dim_spin) );
signal_magnitude = abs(spins_xy);

relative_signal  = 100 * (1 - signal_magnitude(:,:,:,:,1) ./ signal_magnitude(:,:,:,:,2));
difference_signal  = signal_magnitude(:,:,:,:,2) - signal_magnitude(:,:,:,:,1);

figure(1)
subplot(1,2,1)
cla
semilogx(vessel_radius, squeeze(relative_signal(1,1,1,:)));
hold on
semilogx(vessel_radius, squeeze(relative_signal(1,2,1,:))); xlabel('Vessel radius (um)'); ylabel('100 * (1 - S_r_e_s_t / S_a_c_t )'); 
hold off
legend('SE', 'STE')

subplot(1,2,2)
cla
semilogx(vessel_radius, squeeze(difference_signal(1,1,1,:)));
hold on
semilogx(vessel_radius, squeeze(difference_signal(1,2,1,:))); xlabel('Vessel radius (um)'); ylabel('S_a_c_t - S_r_e_s_t'); 
hold off
legend('SE', 'STE')

figure(2)
cla
plot(squeeze(signal_magnitude(:,2,:,:,1)))
hold on
plot(squeeze(signal_magnitude(:,2,:,:,2)))

plot(squeeze(signal_magnitude(:,1,:,:,1)))
plot(squeeze(signal_magnitude(:,1,:,:,2)))
hold off
legend('STE Rest', 'STE Act', 'SE Rest', 'SE Act')

%% random-walk trajectory
clear
clc;
filename = '../../outputs/gre_trajectory_xyz1_fieldmap_0.dat';
[xyz_all, dims, hdr_extra] = read_spinwalk(filename);
size(xyz_all)
xyz_all = xyz_all * 1e6;
n_spins = size(xyz_all, 3);
for s=round(linspace(1, n_spins, 10))
    xyz = xyz_all(:,:,s,2);
    plot3(xyz(1,:), xyz(2,:), xyz(3,:))
    hold on;
end
hold off

% nonzeros(xyz_all)

%% diffusion
clear
clc

file_fieldmap = '../../field_maps/restricted_diffusion_model.dat';

file_m1{1} = '../../outputs/dwi_base_m1_restricted_diffusion_model.dat';
file_m1{2} = '../../outputs/dwi_x_m1_restricted_diffusion_model.dat';
file_m1{3} = '../../outputs/dwi_y_m1_restricted_diffusion_model.dat';
file_m1{4} = '../../outputs/dwi_z_m1_restricted_diffusion_model.dat';

file_T{1} = '../../outputs/dwi_base_xyz1_restricted_diffusion_model.dat';
file_T{2} = '../../outputs/dwi_x_xyz1_restricted_diffusion_model.dat';
file_T{3} = '../../outputs/dwi_y_xyz1_restricted_diffusion_model.dat';
file_T{4} = '../../outputs/dwi_z_xyz1_restricted_diffusion_model.dat';


tissue_type = 1; % 0 = extra-vascular, 1 = intra-vascular, [0,1] = combined
[~, mask, fov] = read_fieldmap(file_fieldmap);

s = zeros(numel(file_m1), 1);
for i=1:numel(file_m1)
    [m1, dims, scales] = read_spinwalk(file_m1{i});
    [xyz1, ~, ~] = read_spinwalk(file_T{i}); 
    mxyz_f = filter_spinwalk(m1(:,end,:,1), xyz1(:,end,:,1), tissue_type, mask, fov);
    s(i)  = abs(complex(sum(mxyz_f(1,:)), sum(mxyz_f(2,:))));
end

s

%%
clear
clc
file_fieldmap = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/field_maps/restricted_diffusion_model.dat';
file_m1 = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/dwi_base_m1_restricted_diffusion_model.dat';
file_T = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/dwi_base_T_restricted_diffusion_model.dat';

tissue_type = 1; % 0 = extra-vascular, 1 = intra-vascular, [0,1] = combined
[~, mask, fov] = read_fieldmap(file_fieldmap);

[m1, dims, scales] = read_spinwalk(file_m1);
[T, ~, ~] = read_spinwalk(file_T); 
% mxyz_f = filter_spinwalk(m1(:,end,:,1), xyz1(:,end,:,1), tissue_type, mask, fov);
% mxyz_f = mxyz_f';






















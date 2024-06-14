clc
clear

fname{1}{1} = '../../outputs/gre_fieldmap_0.h5';
fname{1}{2} = '../../outputs/gre_fieldmap_1.h5';
fname{2}{1} = '../../outputs/se_fieldmap_0.h5';
fname{2}{2} = '../../outputs/se_fieldmap_1.h5';
fname{3}{1} = '../../outputs/ssfp_fieldmap_0.h5';
fname{3}{2} = '../../outputs/ssfp_fieldmap_1.h5';

dim_xyz  = 1;
dim_echo = 2;
dim_spin = 3;
dim_vessel_size = 4;
rad_ref_um = 53.367;

tissue_type = 0; % 0 = extra-vascular, 1 = intra-vascular, [0,1] = combined

signal_magnitude = cell(numel(fname), numel(fname{1}));
for seq = 1:numel(fname)   
    for i=1:numel(fname{seq})
        m1 = h5read(fname{seq}{i}, '/M');
        scales = h5read(fname{seq}{i}, '/scales');
        T = h5read(fname{seq}{i}, '/T');  

        m1_t = zeros(numel(scales), 1);
        for s=1:numel(scales)
            ind = ismember(T(1,end,:,s), tissue_type);
            m1_f = m1(:,end-1,ind(:),s);
            m1_t(s) = abs(complex(sum(m1_f(1,:)), sum(m1_f(2,:))));
        end
        signal_magnitude{seq, i} = m1_t;
    end
end

clf
vessel_radius  = rad_ref_um * scales;
for seq = 1:numel(fname)
    relative_signal  = 100 * (1 - signal_magnitude{seq, 1} ./ signal_magnitude{seq, 2});    
    h = semilogx(vessel_radius, relative_signal); xlabel('Vessel radius (um)'); ylabel('BOLD Signal %'); 
    hold on;
    ylim([0, 7])
end
legend('GRE', 'SE', 'SSFP')


%% grase sequence

% fname{4}{1} = '../../outputs/grase_m1_0.dat';
% fname{4}{2} = '../../outputs/grase_m1_1.dat';
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
fname{1} = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/ste_fieldmap_0.h5';
fname{2} = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/ste_fieldmap_1.h5';

dim_echo = 2;
dim_spin = 3;
dim_vessel_size = 4;
rad_ref_um = 53.367;


spins_xyz = [];
for i=1:numel(fname)
    m_xyz   = h5read(fname{i}, '/M');
    scales  = h5read(fname{i}, '/scales');
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

%% linear gradient calculations
clc
fov = [0.0e-3, 0.0e-3, 0.0e-3
       0.9e-3, 0.9e-3, 0.9e-3];
phase_range = 2*pi; % gradient will create this phase range in FoV
l = diff(fov);
time_step = 50e-6;
gamma = 267515315; % rad/s.T
if range(l) ~= 0
    error('fov must be isotropic');
end

G = phase_range / (gamma * l(1) * time_step); % T/m
disp(['Gradient strength = ' num2str(G * 1000) ' mT/m'])
% generate initial XYZ0 positions and save to h5
x = write_xyz0(fov, 101*101*101, '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/field_maps/xyz0_allzero.h5');

%% linear gradient plot
fname = '/DATA/aaghaeifar/Nextcloud/Projects/microvascular/outputs/gradient_fieldmap_allzero.h5';
m_xyz   = h5read(fname, '/M');
m_xy    = squeeze(complex(m_xyz(1,1,:), m_xyz(2,1,:)));
m_xy    = double(m_xy);
sz      = nthroot(numel(m_xy), 3);
m_xy    = reshape(m_xy, [sz, sz, sz])  ;
phase   = angle(m_xy * exp(-i*pi)); % shift the range to [-pi pi]
close all
vin(phase) % 3D plot

%% random-walk trajectory
clear
clc;
filename = '../../outputs/gre_trajectory_fieldmap_0.h5';
xyz_all = h5read(filename, '/XYZ');
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
clc
clear

folder = '../../outputs/dwi';
base = 'dwi_base_restricted_diffusion_model_sphere.h5';

b = 100:100:5000;

m1 = h5read(fullfile(folder, base), '/M'); 
s0 = abs(complex(sum(m1(1,end-1,:)), sum(m1(2,end-1,:))));
s_dwi  = zeros(numel(b), 1);
for i=1:numel(b) 
    f = fullfile(folder, ['dwi_x_' num2str(b(i)) '_restricted_diffusion_model_sphere.h5' ]);
    m1 = h5read(f, '/M'); 
    s_dwi(i) = abs(complex(sum(m1(1,end-1,:)), sum(m1(2,end-1,:))));
end
s = [s0; s_dwi];
b = [0, b];

[xData, yData] = prepareCurveData( double(b(:)), double(s(:)) );

% Virtual b-value
b = 100:100:5000;
f = fullfile(folder, ['dwi_x_' num2str(b(1)) '_restricted_diffusion_model_sphere.h5' ]);
m1 = h5read(f, '/M'); 
xyz = h5read(f, '/XYZ'); 
[theta, rho] = cart2pol(m1(1,end-1,:), m1(2,end-1,:));
G_scale = sqrt(b / b(1));
for i=1:numel(b)     
    % sig = squeeze(rho) .* exp(1j*(squeeze(theta) * G_scale(i)));
    % s_dwi(i) = abs(sum(sig));
    [x, y] = pol2cart(squeeze(theta) * G_scale(i), squeeze(rho));
    s_dwi(i) = abs(complex(sum(x), sum(y)));
end
s2 = [s0; s_dwi];
b2 = [0, b];5

[xData_v, yData_v] = prepareCurveData( double(b2(:)), double(s2(:)) );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.DiffMinChange = 1e-09;
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.StartPoint = [988912.916299669 -0.000461964926600495];
opts.TolFun = 1e-09;
opts.TolX = 1e-09;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
[fitresult_v, gof_v] = fit( xData_v, yData_v, ft, opts );

% Plot fit with data.
plot( fitresult, xData, yData); hold on;
plot(fitresult_v, xData_v, yData_v); hold off;
legend('S vs. b', 'dwi', 'S_v vs. b', 'dwi_v','Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'b', 'Interpreter', 'none' );
ylabel( 'S', 'Interpreter', 'none' );
grid on
disp(['Diff. Const = ' num2str(abs(fitresult.b) * 1e-6) ' m/s']);
disp(['Diff. Const v = ' num2str(abs(fitresult_v.b) * 1e-6) ' m/s']);


%%




















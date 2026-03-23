function [ktraj_meas, ktraj_hash] = TRAJ_reco(path_raw, path_backup, vendor, plot_id)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026
% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V2, 22.03.2026; unify for different vendors

% ----- Input -----
% path_raw:    path of meas data or [] for select via uigetfile
% path_backup: path of pulseq workspace backup; not necessary for Siemens
% vendor:      vendor name

if strcmp(vendor, 'UI') % trajectory data are inverted for United Imaging
    vz = -1;
else
    vz = 1;
end

%% load trajectory rawdata and pulseq backup
[rawdata, ~, PULSEQ] = pulseq_read_meas(path_raw, path_backup, vendor);

%% read dimensions
NCoils = size(rawdata,1);  % number of coils
NR     = PULSEQ.TRAJ.NR;   % number of spiral arms
NRead  = size(rawdata, 3); % number of sampling points
Nav    = PULSEQ.TRAJ.Nav;  % number of averages

%% Duyn method: 10.1006/jmre.1998.1396
if strcmp(PULSEQ.TRAJ.method, 'duyn')

    % average traj data
    rawdata = squeeze(mean(reshape(rawdata, NCoils, Nav, [], NRead), 2));
    
    % PCA coil compression: x slice
    rawdata_x = rawdata(:,1:size(rawdata,2)/2,:);
    rawdata_x = rawdata_x(:,:);
    [~, ~, pc]  = svd(rawdata_x.', 'econ');
    rawdata_x = pc.' * rawdata_x(:,:);
    rawdata_x = reshape(rawdata_x, [NCoils NR*2 NRead] );
    rawdata_x = squeeze(rawdata_x(1,:,:));

    % PCA coil compression: y slice
    rawdata_y = rawdata(:,1+size(rawdata,2)/2:end,:);
    rawdata_y = rawdata_y(:,:);
    [~, ~, pc]  = svd(rawdata_y.', 'econ');
    rawdata_y = pc.' * rawdata_y(:,:);
    rawdata_y = reshape(rawdata_y, [NCoils NR*2 NRead] );
    rawdata_y = squeeze(rawdata_y(1,:,:));
    clear pc rawdata;

    % correct with reference scans
    temp_ind  = 1:2:2*NR;
    rawdata_x = squeeze(rawdata_x(temp_ind, :)) ./ squeeze(rawdata_x(temp_ind+1, :));
    rawdata_y = squeeze(rawdata_y(temp_ind, :)) ./ squeeze(rawdata_y(temp_ind+1, :));
    clear temp_ind rawdata;
    
    % unwrap signal phase
    phase_x = rawdata_x *0;
    phase_y = rawdata_y *0;    
    for j=1:NR
        phase_x(j,:) = unwrap(angle(squeeze(rawdata_x(j,:))));
        phase_y(j,:) = unwrap(angle(squeeze(rawdata_y(j,:))));
    end
    clear rawdata_x rawdata_y;
    
    % get k-space trajectory via slice offset
    ktraj_meas_x = phase_x / PULSEQ.TRAJ.slice_offset /2/pi;
    ktraj_meas_y = phase_y / PULSEQ.TRAJ.slice_offset /2/pi;
    ktraj_meas(1,:,:) = ktraj_meas_x;
    ktraj_meas(2,:,:) = ktraj_meas_y;
    clear phase_x phase_y;

end

%% Robison method: 10.1002/mrm.27583
if strcmp(PULSEQ.TRAJ.method, 'robison')

    % split data
    temp_ind_x1 = 1:4:size(rawdata,2)/2;
    temp_ind_x2 = 3:4:size(rawdata,2)/2;
    temp_ind_x3 = 2:4:size(rawdata,2)/2;
    temp_ind_x4 = 4:4:size(rawdata,2)/2;
    temp_ind_y1 = temp_ind_x1 + size(rawdata,2)/2;
    temp_ind_y2 = temp_ind_x2 + size(rawdata,2)/2;
    temp_ind_y3 = temp_ind_x3 + size(rawdata,2)/2;
    temp_ind_y4 = temp_ind_x4 + size(rawdata,2)/2;
    rawdata_x1  = rawdata(:,temp_ind_x1,:);
    rawdata_x2  = rawdata(:,temp_ind_x2,:);
    rawdata_x3  = rawdata(:,temp_ind_x3,:);
    rawdata_x4  = rawdata(:,temp_ind_x4,:);
    rawdata_y1  = rawdata(:,temp_ind_y1,:);
    rawdata_y2  = rawdata(:,temp_ind_y2,:);
    rawdata_y3  = rawdata(:,temp_ind_y3,:);
    rawdata_y4  = rawdata(:,temp_ind_y4,:);
    clear rawdata temp_ind_x1 temp_ind_x2 temp_ind_x3 temp_ind_x4 temp_ind_y1 temp_ind_y2 temp_ind_y3 temp_ind_y4;

    % average data
    rawdata_x1 = squeeze(mean(reshape(rawdata_x1, NCoils, Nav, [], NRead), 2));
    rawdata_x2 = squeeze(mean(reshape(rawdata_x2, NCoils, Nav, [], NRead), 2));
    rawdata_x3 = squeeze(mean(reshape(rawdata_x3, NCoils, Nav, [], NRead), 2));
    rawdata_x4 = squeeze(mean(reshape(rawdata_x4, NCoils, Nav, [], NRead), 2));
    rawdata_y1 = squeeze(mean(reshape(rawdata_y1, NCoils, Nav, [], NRead), 2));
    rawdata_y2 = squeeze(mean(reshape(rawdata_y2, NCoils, Nav, [], NRead), 2));
    rawdata_y3 = squeeze(mean(reshape(rawdata_y3, NCoils, Nav, [], NRead), 2));
    rawdata_y4 = squeeze(mean(reshape(rawdata_y4, NCoils, Nav, [], NRead), 2));

    % PCA coil compression: x slice
    rawdata_x1 = rawdata_x1(:,:);
    rawdata_x2 = rawdata_x2(:,:);
    rawdata_x3 = rawdata_x3(:,:);
    rawdata_x4 = rawdata_x4(:,:);
    [~, ~, pc] = svd([rawdata_x1, rawdata_x2, rawdata_x3, rawdata_x3].', 'econ');
    rawdata_x1 = pc.' * rawdata_x1(:,:);
    rawdata_x2 = pc.' * rawdata_x2(:,:);
    rawdata_x3 = pc.' * rawdata_x3(:,:);
    rawdata_x4 = pc.' * rawdata_x4(:,:);
    rawdata_x1 = reshape(rawdata_x1, [NCoils NR NRead] );
    rawdata_x2 = reshape(rawdata_x2, [NCoils NR NRead] );
    rawdata_x3 = reshape(rawdata_x3, [NCoils NR NRead] );
    rawdata_x4 = reshape(rawdata_x4, [NCoils NR NRead] );
    rawdata_x1 = squeeze(rawdata_x1(1,:,:));
    rawdata_x2 = squeeze(rawdata_x2(1,:,:));
    rawdata_x3 = squeeze(rawdata_x3(1,:,:));
    rawdata_x4 = squeeze(rawdata_x4(1,:,:));
    clear pc;

    % PCA coil compression: y slice
    rawdata_y1 = rawdata_y1(:,:);
    rawdata_y2 = rawdata_y2(:,:);
    rawdata_y3 = rawdata_y3(:,:);
    rawdata_y4 = rawdata_y4(:,:);
    [~, ~, pc] = svd([rawdata_y1, rawdata_y2, rawdata_y3, rawdata_y3].', 'econ');
    rawdata_y1 = pc.' * rawdata_y1(:,:);
    rawdata_y2 = pc.' * rawdata_y2(:,:);
    rawdata_y3 = pc.' * rawdata_y3(:,:);
    rawdata_y4 = pc.' * rawdata_y4(:,:);
    rawdata_y1 = reshape(rawdata_y1, [NCoils NR NRead] );
    rawdata_y2 = reshape(rawdata_y2, [NCoils NR NRead] );
    rawdata_y3 = reshape(rawdata_y3, [NCoils NR NRead] );
    rawdata_y4 = reshape(rawdata_y4, [NCoils NR NRead] );
    rawdata_y1 = squeeze(rawdata_y1(1,:,:));
    rawdata_y2 = squeeze(rawdata_y2(1,:,:));
    rawdata_y3 = squeeze(rawdata_y3(1,:,:));
    rawdata_y4 = squeeze(rawdata_y4(1,:,:));
    clear pc;

    % correct phases accoring to doi: 10.1002/mrm.27583
    rawdata_xA = rawdata_x1 ./ rawdata_x2;
    rawdata_xB = rawdata_x4 ./ rawdata_x3;
    rawdata_yA = rawdata_y1 ./ rawdata_y2;
    rawdata_yB = rawdata_y4 ./ rawdata_y3;
    clear rawdata_x1 rawdata_x2 rawdata_x3 rawdata_x4 rawdata_y1 rawdata_y2 rawdata_y3 rawdata_y4

    % unwrap signal phases
    if NR==1
        rawdata_xA = rawdata_xA.';
        rawdata_xB = rawdata_xB.';
        rawdata_yA = rawdata_yA.';
        rawdata_yB = rawdata_yB.';
    end
    phase_xA = zeros(size(rawdata_xA));
    phase_xB = zeros(size(rawdata_xB));
    phase_yA = zeros(size(rawdata_yA));
    phase_yB = zeros(size(rawdata_yB));
    for j=1:NR
        phase_xA(j,:) = unwrap(angle(squeeze(rawdata_xA(j,:))));
        phase_xB(j,:) = unwrap(angle(squeeze(rawdata_xB(j,:))));
        phase_yA(j,:) = unwrap(angle(squeeze(rawdata_yA(j,:))));
        phase_yB(j,:) = unwrap(angle(squeeze(rawdata_yB(j,:))));
    end
    clear rawdata_xA rawdata_xB rawdata_yA rawdata_yB;

    % get k-space trajectory via slice offset
    ktraj_meas_x = vz * (phase_xA + phase_xB) / 4 / PULSEQ.TRAJ.slice_offset /2/pi;
    ktraj_meas_y = vz * (phase_yA + phase_yB) / 4 / PULSEQ.TRAJ.slice_offset /2/pi;
    ktraj_meas(1,:,:) = ktraj_meas_x;
    ktraj_meas(2,:,:) = ktraj_meas_y;
    
    % get eddy current phase
    eddy_x = (phase_xA - phase_xB) / 4;
    eddy_y = (phase_yA - phase_yB) / 4;
    
end

%% load calculated trajectory -> generate hash ID
ktraj_calc   = SPI_load_ktraj(PULSEQ.PULSEQ_SPI);
ktraj_hash   = pulseq_get_wave_hash(ktraj_calc(:));
ktraj_calc_x = squeeze(ktraj_calc(1,:,:));
ktraj_calc_y = squeeze(ktraj_calc(2,:,:));
if NR==1
    ktraj_calc_x = ktraj_calc_x.';
    ktraj_calc_y = ktraj_calc_y.';
end

%% plot trajectory

if sum(abs(plot_id))>0

if NR==1
    plot_id = 1;
end

% get time axis
Nread_meas = NRead;
Nread_calc = size(ktraj_calc_x,2); % if rewinder was acquired: Nread_meas > Nread_calc
t_meas = (1:Nread_meas) *PULSEQ.TRAJ.adc.dwell *1e3;  % [ms] time with grad raster steps
t_calc = (1:Nread_calc) *PULSEQ.PULSEQ_SPI.SPI.adc.dwell  *1e3;  % [ms] time with grad raster steps
kmax = max(abs([ktraj_calc_x(:); ktraj_calc_x(:)]));  % kspace limit

% compare trajectories: kx(t) & ky(t)
figure()
tiledlayout(2,1)
ax1 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_calc, ktraj_calc_x(plot_id(j),:), 'b-', 'LineWidth', 5)
    plot(t_meas, ktraj_meas_x(plot_id(j),:), 'r.', 'MarkerSize', 5)
end
xlabel('time t [ms]')
ylabel('kx [1/m]')
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

ax2 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_calc, ktraj_calc_y(plot_id(j),:), 'b-', 'LineWidth', 5)
    plot(t_meas, ktraj_meas_y(plot_id(j),:), 'r.', 'MarkerSize', 5)
end
xlabel('time t [ms]')
ylabel('ky [1/m]')
legend('calc', 'meas')
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

linkaxes([ax1 ax2], 'x')

% compare trajectories: error_x(t) & error_y(t)
figure()
tiledlayout(2,1)
ax1 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_calc, (ktraj_meas_x(plot_id(j),1:Nread_calc)-ktraj_calc_x(plot_id(j),:))/kmax*100, '-', 'LineWidth', 5)
end
xlabel('time t [ms]')
ylabel('error [%]')
ylim([-1 1])
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

ax2 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_calc, (ktraj_meas_y(plot_id(j),1:Nread_calc)-ktraj_calc_y(plot_id(j),:))/kmax*100, '-', 'LineWidth', 5)
end
xlabel('time t [ms]')
ylabel('error [%]')
ylim([-1 1])
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

linkaxes([ax1 ax2], 'x')

% compare trajectories: ky(kx)
figure()
hold on
for j=1:numel(plot_id)
    plot(ktraj_calc_x(plot_id(j),:), ktraj_calc_y(plot_id(j),:), 'b-', 'LineWidth',  5)
    plot(ktraj_meas_x(plot_id(j),:), ktraj_meas_y(plot_id(j),:), 'r.', 'MarkerSize', 5)
end
xlabel('kx [1/m]')
ylabel('ky [1/m]')
hold off
axis image;
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

% eddy currents
if strcmp(PULSEQ.TRAJ.method, 'robison')
figure()
tiledlayout(2,1)
ax1 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_meas, eddy_x(plot_id(j),:)/pi*180, 'r.', 'MarkerSize', 5)
end
xlabel('time t [ms]')
ylabel('eddy current x [deg]')
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

ax2 = nexttile;
hold on
for j=1:numel(plot_id)
    plot(t_meas, eddy_y(plot_id(j),:)/pi*180, 'r.', 'MarkerSize', 5)
end
xlabel('time t [ms]')
ylabel('eddy current y [deg]')
hold off
set(gca,'linewidth', 2, 'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

linkaxes([ax1 ax2], 'x')
end

end

end
% Script to analyse the stress-strain data as obtained by Sander Spiekhout
% using the acoustical camera (Erasmus University, Rotterdam) to obtain
% bubble P.shell parameters: P.shell elasticity (P.chi = A dsigma/dA) from the
% elastic contribution, and P.shell viscosity (kappa_s) from the viscous
% contribution

% University of Twente,
% Charlotte Nawijn, 2021 - 2024

% change when analyzing each data set:
% - data path to load (below)

%% General
clear, clc, close all
warning off

% If using GPU server to do calculations, the plotting will be turned off
path_cd = cd;
if path_cd(1) == 'C' || path_cd(1) == 'D'
    S.calculate_on_server = 0;
else
    S.calculate_on_server = 1;
end

%% Load Sander's data from the acoustical camera setup --------------------
% data is selected according to SNR in dB (SNRVal): >25 dB
% SNR is determined by comparing region in raw signal [12e3:76e3] with [20e4:40e4]

addpath('../Data and experiment prep')

% % vial C (2.5 µm radius, 60 deg C)
% data_name = 'StrainSignals_VialC_35kPa_1.8MHz_Mono_2.5um_16082023_V2';
% load(['../Data and experiment prep/August 2023 experiments/'...
%     'vial C_1.8 MHz_35 kPa/'...
%     data_name '.mat'])

% vial B (2.3 µm radius, 60 deg C)
data_name = 'StrainSignals_VialB_40kPa_1.9MHz_Mono_2.3um_17082023_V2';
load(['../Data and experiment prep/August 2023 experiments/'...
    'vial B_1.9 MHz_40 kPa/'...
    data_name '.mat'])

data_name = erase(data_name, '_V2');
vial_name_temp = extractBetween(data_name, "Vial", "_");
vial_name = ['vial' vial_name_temp{1}];

% rename the data
strain_full = relWfms;   %strainY(:,1:51801)

% Parametrization Tim Segers
addpath('../Stress_strain solver simulation')
F.fit_surf_tens = load('../Stress_strain solver simulation/Parametrisatie oppervlaktespanning (Tim)/fit_SigmaR_04-08-2017.mat');

% other functions
addpath('functions')

% Phase of the low-frequency transducer, measured using fibre-optic hydrophone
addpath('../Fibre-optic hydrophone measurements 20231020')
if contains(data_name, '1.9')
    phase_LF = load("20231020_1d9MHz_1000rep_PRF10Hz.mat");
elseif contains(data_name, '1.8')
    phase_LF = load("20231020_1d8MHz_1000rep_PRF10Hz_2.mat");
end

input_param.phase_LF = phase_LF.phase;

% folders where to save figures
saving_folder_paper = [cd '\..\..\..\Artikelen\'...
    'Stress-strain\figures\Results\raw plot\' erase(data_name, 'StrainSignals_') ...
    '\viscous contribution\'];

saving_folder_figs = [cd '/viscous contribution/figures '...
    erase(data_name, 'StrainSignals_') '/'];


%% Initialize
% selIdx = find(strainVal>0.08 & SNR'>=25 & strainVal < 0.2);% & R0est>2.2e-6 & R0est<2.4e-6);% & strainVal <0.2);  %hannStrainVal is the average strain value over the entire 200 µs LF pulse (dus niet op de piek!)        %abs(1-(refLvl2./refLvl1))<=0.3) is already in the data (abs(1-(refLvl2./refLvl1)) = FluctVal?)
selIdx = find(strainVal>0.1 & SNR'>=25 & strainVal < 0.2);% & R0est>2.2e-6 & R0est<2.4e-6);% & strainVal <0.2);  %hannStrainVal is the average strain value over the entire 200 µs LF pulse (dus niet op de piek!)        %abs(1-(refLvl2./refLvl1))<=0.3) is already in the data (abs(1-(refLvl2./refLvl1)) = FluctVal?)

% array of which bubble(s) to analyse
wfmIndx_array = selIdx(1:end);

%% =======settings (S)===================================================================
S.data_name = data_name;

% booleans for plotting, saving, etc.
if S.calculate_on_server
    S.save_data = 1;
    S.saving_figs = 0;
    S.do_error_minimization = 1;        % to plot the surface tension figures of the results on the server, disable this
    
    S.plot_all = 0;
else
    S.save_data = 1;
    S.saving_figs = 1;
    S.do_error_minimization = 0;        % to plot the surface tension figures of the results on the server, disable this

    S.plot_all = 1;

    S.analyse_viscous = 1;
    S.saving_figs_visc = 1;
end

% How to determine the elastic contribution. 1= from the line at zero
% strain rate of the full visocelatic plane, 0 = from the points at the envelope of 
% the strain filled into the expression at zero strain rate
S.elastic_from_full_viscoelastic = 1;

% 1= filtering the strain around the fundamental and 2nd harmonic. 0 = no filtering
S.filtering = 1;

% 1= analyse full ramp-up and -down in pressure amplitude,
% 2= only ramp-up,
% 3= only ramp-down
S.section_analysis = 1;

% width of the mask used for filtering in Fourier spectrum of strain around
% the fundamental and second harmonic frequencies
S.width_mask = 0.6e6;       % Hz

%number of points to append (mirrored) to append before filtering and
% obtaining derivatives of the strain to get rid of edge effects
S.ApLg = 500;

% how much to cut off of f array before interpolation (w/ griddata) to
% prevent upshoot due to interpolation of few points at strain extremes
S.mask_cutoff = 0.1;  %fraction

%range in strain on both sides over which the error minimization is performed
S.range_EM = 0.04; 

% step size in strain to determine the surface tension curve from the
% parametrization from Segers
S.strain_input_step = 1e-4;

% number of samples over which to do a moving average to smooth the surface
% tension for the fitting
S.smoothingN = round(S.range_EM/S.strain_input_step/4);

% how many times to run the lsqfit algorithm per bubble (attempts)
S.attempts = 100;


%% ======input parameters (P)===================================================================
% general parameters
P.cw = 1500;                              % speed of sound in water (m/s)
P.P0 = 1e5;                               % ambient pressure (P.Pa)
P.mu = 2*1e-3;                            % viscosity of water (P.Pa*s) multiplied by two to account for thermal damping %1e-3;
P.rho = 1e3;                              % density of water (kg/m^3)
P.sig = 0.072;                            % surface tension water (N/m)

% bubble variables
P.dR0 = 0;                         % initial velocity bubble wall (m/s)
P.kap = 1.07;                      % polytropic exponent C3F8 gas
P.chi = 0.5;      % Tim: 0.55      % P.shell stiffness (N/m)
P.ks = 1e-8;                       % P.shell viscosity (P.Pa*s)
P.shell = 1;                       % boolean whether the bubble has a P.shell
if P.shell == 1
    P.sig_0 = 30e-3; %35e-3;  % Tim: 15e-3   % initial surface tension of the P.shell (N/m)
else
    P.sig_0 = P.sig;
end

% ultrasound variables
freq_string = extractBetween(data_name,"kPa_","MHz");
freq_string = replace(freq_string, 'd', '.');
P.fUS = str2double(freq_string{1})*1e6;       % frequency of the transmitted ultrasound (Hz)
pressure_string = extractBetween(data_name, ['Vial' vial_name_temp{1} '_'], 'kPa');
P.Pa = str2double(pressure_string{1})*1e3;                        % acoustic pressure amplitude (P.Pa); 60 kPa for 2.0 µm data set, 30 kPa for 2.8 µm, 25 kPa for 3.3 µm
P.T = 1/P.fUS;                        % period of transmit frequency (sec)
P.Fs = 250e6;                       % sampling frequency of digitizer (Hz)
P.dt = 1/P.Fs;                        % sampling time step (sec)

init_rad_string = extractBetween(data_name,"Mono_","um");
P.R0 = str2double(replace(init_rad_string, 'd', '.')); % initial radius in µm


%% =======driving ultrasound parameters (D)===================================================================
D.time_temp = 0:P.dt:200e-6;                  % time vector (sec)
% P.T = P.T*1e6;                        % from sec to µs

% HF pulse is 200 + 2*37.5 µs, strain is cropped from sample 500 to 64750
% (2 µs after start HF to 16 µs before end HF).
% see 'SS_Processing_cropping_strain_CN_20231127.m' to see how I came to the offsets of D.time_full 
D.time_full = -35.08e-6:P.dt:((size(D.time_temp,2)-1)*P.dt + 21.92e-6);         % vial B (40 kPa, 1.9 MHz, 2.3 µm radius)

% make envelope (pad with zeros to match the time of the strain: +- 4 µs)
D.env_hanning_temp = hann(round(length(D.time_temp)))';
D.envelop_hanning_full = [zeros(1, round(abs(D.time_full(1))/P.dt))...
    D.env_hanning_temp...
    zeros(1, length(D.time_temp) - length(D.env_hanning_temp))...
    zeros(1, round(abs(D.time_full(end)-200e-6)/P.dt))];

D.Pacc_full = P.Pa*D.envelop_hanning_full.*sin(2*pi*P.fUS.*D.time_full);    % incident acoustic pressure (P.Pa)

% set P.T from -5 to 105 µs for the cropped signals (only ramp-up)
D.time_up = -5e-6: P.dt: 105e-6;

% crop the strain and driving pressure to the correct range
if S.section_analysis == 1
    [~, crop_idx_start_driving] = min(abs(D.time_full - 0));
    [~ , crop_idx_end_driving] = min(abs(D.time_full - 200e-6));
elseif S.section_analysis == 2
    [~, crop_idx_start_driving] = min(abs(D.time_full - 0));
    [~ , crop_idx_end_driving] = min(abs(D.time_full - 100e-6));
elseif S.section_analysis == 3
    [~, crop_idx_start_driving] = min(abs(D.time_full - 100e-6));
    [~ , crop_idx_end_driving] = min(abs(D.time_full - 200e-6));
end

D.envelop_hanning = D.envelop_hanning_full(crop_idx_start_driving:crop_idx_end_driving);
D.time = D.time_full(crop_idx_start_driving:crop_idx_end_driving);
D.Pacc = D.Pacc_full(crop_idx_start_driving:crop_idx_end_driving);


%% ====== Time for strain, see 'SS_Processing_cropping_strain_CN_20231127.m' to see how 
% I came to the offsets of D.time_full (based on cross-correlation
% simulated strain and measured strains

% crop the data
time_strain = D.time_full;%-35.5e-6:P.dt:((size(D.time_temp,2)-1)*P.dt + 21.5e-6);        % time vector for strain

if S.section_analysis == 1
    [~, crop_idx_start] = min(abs(time_strain - 0));
    [~ , crop_idx_end] = min(abs(time_strain - 200e-6));
elseif S.section_analysis == 2
    [~, crop_idx_start] = min(abs(time_strain - 0));
    [~ , crop_idx_end] = min(abs(time_strain - 100e-6));
elseif S.section_analysis == 3
    [~, crop_idx_start] = min(abs(time_strain - 100e-6));
    [~ , crop_idx_end] = min(abs(time_strain - 200e-6));
end

%% Analyse every bubble one by one ----------------------------------------
warning('off')

for wfmIndx = wfmIndx_array(1:end)
    close all

    clc
    fprintf('Storing variables: %3.f%%\n', wfmIndx/wfmIndx_array(end)*100)
    if ismember(wfmIndx, wfmIndx_array)

        saving_folder = [cd '/error minimization data/'];

        if S.section_analysis == 1
            saving_folder = [saving_folder 'updown/'];
        elseif S.section_analysis == 2
            saving_folder = [saving_folder 'up/'];
        elseif S.section_analysis ==3
            saving_folder = [saving_folder 'down/'];
        end

        data_name = erase(data_name, '_V2');

        saving_folder = [saving_folder...
            'range_EM_' num2str(S.range_EM) '_mask_cutoff_'...
            num2str(S.mask_cutoff) '/' data_name '/'];
        saving_dir = [saving_folder ...
            'bubble ' num2str(wfmIndx)];
        saving_file_name = [saving_dir '.mat'];

        if S.calculate_on_server || (~S.calculate_on_server && S.do_error_minimization)
            if exist(saving_file_name, 'file') == 2
                continue
            end
        elseif ~S.calculate_on_server && ~S.do_error_minimization         
            if ~(exist(saving_file_name, 'file') == 2)
                continue
            end
        end

       
        if exist(saving_file_name, 'file')
            S_temp = S;

            load(saving_file_name)

            S_temp.lb = S.lb;
            S_temp.ub = S.ub;
            S=S_temp;
            clear S_temp   % to not overwrite the settings
        else
            continue
        end

        %% Store data to use in the viscous contribution analysis
        viscoelastic_contrib_ND_all(wfmIndx, :, :) = viscoelastic_contrib_int;
        elastic_contrib_ND_all(wfmIndx, :) = elastic_contrib;

        strain_all(wfmIndx, :) = strain_vec;
        dstrain_all(wfmIndx, :) = dstrain_vec;

        strain_grid_all(wfmIndx, :, :) = strain_grid;
        dstrain_grid_all(wfmIndx, :, :) = dstrain_grid;

    end
end


%% Based on all strain and strain rate 'vectors', make new arrays based on
% the minimum and maximum values of all bubbles
strain_vec_all = linspace(min(strain_all, [], 'all'), max(strain_all, [], 'all'), 1e3);
dstrain_vec_all = linspace(min(dstrain_all, [], 'all'), max(dstrain_all, [], 'all'), 1e3);

[strain_grid, dstrain_grid] = meshgrid(strain_vec_all, dstrain_vec_all);

axes_limits_strain = [-max(abs(strain_vec_all)), max(abs(strain_vec_all))];
axes_limits_dstrain = [-max(abs(dstrain_vec_all)), max(abs(dstrain_vec_all))];

axes_limits_strain = axes_limits_strain*0.8;
axes_limits_dstrain = axes_limits_dstrain*0.8;


%% Viscous contribution
% Per bubble, determine the viscous contribution and interpolate the
% viscous contribution to the new vectors. Then make a binary mask of the non-NaN values

% initialize the mask of the viscous contribution
if S.analyse_viscous
    visc_contrib_mask = NaN(wfmIndx_array(end), length(strain_vec_all), length(dstrain_vec_all));
    visc_contrib_ND_all = NaN(wfmIndx_array(end), length(strain_vec_all), length(dstrain_vec_all));
else
    viscelas_contrib_ND_all = NaN(wfmIndx_array(end), length(strain_vec_all), length(dstrain_vec_all));
    elas_contrib_ND_all = NaN(wfmIndx_array(end), length(strain_vec_all), length(dstrain_vec_all));
end

for wfmIndx = wfmIndx_array(1:end)

    clc
    fprintf('Determining viscoelastic, elastic and viscous contribution: %3.f%%\n', wfmIndx/wfmIndx_array(end)*100)


    S.plot_all = 0;

    if S.analyse_viscous
        % data
        visc_contrib_temp = squeeze(zeros(size(viscoelastic_contrib_ND_all(wfmIndx,:,:))));  % stored for (strain rate, strain)
        for i = 1: size(viscoelastic_contrib_ND_all, 3)        % for all strain
            for j = 1: size(viscoelastic_contrib_ND_all, 2)       % for all strain rate
                visc_contrib_temp(j,i) = squeeze(viscoelastic_contrib_ND_all(wfmIndx, j, i)) -...
                    elastic_contrib_ND_all(wfmIndx, i);
            end
        end

        % interpolate (linearly) data onto strain and strain rate values based on
        % all bubbles
        visc_contrib_ND_all(wfmIndx, :, :) = griddata(squeeze(strain_grid_all(wfmIndx, :, :)),...
            squeeze(dstrain_grid_all(wfmIndx, :, :)), visc_contrib_temp, ...
            strain_grid, dstrain_grid, 'linear');

        % smoothing data by moving mean
        visc_contrib_ND_filtered = smoothdata2(squeeze(visc_contrib_ND_all(wfmIndx,:, :)),...
            'movmean', 20);

        % make a binary mask per bubble: 0 = NaN values, 1 = non-NaN values
        nonNaN_idx = ~isnan(visc_contrib_ND_filtered);
        visc_contrib_mask(wfmIndx, :, :) = nonNaN_idx;


    else %mean viscoelastic and elastic planes
        viscelas_contrib_ND_all(wfmIndx, :, :) = griddata(squeeze(strain_grid_all(wfmIndx, :, :)),...
            squeeze(dstrain_grid_all(wfmIndx, :, :)), ...
            squeeze(viscoelastic_contrib_ND_all(wfmIndx, :, :)), ...
            strain_grid, dstrain_grid, 'linear');

        elastic_contrib_ND_all_temp = repmat(elastic_contrib_ND_all(wfmIndx, :), 1, 1, length(dstrain_grid_all(wfmIndx, :, :)));
        % stored as: wfmIndx, strain rate, strain;  (strain rate varies
        % over the columns, strain over the rows)
        elastic_contrib_ND_all_temp = permute(elastic_contrib_ND_all_temp, [1 3 2]);
        
        elas_contrib_ND_all(wfmIndx, :, :) = griddata(squeeze(strain_grid_all(wfmIndx, :, :)),...
            squeeze(dstrain_grid_all(wfmIndx, :, :)), ...
            squeeze(elastic_contrib_ND_all_temp), ...
            strain_grid, dstrain_grid, 'linear');
    end

end


%% Based on masks of individual bubbles, determine common mask (mean)
% (with threshold to pick which values of strain and strain rate to include)

if S.analyse_viscous
    clc
    fprintf('\nDetermining mean mask of viscous contribution\n')
    
    % not sure if this gives the right dimension
    mean_mask = squeeze(mean(visc_contrib_mask, 1, 'omitnan'));
    
    clear visc_contrib_mask
    
    figure(100)
    surf(strain_vec_all, dstrain_vec_all, mean_mask,...
        'edgecolor', 'none')
    xlabel('radial strain {\epsilon}', 'interpreter', 'tex')
    ylabel('radial strain rate {\epsilon}', 'interpreter', 'tex')
    c=colorbar();
    colormap(jet)
    clim([0 1])
    title({'binary mask'})%, ['measurement index: ' num2str(wfmIndx)]})
    grid on
    view(0,90)
    xlim(axes_limits_strain)
    ylim(axes_limits_dstrain)
    
    threshold = 0.85; %higher threshold means less of the strain and strain rate values remain, but more bubbles will be included
    
    
    % Make new binary mask based on the mean mask an the threshold 
    mask_sel_idx = mean_mask>=threshold;%double(mean_mask>=threshold);
    [mask_row_idx, mask_col_idx] = find(mean_mask >= threshold); % mean_mask(mask_row_idx, mask_col_idx) are the nonzero values
    [mask_1D_idx] = find(mean_mask >= threshold); % mean_mask(mask_row_idx, mask_col_idx) are the nonzero values
    
    
    %% Bubbles with values outside will not be considered in SVD, check how many
    % will be discared
    bubble_sel_idx = find(sum(isnan(visc_contrib_ND_all(:, mask_sel_idx)),2)==0)';
    
    % how many of the bubbles remain after filtering based on mean mask
    fprintf('Bubbles remaining (threshold %1.2f): %3.2f%%\n', threshold, ...
        100*length(bubble_sel_idx)/length(wfmIndx_array))

end

%% Singular value decomposition (SVD)
if S.analyse_viscous
% Per bubble, make the viscous contribution (2D) one column-array (1D) 
% and put the 1D viscous contributions of all bubbles next to each other
% (bubble per row, columns are 'pixels': set of strain and strain rate)
clc
fprintf('Performing singular value decomposition (SVD)\n')

% data centering, around mean value (not mean feature! = (strain, strain rate) positions)
visc_contrib_ND_mean = mean(visc_contrib_ND_all(bubble_sel_idx, :, :), 1, 'includenan');
visc_contrib_ND_median = median(visc_contrib_ND_all(bubble_sel_idx, :, :), 1, 'includenan');
visc_contrib_ND_std = std(visc_contrib_ND_all(bubble_sel_idx, :, :), 'includenan');

visc_contrib_ND_centered = visc_contrib_ND_all(bubble_sel_idx, :, :) - visc_contrib_ND_mean(:,:,:);

visc_contrib_ND_col = visc_contrib_ND_centered(:, mask_sel_idx);

addpath('..\..\..\Lab\SVD')

% perform SVD
[U, Sigma, V] = svdecon(visc_contrib_ND_col);           % output: sigma = diagonal matrix of singular values

% plot the features and observations (of the first ~10 features?)
figure(1000)
clf
plot(diag(Sigma), '.-', 'linewidth', 2, 'markersize', 16)
xlabel('singular values')
ylabel('prevalence')
grid on
set(gca, 'fontsize', 12)
xlim([0 10])


% for the first couple features, take only that feature (filter out all others)
% and go back to 'viscous contribution' (a.f.o. of strain and strain rate)
% to see the most dominant features in terms of strain and strain rate

V_trans = V';

for sin_val_idx = 1: 4
    sin_val_min = sin_val_idx;
    sin_val_max = sin_val_idx;
    
    remove_arr = [1:sin_val_min - 1, sin_val_max+1: length(Sigma)];     % values of eigenvalues to remove for both PC1 and PC2
    Sigma_filtered = Sigma;
    Sigma_filtered(remove_arr, remove_arr) = 0;

    figure(1000)
    clf
    hold on
    plot(diag(Sigma), '.-', 'linewidth', 2, 'markersize', 16)
    xlabel('singular value index')
    ylabel('prevalence')
    grid on
    plot(diag(Sigma_filtered), '.-', 'color', 'k', 'linewidth', 2, 'markersize', 16)
    set(gca, 'fontsize', 12)
    xlim([0 20])
    legend('unfiltered', 'filtered')

    drawnow

    % reconstruct original matrix
    visc_contrib_ND_filtered = U*Sigma_filtered*V';

    % assign the value to fill 3D space with 2D array
    visc_contrib_ND_filtered_3D = NaN(length(bubble_sel_idx), size(strain_grid,1), size(strain_grid,2));

    for val_idx = 1: length(visc_contrib_ND_filtered)
    
        if mod(val_idx, 100) == 0
            clc
            fprintf('Filtering singular value %d\n', sin_val_idx)
            fprintf('Reshaping: %3.2f%%\n', ...
            val_idx / length(visc_contrib_ND_filtered)* 100)
        elseif val_idx == length(visc_contrib_ND_filtered)
            clc
            fprintf('Filtering singular value %d\n', sin_val_idx)
            fprintf('Reshaping: %3.2f%%\n', ...
            val_idx / length(visc_contrib_ND_filtered)* 100)
        end

        % find corresponding place in 2D (strain and strain rate) space
        loc_idx = mask_1D_idx(val_idx);

        visc_contrib_ND_filtered_3D(:, loc_idx) =...
            visc_contrib_ND_filtered(:, val_idx);

    end

end
end


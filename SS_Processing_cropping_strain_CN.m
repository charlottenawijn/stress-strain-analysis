% Script to analyse the stress-strain data as obtained by Sander Spiekhout
% using the acoustical camera (Erasmus University, Rotterdam) to obtain
% bubble P.shell parameters: P.shell elasticity (P.chi = A dsigma/dA) from the
% elastic contribution, and P.shell viscosity (kappa_s) from the viscous
% contribution

% University of Twente,
% Charlotte Nawijn, 2021 - 2023

% change when analyzing each data set:
% - data path to load (below)

%% General
clear, clc, close all

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

vial_name = extractBetween(data_name, "Vial", "_");

% rename the data
strain_full = relWfms;  

% Parametrization Tim Segers
addpath('../Stress_strain solver simulation')
F.fit_surf_tens = load('../Stress_strain solver simulation/Parametrisatie oppervlaktespanning (Tim)/fit_SigmaR_04-08-2017.mat');

% other functions
addpath('functions')

% Phase of the low-frequency transducer, measured using fibre-optic
% hydrophone  see script: "scope_traces_determine_phase_shift_20231023.m"
addpath('../Fibre-optic hydrophone measurements 20231020')
if contains(data_name, '1.9')
    phase_LF = load("20231020_1d9MHz_1000rep_PRF10Hz.mat");
elseif contains(data_name, '1.8')
    phase_LF = load("20231020_1d8MHz_1000rep_PRF10Hz_2.mat");
end

input_param.phase_LF = phase_LF.phase;

%% Initialize
selIdx = find(strainVal>0.1 & SNR'>=25 & strainVal < 0.2);

% array of which bubble(s) to analyse
wfmIndx_array = selIdx(1:end);


%% =======settings (S)===================================================================
data_name = erase(data_name, '_V2');
S.data_name = data_name;

% booleans for plotting
S.plot_all = 0;


% 1= filtering the strain around the fundamental and 2nd harmonic. 0 = no filtering
S.filtering = 1;

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
pressure_string = extractBetween(data_name, ['Vial' vial_name{1} '_'], 'kPa');
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
% (2 µs after start HF to 16 µs before end HF)    
D.time_full = -35.08e-6:P.dt:((size(D.time_temp,2)-1)*P.dt + 21.92e-6);        % vial B 
D.time_full = -35.08e-6:P.dt:((size(D.time_temp,2)-1)*P.dt + 21.92e-6);         % vial C

% make envelope (pad with zeros to match the time of the strain: +- 4 µs)
D.env_hanning_temp = hann(round(length(D.time_temp)))';
D.envelop_hanning_full = [zeros(1, round(abs(D.time_full(1))/P.dt))...
    D.env_hanning_temp...
    zeros(1, length(D.time_temp) - length(D.env_hanning_temp))...
    zeros(1, round(abs(D.time_full(end)-200e-6)/P.dt))];

D.Pacc_full = P.Pa*D.envelop_hanning_full.*sin(2*pi*P.fUS.*D.time_full);    % incident acoustic pressure (P.Pa)

% set P.T from -5 to 105 µs for the cropped signals (only ramp-up)
D.time_up = -5e-6: P.dt: 105e-6;


% crop the strain and driving pressure to 0 - 200 us
[~, crop_idx_start] = min(abs(D.time_full - 0));
[~ , crop_idx_end] = min(abs(D.time_full - 200e-6));
D.envelop_hanning = D.envelop_hanning_full(crop_idx_start:crop_idx_end);
D.time = D.time_full(crop_idx_start:crop_idx_end);
D.Pacc = D.Pacc_full(crop_idx_start:crop_idx_end);


%% ======== Response (R) ssimulated with Rayleigh-Plesset equation ======================
% Determine the surface tension curve from Segers
input_param.strain_filtered = -0.2: S.strain_input_step: 0.2;
[F] = SS_func_cropping_surface_tension_Segers(P.R0*1e-6, P.sig_0, ...
    input_param, S, F);

RP_input.time = D.time_full;
RP_input.P_ac = D.Pacc_full;
RP_input.sig_0 = P.sig_0;
RP_input.R0 = P.R0*1e-6;
RP_input.strain_par = F.strain_input_crop;
RP_input.surf_tens_par = F.surf_tension_Segers_crop;

tic
initial_conditions = [RP_input.R0, 0];
options = odeset('BDF','on','AbsTol', 1e-7, 'RelTol', 1e-4); %default AbsTol is 1e-6, RelTol 1e-3, BDF off
%relTol: This tolerance measures the error relative to the magnitude of each solution component.
%absTol: This tolerance is a threshold below which the value of the solution (y) becomes unimportant. The value of AbsTol should take into account the scale of the solution components.

[t, x] = ode45(@(t,x) RP_eq_cropping(t, x, P, D, RP_input),...
    RP_input.time, initial_conditions, options);
toc
    
    
R.time = t';
R.radius = x(:,1)';
R.velocity = x(:,2)';


% strain
R.strain = (R.radius - RP_input.R0) / RP_input.R0;

figure();
plot(D.time_full*1e6, R.strain)
grid on
xlabel('time (µs)')
ylabel('simulated ')

%% Analyse every bubble one by one ----------------------------------------
warning('off')

lag_max_corr.samples = NaN(1, max(wfmIndx_array));
lag_max_corr.us = NaN(1, max(wfmIndx_array));

for wfmIndx = wfmIndx_array(1:end)

    if ismember(wfmIndx, wfmIndx_array)

        saving_folder = [cd '/error minimization data/' data_name ];
        saving_dir = [saving_folder ...           
            '/bubble ' num2str(wfmIndx)];
        saving_file_name = [saving_dir '.mat'];


        fprintf('\nMeasurement index: %d\n', wfmIndx)


        %% Filtering spectra around the fundamental and 2nd harmonic ----------
        [strain_filtered_temp, f_fft_appended] = SS_func_strain_filtering(strain_full, wfmIndx, S, P, D);


        %% Differentiate in the Fourier domain --------------------------------
        %(see script 'diff_fourier.m' for tests)
        [input_param.strain_filtered, input_param.dstrain_dt, input_param.d2strain_dt2] ...
            = SS_func_strain_derivatives(strain_full, strain_filtered_temp, f_fft_appended, wfmIndx, S, P, D);
        
        %% Cross-correlate measured strain with simulated strain  ------------------------
        % Cross-correlation measures the similarity between a vector x and shifted 
        % (lagged) copies of a vector y as a function of the lag
        % the lag is thus the shift necessary for the measured strain to
        % match the driving pressure (+ is to later point in time).
        [corr, lags] = xcorr(R.strain/max(R.strain)*max(input_param.strain_filtered), input_param.strain_filtered);
        [max_corr, max_corr_idx] = max(corr);

        lag_max_corr.samples(wfmIndx) = lags(max_corr_idx);      % in samples
        lag_max_corr.us(wfmIndx) = lags(max_corr_idx) / P.Fs;      % in s



    end
end

% lag for the time vector of the driving pressure, so that driving pressure
% and measured strains match in time (based on xcorr simulated strain and
% measured strains)
median(lag_max_corr.us, 'omitnan')

std(lag_max_corr.us, 'omitnan')
mean(lag_max_corr.us, 'omitnan')

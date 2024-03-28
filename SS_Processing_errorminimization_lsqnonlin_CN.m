% Script to analyse the stress-strain data as obtained by Sander Spiekhout
% and Charlotte Nawijn
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
strain_full = relWfms;   %strainY(:,1:51801)

% Parametrization Tim Segers
addpath('../Stress_strain solver simulation/functions')
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

%% Initialize
selIdx = find(strainVal>0.1 & SNR'>=25 & strainVal < 0.2);% & R0est>2.2e-6 & R0est<2.4e-6);% & strainVal <0.2);  %hannStrainVal is the average strain value over the entire 200 µs LF pulse (dus niet op de piek!)        %abs(1-(refLvl2./refLvl1))<=0.3) is already in the data (abs(1-(refLvl2./refLvl1)) = FluctVal?)

% array of which bubble(s) to analyse
wfmIndx_array = selIdx(1:end);

%% =======settings (S)===================================================================
data_name = erase(data_name, '_V2');
S.data_name = data_name;

% booleans for plotting, saving, etc.
if S.calculate_on_server
    S.save_data = 1;
    S.saving_figs = 0;
    S.do_error_minimization = 1;        % to plot the surface tension figures of the results on the server, disable this
    
    S.plot_all = 0;
else
    S.save_data = 0;
    S.saving_figs = 0;
    S.do_error_minimization = 0;        % to plot the surface tension figures of the results on the server, disable this

    S.plot_all = 0;
end

% How to determine the elastic contribution. 1= from the line at zero
% strain rate of the full visocelatic plane, 0 = from the points at the envelope of 
% the strain filled into the expression at zero strain rate
S.elastic_from_full_viscoelastic = 0;

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
P.cw = 1500;                       % speed of sound in water (m/s)
P.P0 = 1e5;                        % ambient pressure (P.Pa)
P.mu = 2*1e-3;                     % viscosity of water (P.Pa*s) multiplied by two to account for thermal damping %1e-3;
P.rho = 1e3;                       % density of water (kg/m^3)
P.sig = 0.072;                     % surface tension water (N/m)

% bubble variables
P.dR0 = 0;                         % initial velocity bubble wall (m/s)
P.kap = 1.07;                      % polytropic exponent C3F8 gas
% P.chi = 0.5;      % Tim: 0.55    % P.shell stiffness (N/m)
% P.ks = 1e-8;                     % P.shell viscosity (P.Pa*s)
P.shell = 1;                       % boolean whether the bubble has a P.shell
if P.shell == 1
    P.sig_0 = 30e-3;  % Tim: 15e-3   % initial surface tension of the P.shell (N/m)
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


%% ====== Time for strain, see 'SS_Processing_cropping_strain_CN_20231127.m' 
% to see how I came to the offsets of D.time_full (based on cross-correlation
% simulated strain and measured strains

% crop the data
time_strain = D.time_full;        % time vector for strain

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

for wfmIndx = fliplr(wfmIndx_array(1: end))
    close all

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
            num2str(S.mask_cutoff) '/' data_name];

        if ~S.elastic_from_full_viscoelastic
            saving_folder = [saving_folder '/elas_dir'];
        end

        saving_dir = [saving_folder ...
            '/bubble ' num2str(wfmIndx)];
        saving_file_name = [saving_dir '.mat'];

        if S.calculate_on_server || (~S.calculate_on_server && S.do_error_minimization)
            if exist(saving_file_name, 'file') == 2
                continue
            end
        elseif ~S.calculate_on_server && ~S.do_error_minimization

            fig_folder = [cd '/error minimization figures/'];

            if S.section_analysis == 1
                fig_folder = [fig_folder 'updown/'];
            elseif S.section_analysis == 2
                fig_folder = [fig_folder 'up/'];
            elseif S.section_analysis ==3
                fig_folder = [fig_folder 'down/'];
            end
           
            data_name = erase(data_name, '_V2');

            fig_folder = [fig_folder 'range_EM_' num2str(S.range_EM) '_mask_cutoff_'...
                num2str(S.mask_cutoff) '/' data_name];

            if ~S.elastic_from_full_viscoelastic
                fig_folder = [fig_folder '/elas_dir'];
            end

            filename_fig = [fig_folder '\bubble_' num2str(wfmIndx) '_bestfit.png'];
            if exist(filename_fig, 'file')
                continue
            end

            if ~(exist(saving_file_name, 'file') == 2)
                continue
            end
        end

        
        %% Filtering spectra around the fundamental and 2nd harmonic ----------
        [strain_filtered_temp, f_fft_appended] = ...
            SS_func_strain_filtering(strain_full, wfmIndx, S, P, D);


        %% Differentiate in the Fourier domain --------------------------------
        %(see script 'diff_fourier.m' for tests)
        [input_param.strain_filtered, input_param.dstrain_dt, input_param.d2strain_dt2] ...
            = SS_func_strain_derivatives(strain_full, strain_filtered_temp, f_fft_appended, wfmIndx, S, P, D);


        % crop the strain and driving to 0 - 200 us
        input_param.strain_filtered = input_param.strain_filtered(crop_idx_start:crop_idx_end);
        input_param.dstrain_dt = input_param.dstrain_dt(crop_idx_start:crop_idx_end);
        input_param.d2strain_dt2 = input_param.d2strain_dt2(crop_idx_start:crop_idx_end);


        % median of this range over the minima and maxima determine the
        % lower and upper limits of the strain vector, respectively
        S.range_vec = round(length(input_param.strain_filtered)*0.01);

        %% Error minimization (3D) --------------------------------------------
        if S.do_error_minimization

            input_param.regime1 = zeros(1, round(S.range_EM/S.strain_input_step));
            input_param.regime2 = P.sig.*ones(1, round(S.range_EM/S.strain_input_step));

            t_timer_start = tic;

            S.plot_all = 0;

            continue_next = NaN(1, S.attempts);

            % for lsqnonlin and lsqcurvefit with 'trust-region-reflective' (default), the
            % function tolerance is relative, while the step and optimality
            % tolerances are absolute. For 'levenberg-marquardt' they are all
            % relative
            opts = optimoptions(@lsqnonlin, 'Display', 'off',... 'final-detailed',...
                'FunctionTolerance', 1e-12, ... %1e-12);%,...
                'OptimalityTolerance', 1e-30, ...%1e-8
                'StepTolerance', 1e-7, 'DiffMinChange', 0.2, ...%1e-22, ...);%1e-10
                'MaxFunctionEvaluations', 5e2, 'MaxIterations', 1e4);%,...
            % 'PlotFcn', 'optimplotx');%'PlotFcn', 'optimplotx', optimplotfval


            initial_guess = NaN(S.attempts, 4);
            input_param.weights = [1e6, 1e2, 5e6, 1e-4];

            input_param.R0_Sander = R0est(wfmIndx); % in m

           
            % initial guess for the fitting parameters, chosen randomly
            % within the lower and upper bounds

            % R0 (m), sig_0 (N/m), shift (driving pressure), pressure
            % amplitude, (damping constant)
            S.lb = [1e-6 0 -P.T*2/3 P.Pa];
            S.ub = [5e-6 72e-3 P.T/3 P.Pa];

            % For debugging
            % S.lb = [2e-6 5e-3 -1e-7, P.Pa];
            % S.ub = S.lb;


            S.lb = S.lb.*input_param.weights;
            S.ub = S.ub.*input_param.weights;

            initial_guess(:,1) = (S.ub(1) - S.lb(1)) .* rand(1, S.attempts) + S.lb(1); %ones(1, S.attempts)*R0_Sander;
            initial_guess(:,2) = (S.ub(2) - S.lb(2)) .* rand(1, S.attempts) + S.lb(2);
            initial_guess(:,3) = (S.ub(3) - S.lb(3)) .* rand(1, S.attempts) + S.lb(3);
            initial_guess(:,4) = (S.ub(4) - S.lb(4)) .* rand(1, S.attempts) + S.lb(4);


            xestimated = NaN(size(initial_guess));
            resnorm = NaN(S.attempts, 1);
            residuals = NaN(S.attempts, 2*round(S.range_EM/S.strain_input_step));
            exitflag = NaN(S.attempts, 1);
            % output = NaN(S.attempts, 1);

            % -------------------------------------------------------------
            % error minimization function over two regimes (maximum no. of points per
            % regime is round(S.range_EM/S.strain_input_step)
            fun = @(x) SS_func_error_lsqnonlin(x, input_param,...
                wfmIndx, S, P, D, F);
            %-------------------------------------------------------------


            parfor i = 1: S.attempts
                [xestimated(i,:), resnorm(i), residuals(i,:), exitflag(i), output(i)] = ...
                    lsqnonlin(fun, initial_guess(i,:), S.lb, S.ub, opts);
            end

            % Move to next bubble if there is no iteration for which there
            % are enough points to do the fitting
            if sum(residuals(:,1) == 1e5) == S.attempts
                fprintf('\nWarning! Could not perform the fitting on this bubble')   
                break
            end

            % find minimum resnorm
            [min_resnorms, min_resnorm_indices] = mink(resnorm, ceil(S.attempts/10));
            [min_resnorm, min_resnorm_idx] = min(resnorm);

            [min_residuals, min_residual_indices] = mink(sum(abs(residuals), 2), ceil(S.attempts/10));
            [min_residual, min_residual_idx] = min(sum(abs(residuals), 2));

            xestimated = xestimated./input_param.weights;

            error_min_R0 = median(xestimated(min_resnorm_indices, 1));
            error_min_sig_0  = median(xestimated(min_resnorm_indices, 2));
            error_min_delay = median(xestimated(min_resnorm_indices, 3));
            error_min_Pa = median(xestimated(min_resnorm_indices, 4));

            % show fitting results, residuals
            t_timer=toc(t_timer_start);

            fprintf(['\nMeasurement index: %d\n' ...
                'R0est = %.3f µm\n' ...
                'Result fit min resnorm (%d attempts)\nR0 = %.3f µm\n' ...
                'sig0 = %.3f mN/m\ndelay = %.3f µs (%.3f pi)\n'...
                'Pa = %.3f kPa\n',...
                'Total elapsed time is %.2f minute(s)\n'], ...
                wfmIndx, R0est(wfmIndx)*1e6, S.attempts, ...
                error_min_R0*1e6, error_min_sig_0*1e3, error_min_delay*1e6, ...
                error_min_delay/(P.T)*2, error_min_Pa*1e-3,...vestimated(min_resnorm_idx,4),
                t_timer/60)

            warning('on')

        else
            saving_folder = [cd '/error minimization data/'];

            if S.section_analysis == 1
                saving_folder = [saving_folder 'updown/'];
            elseif S.section_analysis == 2
                saving_folder = [saving_folder 'up/'];
            elseif S.section_analysis == 3
                saving_folder = [saving_folder 'down/'];
            end

            saving_folder = [saving_folder...
                'range_EM_' num2str(S.range_EM) '_mask_cutoff_'...
                num2str(S.mask_cutoff) '/' data_name];
             
            filename = [saving_folder '\bubble ' num2str(wfmIndx) '.mat'];

            % filename_fig = [pwd '\error minimization figures\'...
            %     data_name '\bubble_' num2str(wfmIndx) '_bestfit.png'];
            if exist(filename, 'file') % && ~exist(filename_fig, 'file')
                S_temp = S;

                load(filename)

                S_temp.lb = S.lb;
                S_temp.ub = S.ub;
                S=S_temp;
                clear S_temp   % to not overwrite the settings
            else
                continue
            end
        end


        %% Apply the R0, sigma(R0), delay, pressure amplitude found from error minimization,
        % and determine the surface tension curve

        [strain_vec, dstrain_vec, strain_grid, dstrain_grid,...
            strain_vec_crop, surf_tension, strain_input_crop, surf_tension_crop,...
            surf_tension_smooth, surf_tension_old,...
            viscoelastic_contrib_int, viscoelastic_contrib_int_old,...
            elastic_contrib, elastic_contrib_old]...
            = SS_func_surface_tension(error_min_R0,...
            error_min_sig_0, error_min_delay, error_min_Pa, ...
            input_param, wfmIndx, S, P, D, data_name);

        %% Determine the surface tension curve from Segers
        [F] = SS_func_surface_tension_Segers(error_min_R0, error_min_sig_0, ...
            input_param, S, F);

   

    end
end

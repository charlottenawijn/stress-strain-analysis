% Script to analyse the results of the error minimization of the
% stress-strain data as obtained by Sander Spiekhout
% using the acoustical camera (Erasmus University, Rotterdam) to obtain
% bubble shell parameters: shell elasticity (chi = A dsigma/dA) from the
% elastic contribution, and shell viscosity (kappa_s) from the viscous
% contribution

% University of Twente,
% Charlotte Nawijn, 2021 - 2024

clear, clc, close all

global sig_0 R0 fit_surf_tens


%% Settings
% which values of the maximum strain (thus not Sander's strainVal
% parameter) to consider. The processing script (stress-strain analysis
% already contains a selection: 0.1 < strainVal < 0.2)
S_analysis.strainValMin = 0.1;
S_analysis.strainValMax = 0.2;

S_analysis.max_error = inf;

S_analysis.SNR_threshold = 25;

% load data Sander
% 2.5 µm radius bubbles, heated (vial C)
dataset_name_general{1} = 'StrainSignals_VialC_35kPa_1.8MHz_Mono_2.5um_16082023_V2';
dataset_folder{1} = 'vial C_1.8 MHz_35 kPa';

% 2.3 µm radius bubbles, heated (vial B)
dataset_name_general{2} = 'StrainSignals_VialB_40kPa_1.9MHz_Mono_2.3um_17082023_V2';
dataset_folder{2} = 'vial B_1.9 MHz_40 kPa';


S_analysis.mask_cutoff = 0.1;       % select the version corresponding to which cutoff of the processed data to show
S_analysis.range_EM = 0.04;

% 1= analyse full ramp-up and -down in pressure amplitude,
% 2= only ramp-up,
% 3= only ramp-down
S_analysis.section_analysis = 1;

% How to determine the elastic contribution. 1= from the line at zero
% strain rate of the full visocelatic plane, 0 = from the points at the envelope of 
% the strain filled into the expression at zero strain rate
S.elastic_from_full_viscoelastic = 1;




%% initialize
[surf_tens_EM_array, elasticity_EM_array, strain_vec_EM_array, initial_surface_tension_EM_array,...
    pressure_amplitude_EM_array] = deal({});
number_of_bubbles = NaN(1,1);


addpath([pwd '\..\Stress_strain solver simulation'])
fit_surf_tens = load('..\Stress_strain solver simulation\Parametrisatie oppervlaktespanning (Tim)\fit_SigmaR_04-08-2017.mat');


%% Load data into a 3D array
for data_idx = 2 %1:2

    data_name = erase(dataset_name_general{data_idx}, '_V2');

    % saving path for data
    saving_path_data = [cd '/error minimization data/'];

    if S_analysis.section_analysis == 1
        saving_path_data = [saving_path_data 'updown/'];
    elseif S_analysis.section_analysis == 2
        saving_path_data = [saving_path_data 'up/'];
    elseif S_analysis.section_analysis ==3
        saving_path_data = [saving_path_data 'down/'];
    end
    
    saving_path_data = [saving_path_data...
            'range_EM_' num2str(S_analysis.range_EM) '_mask_cutoff_'...
            num2str(S_analysis.mask_cutoff) '/' data_name];
    if ~S.elastic_from_full_viscoelastic
        saving_path_data = [saving_path_data '/elas_dir'];
    end

    if ~exist(saving_path_data)
        mkdir(saving_path_data)
    end

    vial_name = extractBetween(dataset_name_general{data_idx}, 'Vial', '_');
    vial_name = vial_name{1};

    % saving path for figures
    saving_path_fig = replace(saving_path_data, ...
        'error minimization data', 'analysis plots');

    if ~exist(saving_path_fig)
        mkdir(saving_path_fig)
    end


    % save settings
    save([saving_path_fig '\settings.mat'], 'S_analysis')

    % saving figures for paper
    saving_paper_path{data_idx} = [cd '\..\..\..\Artikelen\Stress-strain\figures\Results\'...
        'raw plot\' erase(erase(dataset_name_general{data_idx}, 'StrainSignals_'), '_V2') '\'];
    if S_analysis.section_analysis == 1
        saving_paper_path{data_idx} = [saving_paper_path{data_idx} 'updown/'];
    elseif S_analysis.section_analysis == 2
        saving_paper_path{data_idx} = [saving_paper_path{data_idx} 'up/'];
    elseif S_analysis.section_analysis ==3
        saving_paper_path{data_idx} = [saving_paper_path{data_idx} 'down/'];
    end

    if ~S.elastic_from_full_viscoelastic
        saving_paper_path{data_idx} = [saving_paper_path{data_idx} 'elas_dir/'];
    end


    if ~exist(saving_paper_path{data_idx}, 'dir')
        mkdir(saving_paper_path{data_idx})
    end

    dataset_dir = dataset_name_general{data_idx};

    load(['..\Data and experiment prep\August 2023 experiments\'...
        dataset_folder{data_idx} '\' dataset_dir])


    number_of_bubbles(data_idx) = size(relWfms,1);        %without making assumptions on strainVal

    % clear wfms hanndR
    R0_diff = NaN(1, number_of_bubbles(data_idx));


    selIdx{data_idx} = [];
    
    % initialize
    [surf_tens_EM, surf_tens_EM_os, elasticity_surf_tens, strain_central, strain_filtered_EM,...
        surf_tens_par, strain_vec_par, radius_from_strain_EM, visc_elas_EM, ...
        strain_vec_EM, strain_vec_full_EM, dstrain_vec_EM] = deal(cell(1, number_of_bubbles(data_idx)));

    [error_EM_min, initial_radii_EM, initial_surface_tension_EM, elasticity_init, elasticity_init_Segers, delay_EM, ...
        pressure_amplitude_EM] = deal(NaN(1, number_of_bubbles(data_idx)));

    RMSE = NaN(1, number_of_bubbles(data_idx));


    pressure_before_dip = NaN(1, number_of_bubbles(data_idx));
    time_before_dip = NaN(1, number_of_bubbles(data_idx));
    EC_ratio = NaN(1, number_of_bubbles(data_idx));

    for wfmIndx = 1: number_of_bubbles(data_idx)
        filename = [saving_path_data '\bubble ' num2str(wfmIndx) '.mat'];
        if exist(filename, 'file') == 2 %&& ismember(wfmIndx, selIdx{data_idx})%(1:44))

            close all

            % load data
            load(filename)

            R0 = error_min_R0;
            sig_0 = error_min_sig_0;


            mean_strain_step(wfmIndx) = mean(diff(strain_vec), 'omitnan');

            R0_diff(wfmIndx) = abs(R0est(wfmIndx) - R0);        % difference in sizes in m


            if strainVal(wfmIndx) > S_analysis.strainValMin && strainVal(wfmIndx) < S_analysis.strainValMax ...
                && SNR(wfmIndx) >= S_analysis.SNR_threshold
                    selIdx{data_idx} = [selIdx{data_idx} wfmIndx];
            end

            %%% expansion over compression ratio from the filtered strain
            expansion = mean(maxk(input_param.strain_filtered, 5));
            compression = mean(mink(input_param.strain_filtered, 5));
            EC_ratio(wfmIndx) = expansion/abs(compression);
            %%%

            visc_elas_EM{data_idx, wfmIndx} = viscoelastic_contrib_int;

            error_EM_min(data_idx, wfmIndx) = median(min_resnorms);
            
            initial_radii_EM(data_idx, wfmIndx) = error_min_R0;
            initial_surface_tension_EM(data_idx, wfmIndx) = error_min_sig_0;
            delay_EM(data_idx, wfmIndx) = error_min_delay;
            pressure_amplitude_EM(data_idx, wfmIndx) = error_min_Pa;

            strain_filtered_EM{data_idx, wfmIndx} = input_param.strain_filtered;


            % surface tension from error minimization
            surf_tens_EM{data_idx, wfmIndx} = surf_tension_crop; %surf_tension;
            strain_vec_EM{data_idx, wfmIndx} = strain_vec_crop; %strain_vec;
            
            strain_vec_full_EM{data_idx, wfmIndx} = strain_vec;
            dstrain_vec_EM{data_idx, wfmIndx} = dstrain_vec;


            % interpolate (linearly) to match dimension of measured curve (Segers)
            surf_tens_EM_os{data_idx, wfmIndx} = interp1(strain_vec_EM{data_idx, wfmIndx},...         % os=oversampled
                surf_tens_EM{data_idx, wfmIndx}, F.strain_input_crop);

        
            % calculate RMS error of the EM result w.r.t. measured curve (Segers)
            surf_tens_par{data_idx, wfmIndx} = F.surf_tension_Segers_crop;
            strain_vec_par{data_idx, wfmIndx} = F.strain_input_crop;
            

            RMSE(data_idx, wfmIndx) = sqrt(mean((surf_tens_EM_os{data_idx, wfmIndx}*1e3 -...
                surf_tens_par{data_idx, wfmIndx}*1e3).^2, 'omitnan'));


            %% elasticity: from slope near zero strain; fit linear line
            % elasticity of the back-calculated surface tension
            [~, zero_strain_idx] = min(abs(strain_vec_EM{data_idx, wfmIndx}));
            range_idx = round(0.01 / mean(diff(strain_vec_EM{data_idx, wfmIndx})));      % < 5 kPa would give a maximum strain of ~0.008

            start_idx = zero_strain_idx-range_idx;
            end_idx = zero_strain_idx+range_idx;
            lin_fit = polyfit(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx),...
                surf_tens_EM{data_idx, wfmIndx}(start_idx: end_idx), 1);
           elasticity_init(data_idx, wfmIndx) = (1+mean(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx)))/2 .* lin_fit(1);        % lin_fit(1) is slope (derivative of surface tension vs radius plot)

            % for the entire strain region      
            smooth_range = length(dstrain_vec_EM{data_idx, wfmIndx})/10;
            surf_tens_smooth = smooth(surf_tens_EM{data_idx, wfmIndx}, smooth_range);

            elasticity_surf_tens{data_idx, wfmIndx}(1: range_idx) = NaN;
            elasticity_surf_tens_smooth{data_idx, wfmIndx}(1: range_idx) = NaN;


            for i = 1+range_idx: length(surf_tens_EM{data_idx, wfmIndx}) - range_idx
                clear lin_fit lin_fit_smooth

                start_idx = i-range_idx;
                end_idx = i+range_idx;
                lin_fit = polyfit(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx), ...
                    surf_tens_EM{data_idx, wfmIndx}(start_idx: end_idx), 1);

                lin_fit_smooth = polyfit(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx), ...
                    surf_tens_smooth(start_idx: end_idx), 1);

                elasticity_surf_tens{data_idx, wfmIndx}(i) = (1+mean(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx)))/2 ...
                    .* lin_fit(1);           % in N/m
                elasticity_surf_tens_smooth{data_idx, wfmIndx}(i) = (1+mean(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx)))/2 ...
                    .* lin_fit_smooth(1);           % in N/m

            end

            elasticity_surf_tens{data_idx, wfmIndx}( 1+ ...
                length(surf_tens_EM{data_idx, wfmIndx}) - range_idx : ...
                length(surf_tens_EM{data_idx, wfmIndx})) = NaN;
            elasticity_surf_tens_smooth{data_idx, wfmIndx}( 1+ ...
                length(surf_tens_EM{data_idx, wfmIndx}) - range_idx : ...
                length(surf_tens_EM{data_idx, wfmIndx})) = NaN;


            [~, zero_strain_idx] = min(abs(F.strain_input_crop));
            range_idx = round(0.01 / mean(diff(F.strain_input_crop)));      % < 5 kPa would give a maximum strain of ~0.008

            start_idx = zero_strain_idx-range_idx;
            end_idx = zero_strain_idx+range_idx;
            lin_fit = polyfit(F.strain_input_crop(start_idx: end_idx), ...
                F.surf_tension_Segers_crop(start_idx: end_idx), 1);
            elasticity_init_Segers(data_idx, wfmIndx) = (1+mean(F.strain_input_crop(start_idx: end_idx)))/2 .* lin_fit(1);        % lin_fit(1) is slope (derivative of surface tension vs radius plot)

            surf_tens_init_Segers(data_idx, wfmIndx) = F.surf_tension_Segers_crop(zero_strain_idx);
        end
    end

    % selIdx{data_idx} = selIdx{data_idx}(selIdx{data_idx}>118);


    %% concatenate
    % strain_vec_EM_array{data_idx} = strain_vec_EM_temp(selIdx{data_idx});
    surf_tens_EM_array{data_idx} = cat(1, surf_tens_EM{data_idx, selIdx{data_idx}});
    elasticity_EM_array{data_idx} = cat(1, elasticity_surf_tens{data_idx, selIdx{data_idx}});
    elasticity_smooth_EM_array{data_idx} = cat(1, elasticity_surf_tens_smooth{data_idx, selIdx{data_idx}});
    strain_vec_EM_array{data_idx} = cat(1, strain_vec_EM{data_idx,selIdx{data_idx}});

    initial_radii_EM_array{data_idx} = initial_radii_EM(data_idx,selIdx{data_idx});
    delay_EM_array{data_idx} = delay_EM(data_idx,selIdx{data_idx});
    error_EM_array{data_idx} = error_EM_min(data_idx, selIdx{data_idx});
    initial_surface_tension_EM_array{data_idx} = initial_surface_tension_EM(data_idx, selIdx{data_idx});
    initial_elasticity_EM_array{data_idx} = elasticity_init(data_idx, selIdx{data_idx});
    pressure_amplitude_EM_array{data_idx} = pressure_amplitude_EM(data_idx, selIdx{data_idx});

    RMSE_array{data_idx} = RMSE(data_idx,selIdx{data_idx});


%%
median(initial_surface_tension_EM_array{data_idx})*1e3
std(initial_surface_tension_EM_array{data_idx})*1e3

median(initial_elasticity_EM_array{data_idx})
std(initial_elasticity_EM_array{data_idx})

clear R0


%% ! Plot all surface tensions together (with non-dimensional strain on x-axis)

% match the indices of the strain with the values between bubbles
min_strain_EM = min(min([strain_vec_EM_array{data_idx}]));
min_strain_EM=min(min_strain_EM, -0.2522);

max_strain_EM = max(max(strain_vec_EM_array{data_idx}));
strain_full_array = min_strain_EM: 1e-4: max_strain_EM;    % array from min to max strain for all bubbles

% k >= the max, min and mean value of all surface tensions at a certain
% strain
count = zeros(1,length(strain_full_array));
surf_tens_temp = [];

for j = 1: length(strain_full_array)            % for all strain values
    strain_temp = strain_full_array(j);

    for i = 1: size(strain_vec_EM_array{data_idx}, 1)       %per bubble       
        % find if there is this strain for each bubble
        % [~, k] = find(abs(strain_vec_EM_array{data_idx}(i,:) - strain_temp) <= 1e-3, 1, 'first');
        [~, k] = find(abs(strain_vec_EM_array{data_idx}(i,:) - strain_temp) <= 2e-3, 1, 'first');
        if ~isempty(k)
            % all surface tension values at this strain
            surf_tens_temp(j,i) = surf_tens_EM_array{data_idx}(i,k);
            count(j) = count(j)+1;
        else
            surf_tens_temp(j,i) = NaN;
        end
    end
end


%% determine mean, median, standard deviation of the surface tension
mean_surf_tens_full_array = mean(surf_tens_temp, 2, 'omitnan')';
median_surf_tens_full_array = median(surf_tens_temp, 2, 'omitnan')';
min_surf_tens_full_array = min(surf_tens_temp, [], 2)';
max_surf_tens_full_array = max(surf_tens_temp, [], 2)';
std_surf_tens_full_array = std(surf_tens_temp, 0, 2, 'omitnan')';
std_region_min = mean_surf_tens_full_array - std_surf_tens_full_array;
std_region_plus = mean_surf_tens_full_array + std_surf_tens_full_array;

%indices between which min and plus are not the same (gives trouble with
%patch
% start_idx_std = find(std_region_plus > std_region_min, 1, 'first');
% end_idx_std = find(std_region_plus > std_region_min, 1, 'last');
start_idx_std = find(std_region_plus >= std_region_min, 1, 'first');
end_idx_std = find(std_region_plus >= std_region_min, 1, 'last');


std_region_min(isnan(std_region_min)) = 0;
std_region_plus(isnan(std_region_plus)) = 0;



%measured curve (Segers)
% surface tension from measured curve (Segers)
R0 = median(initial_radii_EM_array{data_idx}, 'omitnan');
sig_0 = median(initial_surface_tension_EM_array{data_idx}, 'omitnan');

R_input_param_step = 1e-4;
% R_input_param_par = min(min([strain_filtered_EM{:}] * R0 + R0))*1e6:...
%     R_input_param_step: max(max([strain_filtered_EM{:}] * R0 + R0))*1e6; % in µm
R_input_param_par = (min(min([strain_vec_EM_array{data_idx}])) * R0 + R0)*1e6:...
    R_input_param_step: (max_strain_EM * R0 + R0)*1e6; % in µm
[A0c, sigma_par_input_param] = A0cor3(R_input_param_par);

% shift input_param surface tension sigma_par_input_param to overlap with initial
% surface tension (caused by difference in bubble radius, but shape
% measured curve (Segers) is the same for all sizes (see mail Tim 20/7/2021)
shift_sigma_um = R0*1e6 - R_input_param_par(find(sigma_par_input_param >...
    mean(initial_surface_tension_EM_array{data_idx}, 'omitnan'), 1, 'first') -1);
shift_sigma = round(shift_sigma_um / R_input_param_step);
if shift_sigma >= 0
    sigma_par_input_param_shifted = zeros(size(sigma_par_input_param));
    sigma_par_input_param_shifted(shift_sigma+1: end)...
        = sigma_par_input_param(1: end-shift_sigma);      % shift to the right
elseif shift_sigma < 0
    sigma_par_input_param_shifted = zeros(size(sigma_par_input_param));
    sigma_par_input_param_shifted(1: shift_sigma+end)...
        = sigma_par_input_param(-shift_sigma+1: end);      % shift to the right
    sigma_par_input_param_shifted(end+shift_sigma:end) = sigma_par_input_param_shifted(end+shift_sigma-1);
end
surf_tens_par_combi = sigma_par_input_param_shifted;


% median instead of mean
std_region_min = median_surf_tens_full_array - std_surf_tens_full_array;
std_region_plus = median_surf_tens_full_array + std_surf_tens_full_array;
%indices between which min and plus are not the same (gives trouble with
%patch
start_idx_std = find(std_region_plus >= std_region_min, 1, 'first');
end_idx_std = find(std_region_plus >= std_region_min, 1, 'last');



std_region_min(isnan(std_region_min)) = 0;
std_region_plus(isnan(std_region_plus)) = 0;



%% elasticity of median surface tension
% initial elasticity of median: from slope near zero strain
[~, zero_strain_idx] = min(abs(strain_full_array));
range_idx = round(0.01 / mean(diff(strain_full_array)));      % < 5 kPa would give a maximum strain of ~0.008
surf_tens_smooth = smooth(median_surf_tens_full_array, 100);

clear lin_fit
start_idx = zero_strain_idx - range_idx;
end_idx = zero_strain_idx + range_idx;
     
lin_fit = polyfit(strain_full_array(start_idx: end_idx), ...
         median_surf_tens_full_array(start_idx: end_idx), 1);

elasticity_init_med_surf_tens(data_idx) = (1+mean(strain_full_array(start_idx: end_idx)))/2 ...
    .* lin_fit(1);           % in N/m



% elasticity of median surface tension for entire range of strain
elasticity_med_surf = [];

range_idx = round(0.01 / mean(diff(strain_full_array)));      % < 5 kPa would give a maximum strain of ~0.008

elasticity_med_surf(1: range_idx) = NaN;

 for i = range_idx+1: length(strain_full_array)-range_idx
     clear lin_fit
     start_idx = i-range_idx;
     end_idx = i+range_idx;

     lin_fit = polyfit(strain_full_array(start_idx: end_idx), ...
         median_surf_tens_full_array(start_idx: end_idx), 1);
     elasticity_med_surf(i) = (1+mean(strain_full_array(start_idx: end_idx)))/2 ...
         .* lin_fit(1);        % lin_fit(1) is slope (derivative of surface tension vs radius plot)
 end


elasticity_med_surf(1 + length(strain_full_array)-range_idx : length(strain_full_array)) = NaN;

%% elasticity of Tim's curve
elasticity_Segers = [];

strain_Segers = (R_input_param_par - R0*1e6)/ (R0*1e6);

range_idx = round(0.01 / mean(diff(strain_Segers)));      % < 5 kPa would give a maximum strain of ~0.008

elasticity_Segers(1: range_idx) = NaN;

 for i = range_idx+1: length(strain_Segers)-range_idx
     clear lin_fit
     start_idx = i-range_idx;
     end_idx = i+range_idx;

     lin_fit = polyfit(strain_Segers(start_idx: end_idx), ...
         surf_tens_par_combi(start_idx: end_idx), 1);
     elasticity_Segers(i) = (1+mean(strain_Segers(start_idx: end_idx)))/2 .* lin_fit(1);        % lin_fit(1) is slope (derivative of surface tension vs radius plot)

 end

elasticity_Segers(1 + length(strain_Segers)-range_idx : length(strain_Segers)) = NaN;


%% ! Plot all elasticities together (with non-dimensional strain on x-axis)

% match the indices of the strain with the values between bubbles
%     strain_vec_EM_array = cat(3,strain_vec_EM{data_idx,selIdx_new});
min_strain_EM = min(min([strain_vec_EM_array{data_idx}]));
min_strain_EM=min(min_strain_EM, -0.2522);

max_strain_EM = max(max(strain_vec_EM_array{data_idx}));
strain_full_array = min_strain_EM: 1e-4: max_strain_EM;    % array from min to max strain for all bubbles

% k >= the max, min and mean value of all surface tensions at a certain
% strain
count = zeros(1,length(strain_full_array));
elasticity_temp = [];

for j = 1: length(strain_full_array)            % for all strain values
    strain_temp = strain_full_array(j);

    for i = 1: size(strain_vec_EM_array{data_idx}, 1)       %per bubble       
        % find if there is this strain for each bubble
        % [~, k] = find(abs(strain_vec_EM_array{data_idx}(i,:) - strain_temp) <= 1e-3, 1, 'first');
        [~, k] = find(abs(strain_vec_EM_array{data_idx}(i,:) - strain_temp) <= 2e-3, 1, 'first');

        if ~isempty(k)
            % all surface tension values at this strain
            elasticity_temp(j,i) = elasticity_EM_array{data_idx}(i,k);
            elasticity_smooth_temp(j,i) = elasticity_smooth_EM_array{data_idx}(i,k); 

            count(j) = count(j)+1;
        else
            elasticity_temp(j,i) = NaN;
            elasticity_smooth_temp(j,i) = NaN;
        end
    end
end


mean_elas_full_array = mean(elasticity_temp, 2, 'omitnan')';
median_elas_full_array = median(elasticity_temp, 2, 'omitnan')';
min_elas_full_array = min(elasticity_temp, [], 2)';
max_elas_full_array = max(elasticity_temp, [], 2)';
std_elas_full_array = std(elasticity_temp, 0, 2, 'omitnan')';


%% histograms with results error minimization, also error as function of the,
%% initial radius R0
% % size estimation directly from acoustical camera (Sander)
[counts_AC, centers_AC] = hist(R0est(selIdx{data_idx}), round(1+3.322*log(length(R0est(selIdx{data_idx})))));

% from error minimization
[counts(data_idx,:), centers(data_idx, :)] = hist(initial_radii_EM_array{data_idx}, ...
    round(1+3.322*log(length(initial_radii_EM_array{data_idx}))));

% from Coulter counter
coulter_data_path = '../Data and experiment prep/August 2023 experiments/Coulter data/with subtraction of lipids/Acoustical Camera bubbles 202307end Coulter 20230817';
coulter_all_files = dir(fullfile([coulter_data_path '/vial' vial_name '*_50ul_*.xls']));
coulter_start_idx = str2num(coulter_all_files(1).name(end-10:end-8));
coulter_end_idx = coulter_start_idx+2; %str2num(coulter_all_files(end).name(end-10:end-8));

[bin_radius, bin_count] = deal(cell(1,1));

for c = 1: length(coulter_all_files)
    coulter_file_name = [coulter_data_path '\' coulter_all_files(c).name];
    if exist(coulter_file_name) && ~contains(coulter_file_name, 'bground') && ~contains(coulter_file_name, 'destroyed')
        Coulter_data_temp = readtable([fileparts(cd) coulter_file_name(3:end)], 'FileType', 'text');

        if sum(strcmp('Diff_', Coulter_data_temp.Properties.VariableNames))>0
            bin_radius{c} = Coulter_data_temp.BinDiameter/2;
            bin_count{c} = Coulter_data_temp.Diff_;
        else
            Coulter_data_temp = table2array(Coulter_data_temp);
            bin_radius{c} = Coulter_data_temp(1:length(bin_radius{1}),2)./2;
            bin_count{c} = Coulter_data_temp(1:length(bin_radius{1}),3);
        end

    end
end

% free lipids (Coultered after 2x 1 minute in the ultrasonic bath)
[bin_radius_lipids, bin_count_lipids] = deal(cell(1,1));

for c=1: length(coulter_all_files)
    coulter_file_name = [coulter_data_path '\' coulter_all_files(c).name];
    if contains(coulter_file_name, 'destroyed')
        Coulter_data_temp = readtable([fileparts(cd) coulter_file_name(3:end)], 'FileType', 'text');

        if sum(strcmp('Diff_', Coulter_data_temp.Properties.VariableNames))>0
            bin_radius_lipids{c} = Coulter_data_temp.BinDiameter/2;
            bin_count_lipids{c} = Coulter_data_temp.Diff_;
        else
            Coulter_data_temp = table2array(Coulter_data_temp);
            bin_radius_lipids{c} = Coulter_data_temp(:,2)./2;
            bin_count_lipids{c} = Coulter_data_temp(1:length(bin_radius{1}),3);
        end
    end
end
bin_count_lipids_temp = [bin_count_lipids{:}];
bin_count_lipids = reshape(bin_count_lipids_temp, [], 3);


% median size and standard deviation of Coulter
Coulter_count_temp = mean([bin_count{:}], 2)-mean(bin_count_lipids, 2);
Coulter_count = Coulter_count_temp(Coulter_count_temp>=0);
Coulter_radii = bin_radius{1}(Coulter_count_temp>=0);

for i = 1: length(Coulter_radii)
    Coulter_all{i} = repmat(Coulter_radii(i), 1, round(Coulter_count(i)*1000));
end
Coulter_med_size = median([Coulter_all{:}]); % in µm
Coulter_mean_size = mean([Coulter_all{:}]); % in µm
Coulter_std_size = std([Coulter_all{:}]); % in µm


end

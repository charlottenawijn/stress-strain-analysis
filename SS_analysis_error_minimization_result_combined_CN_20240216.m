% Script to analyse the results of the error minimization of the
% stress-strain data as obtained by Sander Spiekhout
% using the acoustical camera (Erasmus University, Rotterdam) to obtain
% bubble shell parameters: shell elasticity (chi = A dsigma/dA) from the
% elastic contribution, and shell viscosity (kappa_s) from the viscous
% contribution

% University of Twente,
% Charlotte Nawijn, 2021 - 2023

clear, clc, close all

global sig_0 R0 fit_surf_tens


%% Settings
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

saving_figs = 0;
saving_figs_paper_combined = 0;
saving_figs_paper_examples = 0;
addpath('..\..\..\Artikelen\Stress-strain')
  

plot_individual = 0;        % in figure 4, plot every line individually to make a video, see line below
make_vid = 0;

markersize = 20;
plot_position = [550 300 550 420];%[574 439 650 420];

% which values of the maximum strain (thus not Sander's strainVal
% parameter) to consider. The processing script (stress-strain analysis
% already contains a selection: 0.1 < strainVal < 0.2)
S_analysis.strainValMin = 0.1;
S_analysis.strainValMax = 0.2;%0.2; % more strict to ensure that the ramp up is not done too quickly?, to have good sampling!

S_analysis.max_error = inf;%0.05;

S_analysis.SNR_threshold = 25;%25;%27.5;

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
S.elastic_from_full_viscoelastic = 0;




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
    % selIdx{data_idx} = find(strainVal>0.08 & SNR'>=25);% & strainVal < 0.2 & R0est>2.2e-6 & R0est<2.4e-6);% & strainVal <0.2);  %hannStrainVal is the average strain value over the entire 200 µs LF pulse (dus niet op de piek!)        %abs(1-(refLvl2./refLvl1))<=0.3) is already in the data (abs(1-(refLvl2./refLvl1)) = FluctVal?)

    
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


            % load fit to surface tension (parametrization) from Segers
            % fit_surf_tens = load('..\Stress_strain solver simulation\Parametrisatie oppervlaktespanning (Tim)\fit_SigmaR_04-08-2017.mat');


            % error_min_R0 = median(vestimated(vestimated(:,1)<2.4, 1))*1e-6;           
            % error_min_sig_0 = median(vestimated(vestimated(:,1)<2.4, 2))*1e-3;
            % error_min_delay = median(vestimated(vestimated(:,1)<2.4, 3));


            R0 = error_min_R0;
            sig_0 = error_min_sig_0;




            mean_strain_step(wfmIndx) = mean(diff(strain_vec), 'omitnan');

            R0_diff(wfmIndx) = abs(R0est(wfmIndx) - R0);        % difference in sizes in m


            % if max(strain_vec, [], 2) > strain_min && max(strain_vec, [], 2) < strain_max    % thus of the strain BEFORE cutoff
            %     selIdx{data_idx} = [selIdx{data_idx} wfmIndx];
            % end
            if strainVal(wfmIndx) > S_analysis.strainValMin && strainVal(wfmIndx) < S_analysis.strainValMax ...
                && SNR(wfmIndx) >= S_analysis.SNR_threshold
                    selIdx{data_idx} = [selIdx{data_idx} wfmIndx];
            end

            %%% expansion over compression ratio from the filtered strain
            expansion = mean(maxk(input_param.strain_filtered, 5));
            compression = mean(mink(input_param.strain_filtered, 5));
            EC_ratio(wfmIndx) = expansion/abs(compression);
            %%%

            % if exist('viscoelastic_contrib_int', 'var')
                visc_elas_EM{data_idx, wfmIndx} = viscoelastic_contrib_int;
            % else
            %     visc_elas_EM{data_idx, wfmIndx} = viscous_contrib_int;
            % end
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

            % radius_from_strain_EM{1,wfmIndx} = radius_from_strain;

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
            % deriv_init_surf_tens = diff([surf_tens_EM{data_idx, wfmIndx}(zero_strain_idx-range_idx)...
            %     surf_tens_EM{data_idx, wfmIndx}(zero_strain_idx+range_idx)]) ...
            %     ./ (diff([strain_vec_EM{data_idx, wfmIndx}(zero_strain_idx-range_idx)...
            %    strain_vec_EM{data_idx, wfmIndx}(zero_strain_idx+range_idx)]*R0 + R0));
            % elasticity_init(data_idx, wfmIndx) = R0./2 .* deriv_init_surf_tens;           % in N/m
   
            start_idx = zero_strain_idx-range_idx;
            end_idx = zero_strain_idx+range_idx;
            lin_fit = polyfit(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx),...
                surf_tens_EM{data_idx, wfmIndx}(start_idx: end_idx), 1);
           elasticity_init(data_idx, wfmIndx) = (1+mean(strain_vec_EM{data_idx, wfmIndx}(start_idx: end_idx)))/2 .* lin_fit(1);        % lin_fit(1) is slope (derivative of surface tension vs radius plot)

           % figure; plot(strain_vec_EM{data_idx, wfmIndx}*R0 + R0 , surf_tens_EM{data_idx, wfmIndx})
           % hold on;
           % plot(strain_vec_EM{data_idx, wfmIndx}*R0 + R0, lin_fit(1)*(strain_vec_EM{data_idx, wfmIndx}*R0 + R0) + lin_fit(2), 'r-')

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

                % elasticity_surf_tens{data_idx, wfmIndx}(i) = R0./2 .* lin_fit(1);           % in N/m
                % elasticity_surf_tens_smooth{data_idx, wfmIndx}(i) = R0./2 .* lin_fit_smooth(1);           % in N/m

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


            
            % figure(2001)
            % clf
            % plot(strain_vec_EM{data_idx, wfmIndx}, ...
            %     elasticity_surf_tens{data_idx, wfmIndx}, 'linewidth', 1.5)
            % grid on
            % xlabel('\epsilon')
            % ylabel('\chi (N/m)')
            % yline(elasticity_init(data_idx, wfmIndx), 'r-', 'linewidth', 1.5)
            % hold on
            % plot(strain_vec_EM{data_idx, wfmIndx}, ...
            %     elasticity_surf_tens_smooth{data_idx, wfmIndx}, 'linewidth', 1.5)

            %%
            % (initial) elasticity of Segers curve
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
    % surf_tens_EM_temp = surf_tens_EM{data_idx, :};
    % surf_tens_EM_array{data_idx} = surf_tens_EM_temp(selIdx{data_idx});
    % 
    % strain_vec_EM_temp = strain_vec_EM{data_idx};
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

% save([saving_path_fig '\settings.mat'], 'S_analysis', 'S')


%%
median(initial_surface_tension_EM_array{data_idx})*1e3
std(initial_surface_tension_EM_array{data_idx})*1e3

median(initial_elasticity_EM_array{data_idx})
std(initial_elasticity_EM_array{data_idx})

%% save

% save([saving_path_data '/initial_radii.mat'], 'initial_radii_EM_array')
% save([saving_path_data '/initial_surface_tension.mat'], 'initial_surface_tension_EM_array')
% save([saving_path_data '/initial_elasticity.mat'], 'initial_elasticity_EM_array')
% save([saving_path_data '/delay.mat'], 'delay_EM_array')
% save([saving_path_data '/pressure_amplitude.mat'], 'pressure_amplitude_EM_array')
% save([saving_path_data '/error_final.mat'], 'error_EM_array')

clear R0

%% saving figure of strain and surface tension for paper
if saving_figs_paper_examples && data_idx == 2

    fig_size = 0.6;%0.5;

    wfmIndx_sel = fliplr([4 32 72]);

    for wfmIndx = wfmIndx_sel
        
        % strain
        fig_strain = figure(1000+wfmIndx);
        clf
        plot(D.time*1e6, strain_filtered_EM{data_idx, wfmIndx}, 'linewidth', 1)%, 'color', 'k')
        ylabel('\epsilon')
        xlabel('time (µs)')
        % title('radial strain')
        grid on
        % xlim([min(D.time_full*1e6) max(D.time_full*1e6)])
        xlim([min(D.time*1e6) max(D.time*1e6)])
        ylim([-0.25 0.25])
        xticks(0:50:200)
        yticks(-0.3:0.1:0.3)
        set(gca,'box','on')

        [fig_strain, ~] = Plotting_size_settings(fig_strain, gca, fig_size);
          
        % saving
            % filepath = ['C:\Users\NawijnCL\OneDrive - Universiteit Twente\Artikelen\'...
            %     'Stress-strain\figures'...
            %     '\Results\' erase(dataset_name_general{data_idx}, 'StrainSignals_')...
            %     '\examples strain surface tension'];
            filepath = [saving_paper_path{data_idx} 'examples strain surface tension'];

            if ~exist(filepath)
                mkdir(filepath)
            end
            
            saving_path = ['\bubble_' num2str(wfmIndx) '_strain'];
            
            saveas(fig_strain, [filepath, saving_path, '.fig'])
            exportgraphics(fig_strain, [filepath, saving_path, '.pdf'], 'ContentType','vector')
            exportgraphics(fig_strain, [filepath, saving_path, '.eps'], 'ContentType','vector')
            exportgraphics(fig_strain, [filepath, saving_path, '.png'])
            saveas(fig_strain, [filepath, saving_path, '.svg'])

      % viscoelastic contribution
        fig_viscelas = figure(1001+wfmIndx);
        clf
        p_surf = surf(strain_vec_full_EM{data_idx, wfmIndx}, dstrain_vec_EM{data_idx, wfmIndx}, ...
            visc_elas_EM{data_idx, wfmIndx},...
            'edgecolor', 'none');
        xlabel('\epsilon')
        ylabel('$\dot{\epsilon}$', 'interpreter', 'latex')
        grid on

        addpath('functions\cbrewer')
        % title('viscoelastic contribution')
        xlim([-0.22 0.22])
        ylim([-1.4 1.4])
        xticks(-0.2:0.1:0.2)
        yticks(-2:1:2)
        view(0,90)
        set(gca,'box','on')

        [cmap] = cbrewer('div', 'RdBu', 100);
        colormap(cmap)
        clim([-max(abs([visc_elas_EM{data_idx, wfmIndx}]),[],'all') ...
            max(abs([visc_elas_EM{data_idx, wfmIndx}]),[],'all')])
        [fig_viscelas, ~] = Plotting_size_settings(fig_viscelas, gca, fig_size);

          % saving           
          saving_path = ['\bubble_' num2str(wfmIndx) '_viscoelastic_contrib'];

          saveas(fig_viscelas, [filepath, saving_path, '.fig'])
          exportgraphics(fig_viscelas, [filepath, saving_path, '.pdf'], 'ContentType','vector')
          exportgraphics(fig_viscelas, [filepath, saving_path, '.eps'], 'ContentType','vector')
          exportgraphics(fig_viscelas, [filepath, saving_path, '.png'])
          saveas(fig_viscelas, [filepath, saving_path, '.svg'])



          % add colorbar and save again
          c=colorbar();
          c.Title.String = '$\tilde{P}_\mathrm{VE}$';
          c.Title.Interpreter = 'latex';
          c.Title.FontName = 'Arial';
          c.Position(4) = 0.65;
          c.Position(1) = 0.8;
          c.Position(2) = 0.2;
          set(gca, 'visible', 'off')
          p_surf.Visible = 'off';

        % set(gca, 'fontsize', 10)

        
        % saving           
            saving_path = ['\bubble_' num2str(wfmIndx) '_viscoelastic_contrib_colorbar'];
            
            saveas(fig_viscelas, [filepath, saving_path, '.fig'])
            exportgraphics(fig_viscelas, [filepath, saving_path, '.pdf'], 'ContentType','vector')
            exportgraphics(fig_viscelas, [filepath, saving_path, '.eps'], 'ContentType','vector')
            exportgraphics(fig_viscelas, [filepath, saving_path, '.png'])
            saveas(fig_viscelas, [filepath, saving_path, '.svg'])


        % surface tension
        fig_surf_tens = figure(1002+wfmIndx);
        clf
        plot(strain_vec_EM{data_idx, wfmIndx}, surf_tens_EM{data_idx, wfmIndx}*1e3,...
            'linewidth', 1)
        hold on
        plot(strain_vec_par{data_idx, wfmIndx}, surf_tens_par{data_idx, wfmIndx}*1e3, '--',...
            'color', [0.8500 0.3250 0.0980], 'linewidth', 1.33)
        xlabel('\epsilon')
        ylabel('\sigma (mN/m)')
        % title('surface tension')
        grid on
        ylim([-20 100])
        % xlim([-0.16 0.16])
        % yticks(0:50:100)
        yticks(-40: 40: 100)
        % xticks(-0.2:0.1:0.2)
        % yticks(0:20:100)
        xlim([-0.22 0.22])
        xticks(-0.3:0.1:0.3)
        if wfmIndx == wfmIndx_sel(1)
            legend({'data single bubble', ['similar population' newline '(Segers et al.)']}, 'location', 'northwest')
        end      
        set(gca,'box','on')

    
        [fig_surf_tens, ~] = Plotting_size_settings(fig_surf_tens, gca, fig_size);
          
        % set(gca, 'fontsize', 10)

        
        % saving           
            saving_path = ['\bubble_' num2str(wfmIndx) '_surface_tension'];
            
            saveas(fig_surf_tens, [filepath, saving_path, '.fig'])
            exportgraphics(fig_surf_tens, [filepath, saving_path, '.pdf'], 'ContentType','vector')
            exportgraphics(fig_surf_tens, [filepath, saving_path, '.eps'], 'ContentType','vector')
            exportgraphics(fig_surf_tens, [filepath, saving_path, '.png'])
            saveas(fig_surf_tens, [filepath, saving_path, '.svg'])

    end

end


%% ! Plot all surface tensions together (with non-dimensional strain on x-axis)

% match the indices of the strain with the values between bubbles
%     strain_vec_EM_array = cat(3,strain_vec_EM{data_idx,selIdx_new});
%% 
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

% if no data, make count 0
strain_full_array_plot = [min_strain_EM strain_full_array max_strain_EM -min_strain_EM];
count = [0 count 0 0];



%% show how many data points (bubbles) we have at each value of strain
fig_count_strain=figure(1);
clf
plot(strain_full_array_plot, count, 'linewidth', 1.5)
xlabel('\epsilon')
ylabel('count')
grid on
% xlim([-0.22 0.22])
xlim([min_strain_EM -min_strain_EM])
xticks(-0.4:0.1:0.4)
set(gca, 'box', 'on')

[fig_count_strain, ~] = Plotting_size_settings(fig_count_strain, gca, 0.6);

if saving_figs_paper_combined  
    % saving
    saving_path = [saving_paper_path{data_idx} 'count_per_strain'];

    
    saveas(fig_count_strain, [saving_path, '.fig'])
    exportgraphics(fig_count_strain, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_count_strain, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_count_strain, [saving_path, '.png'])
    saveas(fig_count_strain, [saving_path, '.svg'])
end

if saving_figs
    saveas(gcf, [saving_path_fig...
        '/count vs strain'], 'png')
    saveas(gcf, [saving_path_fig...
        '/count vs strain'], 'fig')
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


figure(5)
clf
set(gcf, 'Position', plot_position)
hold on
set(gca, 'fontsize', 14)
xlabel('strain (R-R_0)/R_0')
ylabel('surface tension (mN/m)')
title(['surface tension as a function of strain (N=' num2str(length(selIdx{data_idx})) ')'])
grid on
plot(strain_full_array, mean_surf_tens_full_array*1e3, 'linewidth', 2)
p=patch([strain_full_array(start_idx_std:end_idx_std)...
    fliplr(strain_full_array(start_idx_std:end_idx_std))],...
    [std_region_min(start_idx_std:end_idx_std)*1e3 ...
    fliplr(std_region_plus(start_idx_std:end_idx_std)*1e3)],...
    [0 0.4470 0.7410]);

p.FaceAlpha = 0.5;
p.EdgeColor = 'none';
plot((R_input_param_par - R0*1e6)/ (R0*1e6), surf_tens_par_combi*1e3, '--',...
    'color', [0.8500 0.3250 0.0980], 'linewidth', 1.75)
% plot(strain_full_array, min_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% plot(strain_full_array, max_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% legend('mean', 'std', 'min', 'max', 'location','northwest')
legend('mean', 'standard deviation', 'measured curve (Segers)', 'location', 'northwest')
% xlim([-0.15 0.15])
if saving_figs
    saveas(gcf, [saving_path_fig...
        '/surface tension mean+-std vs strain'], 'png')
    saveas(gcf, [saving_path_fig...
        '/surface tension mean+-std vs strain'], 'fig')
end

% median instead of mean
std_region_min = median_surf_tens_full_array - std_surf_tens_full_array;
std_region_plus = median_surf_tens_full_array + std_surf_tens_full_array;
%indices between which min and plus are not the same (gives trouble with
%patch
start_idx_std = find(std_region_plus >= std_region_min, 1, 'first');
end_idx_std = find(std_region_plus >= std_region_min, 1, 'last');



std_region_min(isnan(std_region_min)) = 0;
std_region_plus(isnan(std_region_plus)) = 0;

%% plot surface tension vs strain for paper
fig_size = 1;% 0.8;

fig_surf_tens_med = figure(50);
clf
set(gcf, 'Position', plot_position)
hold on
set(gca, 'fontsize', 14)
xlabel('\epsilon')
% xlabel('radial strain \epsilon')
ylabel('\sigma (mN/m)')
% title(['surface tension (n=' num2str(length(selIdx{data_idx})) ')'])
grid on
% plot(strain_full_array, mean_surf_tens_full_array*1e3, 'linewidth', 2)
% plot(strain_full_array, median_surf_tens_full_array*1e3, 'linewidth', 1.5)%, 'color', 'k')
plot(strain_full_array, median_surf_tens_full_array*1e3, 'linewidth', 1.5)%, 'color', 'k')
p=patch([strain_full_array(start_idx_std:end_idx_std)...
    fliplr(strain_full_array(start_idx_std:end_idx_std))],...
    [std_region_min(start_idx_std:end_idx_std)*1e3 ...
    fliplr(std_region_plus(start_idx_std:end_idx_std)*1e3)],...
    [0 0.4470 0.7410]);
p.FaceAlpha = 0.35;
p.EdgeColor = 'none';

strain_Segers = (R_input_param_par - R0*1e6)/ (R0*1e6);

plot(strain_Segers, surf_tens_par_combi*1e3, '--',...
    'color', [0.8500 0.3250 0.0980], 'linewidth', 1.5)
% plot(strain_full_array, min_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% plot(strain_full_array, max_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% legend('mean', 'std', 'min', 'max', 'location','northwest')
legend({'median', ['standard deviation'], ...
    ['similar population' newline '(Segers et al.)']},...
    'location', 'northwest')
% xlim([-0.22 0.22])
xlim([min_strain_EM -min_strain_EM])
ylim([-20 100])
xticks(-0.2:0.1:0.2)
yticks(-20: 20: 100)
set(gca,'box','on')

% change plot size, font size etc.
[fig_surf_tens_med, ~] = Plotting_size_settings(fig_surf_tens_med, gca, fig_size);


if saving_figs_paper_combined

    % set(gca, 'fontsize', 10)

    saving_path = [saving_paper_path{data_idx} 'surface_tension_median'];

    saveas(fig_surf_tens_med, [saving_path, '.fig'])
    exportgraphics(fig_surf_tens_med, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_surf_tens_med, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_surf_tens_med, [saving_path, '.png'])
    saveas(fig_surf_tens_med, [saving_path, '.svg'])

end

% xlim([-0.3 0.3])
% ylim([-40 120])
if saving_figs
    saveas(gcf, [saving_path_fig...
        '/surface tension median+-std vs strain'], 'png')
    saveas(gcf, [saving_path_fig...
        '/surface tension median+-std vs strain'], 'fig')
end



%% Standard deviation surface tension vs. strain
figure(500)
set(gcf, 'Position', plot_position)
hold on
plot(strain_full_array, std_surf_tens_full_array)
title(['standard deviation surface tension (N=' num2str(length(selIdx{data_idx})) ')'])
set(gca, 'fontsize', 14)
xlabel('radial strain (R-R_0)/R_0')
ylabel('standard deviation in surface tension (mN/m)')
xlim([-0.3 0.3])


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


 %% plot elasticity (of median surface tension) vs strain for paper
fig_elas_med = figure(1000);
clf
set(gcf, 'Position', plot_position)
hold on
set(gca, 'fontsize', 14)
xlabel('\epsilon')
ylabel('\chi (N/m)')
% title('elasticity')
grid on
plot(strain_full_array, elasticity_med_surf, 'linewidth', 1)%, 'color', 'k')
% plot(strain_central, elasticity_whole, 'linewidth', 1.5)%, 'color', 'k')
% plot(strain_central_Segers, elasticity_Segers, '--',...
%     'color', [0.8500 0.3250 0.0980], 'linewidth', 2)
plot(strain_Segers, elasticity_Segers, '--',...
    'color', [0.8500 0.3250 0.0980], 'linewidth', 1.5)
% plot(strain_full_array, min_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% plot(strain_full_array, max_surf_tens_full_array*1e3, '--', 'linewidth', 1)
% xlim([-0.22 0.22])
xlim([min_strain_EM -min_strain_EM])
% ylim([-20 100])
ylim([-0.5 0.8])
yticks(-0.8: 0.4: 0.8)
set(gca,'box','on')
% legend({'of median', ['of similar population' newline '(Segers et al.)']},...
%     'location', 'best')
legend({'median', ...
    ['similar population' newline '(Segers et al.)']},...
    'location', 'southwest')
        set(gca,'box','on')


% change plot size, font size etc.
fig_size = 1;%0.8;
[fig_elas_med, ~] = Plotting_size_settings(fig_elas_med, gca, fig_size);


if saving_figs_paper_combined
 
    saving_path = [saving_paper_path{data_idx} 'elasticity'];

    saveas(fig_elas_med, [saving_path, '.fig'])
    exportgraphics(fig_elas_med, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_elas_med, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_elas_med, [saving_path, '.png'])
    saveas(fig_elas_med, [saving_path, '.svg'])

end



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

std_region_min = median_elas_full_array - std_elas_full_array;
std_region_plus = median_elas_full_array + std_elas_full_array;
%indices between which min and plus are not the same (gives trouble with
%patch
start_idx_std = find(std_region_plus >= std_region_min, 1, 'first');
end_idx_std = find(std_region_plus >= std_region_min, 1, 'last');

std_region_min(isnan(std_region_min)) = 0;
std_region_plus(isnan(std_region_plus)) = 0;


%% plot elasticity vs. strain for paper
fig_elas = figure(2000);
clf
set(gcf, 'Position', plot_position)
hold on
set(gca, 'fontsize', 14)
xlabel('\epsilon')
% xlabel('radial strain \epsilon')
ylabel('\chi (N/m)')
% title(['surface tension (n=' num2str(length(selIdx{data_idx})) ')'])
grid on
plot(strain_full_array, median_elas_full_array, 'linewidth', 1.5)%, 'color', 'k')
p=patch([strain_full_array(start_idx_std:end_idx_std)...
    fliplr(strain_full_array(start_idx_std:end_idx_std))],...
    [std_region_min(start_idx_std:end_idx_std) ...
    fliplr(std_region_plus(start_idx_std:end_idx_std))],...
    [0 0.4470 0.7410]);
p.FaceAlpha = 0.35;
p.EdgeColor = 'none';

plot(strain_Segers, elasticity_Segers, '--',...
    'color', [0.8500 0.3250 0.0980], 'linewidth', 1.5)

% plot(strain_full_array, median(elasticity_smooth_temp, 2, 'omitnan')', 'linewidth', 1)

% legend({'median', ['standard deviation'], ...
%     ['similar population' newline '(Segers et al.)']},...
%     'location', 'southwest')
% xlim([-0.22 0.22])
xlim([min_strain_EM -min_strain_EM])
% ylim([-20 100])
xticks(-0.2:0.1:0.2)
% yticks(-20: 20: 100)
ylim([-0.5 0.8])
yticks(-0.8: 0.2: 0.8)
set(gca,'box','on')

% change plot size, font size etc.
fig_size = 1;%0.8;
[fig_elas, ~] = Plotting_size_settings(fig_elas, gca, fig_size);


if saving_figs_paper_combined

    saving_path = [saving_paper_path{data_idx} 'median_elasticity'];

    saveas(fig_elas, [saving_path, '.fig'])
    exportgraphics(fig_elas, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_elas, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_elas, [saving_path, '.png'])
    saveas(fig_elas, [saving_path, '.svg'])

end

if saving_figs
    saveas(gcf, [saving_path_fig...
        '/elasticity median+-std vs strain'], 'png')
    saveas(gcf, [saving_path_fig...
        '/elasticity median+-std vs strain'], 'fig')
end

%% Make video of all individual surface tension curves

if make_vid
    %     vidfile = VideoWriter([saving_path_fig '/surface tension individual.mp4'], 'MPEG-4');
    vidfile = VideoWriter('analysis plots\surface tension individual.mp4', 'MPEG-4');
    open(vidfile);

    close(figure(40))
    figure(40)
    set(gcf, 'Position', plot_position)

    xlabel('radial strain (R-R_0)/R_0')
    ylabel('surface tension (mN/m)')
    title('surface tension vs strain for all bubbles')
    grid on


    for i = 1:data_idx
        for j = 1: size(strain_vec_EM_array{i},1)
            plot(strain_vec_EM_array{i}(j,:), surf_tens_EM_array{i}(j,:)*1e3, 'linewidth', 2)
            hold on



            % find max strain
            [max_strain, i_max] = max(strain_vec_EM_array{i}(j, :));
            plot(max_strain, surf_tens_EM_array{i}(j,i_max)*1e3, 'r.', 'markersize', 15)
            %         xlim([-0.15 0.15])

            plot((R_input_param_par - R0*1e6)/ (R0*1e6), surf_tens_par_combi*1e3, '--k',...
                'linewidth', 1.75)


            %         xlabel('strain (R-R_0)/R_0')
            %         ylabel('surface tension (mN/m)')
            title(['surface tension vs strain of bubble ' num2str(selIdx{data_idx}(j))])
            xlabel('radial strain (R-R_0)/R_0')
            ylabel('surface tension (mN/m)')
            grid on
            legend('data', 'max strain', 'measured curve (Segers)', 'location', 'northwest')
            set(gcf,'color','w')

            ylim([1.2*min(surf_tens_EM_array{i}(:)*1e3, [], 'all', 'omitnan')...
                1.2*max(surf_tens_EM_array{i}(:)*1e3, [], 'all', 'omitnan')])
            xlim([1.2*min(strain_vec_EM_array{i}(:), [], 'all', 'omitnan')...
                1.2*max(strain_vec_EM_array{i}(:), [], 'all', 'omitnan')])

            pause(0.001)
            drawnow

            im(j) = getframe(gcf);
            writeVideo(vidfile, im(j));

            hold off

        end
    end

    writeVideo(vidfile, getframe(gcf));
    close(vidfile)
    status = movefile(['analysis plots\' vidfile.Filename], [saving_path_fig '/surface tension individual.mp4']);
    if status == 0
        disp('transfer video unsuccessful')
    end
end


if plot_individual
    figure(4)
    set(gcf, 'Position', plot_position)
    xlabel('radial strain (R-R_0)/R_0')
    ylabel('surface tension (mN/m)')
    title('surface tension vs strain for all bubbles')
    grid on

    plot((R_input_param_par - R0*1e6)/ (R0*1e6), surf_tens_par_combi*1e3, '--k',...
        'linewidth', 1.75)

    for i = 1:data_idx
        for j = 1: size(strain_vec_EM_array{i},1)
            %         figure(400+j)
            %         hold on
            plot(strain_vec_EM_array{i}(j,:), surf_tens_EM_array{i}(j,:)*1e3)
            pause(0.01)
            % find max strain
            [max_strain, i_max] = max(strain_vec_EM_array{i}(j,~isnan(surf_tens_EM_array{i}(j,:))));
            plot(max_strain, surf_tens_EM_array{i}(j,i_max)*1e3, 'r.', 'markersize', 15)
            %         xlim([-0.15 0.15])
            drawnow
            %         xlabel('strain (R-R_0)/R_0')
            %         ylabel('surface tension (mN/m)')
            title('surface tension vs strain for all bubbles')
            grid on

            if make_vid
                im(j) = getframe(gcf);
                writeVideo(vidfile, im(j));
            end
        end
    end


    plot((R_input_param_par - R0*1e6)/ (R0*1e6), surf_tens_par_combi*1e3, '--k',...
        'linewidth', 1.75)


    if saving_figs
        saveas(gcf, [saving_path_fig...
            '/surface tension vs strain for all bubbles'], 'png')
        saveas(gcf, [saving_path_fig...
            '/surface tension vs strain for all bubbles'], 'fig')
    end

end




%% histograms with results error minimization, also error as function of the,
%% initial radius R0
% % size estimation directly from acoustical camera (Sander)
% [counts_AC, centers_AC] = hist(R0est(R0est<5e-6), round(1+3.322*log(length(R0est(R0est<5e-6)))));
% [counts_AC, centers_AC] = hist(R0est, round(1+3.322*log(length(R0est))));
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

        % plot(bin_radius{c}, bin_count{c}, 'linewidth', 1.25)%'DisplayName', [num2str(c+coulter_start_idx-1) ''])
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
Coulter_med_size = median([Coulter_all{:}]); % in um
Coulter_std_size = std([Coulter_all{:}]);

% figure(7)
% clf(7)
% set(gcf, 'Position', [350 250 1000 600])
% hold on
% title({'size distribution', dataset_dir}, 'interpreter', 'none')
% subplot(3,1,1)
% bar(centers*1e6, counts)
% xlabel('radius (µm)')
% ylabel('count')
% grid on
% xlim([1 5])
% title('error minimization (lsqcurvefit)')
% subplot(3,1,2)
% 
% bar(bin_radius{1}, mean([bin_count{:}], 2)-mean(bin_count_lipids, 2), 'linewidth', 2)
% 
% xlim([1 5])
% ylim([0 inf])
% ylabel('count')
% grid on
% title('Coulter')
% subplot(3,1,3)
% bar(centers_AC*1e6, counts_AC)
% xlabel('radius (µm)')
% ylabel('count')
% grid on
% xlim([1 5])
% title('acoustical camera sizing (<5 µm)')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/hist size distribution separate'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/hist size distribution separate'], 'fig')
% end



%% Size distribution: plot in a single plot
color = colororder;

size_fig = 0.6;

fig_size_distrib = figure(51);
clf(51)
hold on
set(gcf, 'Position', plot_position)
% if size_fig == 0.5
   plot(centers(data_idx, :)*1e6, counts(data_idx,:)/max(counts(data_idx,:)), 'color', 'k', 'linewidth', 1.5)
% elseif size_fig == 1
   % plot(centers(data_idx, :)*1e6, counts(data_idx,:)/max(counts(data_idx,:)), 'color', 'k', 'linewidth', 2)
% end
hold on
% plot(centers_AC*1e6, counts_AC/max(counts_AC))
xlabel('radius (µm)')
ylabel('count (normalized)')
% if size_fig == 0.5
    plot(Coulter_radii, Coulter_count./max(Coulter_count), '--', 'color', color(2,:),...
        'linewidth', 1)
% elseif size_fig == 1
%     plot(Coulter_radii, Coulter_count./max(Coulter_count), '--', 'color', color(2,:),...
%         'linewidth', 1.5)
% end
xlim([1 5])
ylim([0 1])
grid on
title('size distribution')
% title({'size distribution', dataset_name}, 'interpreter', 'none')
% legend({'EM', 'AC sizing', ['Coulter' newline '(free lipids subtr.)']})
legend({'fitting', ['Coulter' newline 'counter']}, 'location', 'northeast')
% set(gca,'fontsize', 14)
yticks(0:0.25:1)

% change plot size, font size etc.
[fig_size_distrib, ~] = Plotting_size_settings(fig_size_distrib, gca, size_fig);

if saving_figs_paper_combined
   
    saving_path = [saving_paper_path{data_idx} 'size_distribution'];

    saveas(fig_size_distrib, [saving_path, '.fig'])
    exportgraphics(fig_size_distrib, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_size_distrib, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_size_distrib, [saving_path, '.png'])
    saveas(fig_size_distrib, [saving_path, '.svg'])

end


if saving_figs
    saveas(gcf, [saving_path_fig...
        '/hist size distribution'], 'png')
    saveas(gcf, [saving_path_fig...
        '/hist size distribution'], 'fig')
end


%% delay (bubble position) 
[counts_phase, centers_phase] = hist(delay_EM_array{data_idx},...
    round(1+3.322*log(length(delay_EM_array{data_idx}))));

size_fig = 0.6;


fig_delay = figure(52);
set(gcf, 'Position', plot_position)
hold on
% bar(centers_phase*1e6, counts_phase)
% if size_fig == 0.5
    plot(centers_phase*1e6, counts_phase/max(counts_phase), 'color', 'k', 'linewidth', 1.5)
% elseif size_fig == 1
    % plot(centers_phase*1e6, counts_phase/max(counts_phase), 'color', 'k', 'linewidth', 2)
% end
% text(centers_phase*1e6, counts_phase, num2str(counts_phase'),'vert','bottom','horiz','center');
xlabel('delay position (µs)')
ylabel('count (normalized)')
grid on
title('delay')
set(gca, 'fontsize', 14)
% xlim([S.lb(3)/input_param.weights(3)*1e6, S.ub(3)/input_param.weights(3)*1e6])
xlim([-0.4, 0.2])
ylim([0 1])
yticks(0:0.25:1)

% change plot size, font size etc.
[fig_delay, ~] = Plotting_size_settings(fig_delay, gca, size_fig);

if saving_figs_paper_combined
   
    saving_path = [saving_paper_path{data_idx} 'delay'];

    saveas(fig_delay, [saving_path, '.fig'])
    exportgraphics(fig_delay, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_delay, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_delay, [saving_path, '.png'])
    saveas(fig_delay, [saving_path, '.svg'])

end

if saving_figs
    saveas(gcf, [saving_path_fig...
        '/hist phase delay driving pressure'], 'png')
    saveas(gcf, [saving_path_fig...
        '/hist phase delay driving pressure'], 'fig')
end



%% initial surface tension
[counts_surf_tens, centers_surf_tens] =...
    hist(initial_surface_tension_EM_array{data_idx}, ...
    round(1+3.322*log(length(initial_surface_tension_EM_array{data_idx}))));

% if no data, make count 0
centers_surf_tens = [0 centers_surf_tens(1) centers_surf_tens centers_surf_tens(end) 72e-3];
counts_surf_tens = [0 0 counts_surf_tens 0 0];

% size figures histograms
size_fig = 0.6;%0.5;


fig_init_surf_tens = figure(53);
clf
set(gcf, 'Position', plot_position)
hold on
% bar(centers_surf_tens*1e3, counts_surf_tens/max(counts_surf_tens))
% if size_fig == 0.5 | size_fig==0.4
    plot(centers_surf_tens*1e3, counts_surf_tens/max(counts_surf_tens), 'linewidth', 1.5)%1.333)
% elseif size_fig == 1
    % plot(centers_surf_tens*1e3, counts_surf_tens/max(counts_surf_tens), 'linewidth', 2)
% end
% text(centers_surf_tens*1e3, counts_surf_tens, num2str(counts_surf_tens'),'vert','bottom','horiz','center');
xlabel('\sigma_0 (mN/m)')%initial surface tension (mN/m)')
ylabel('count (normalized)')
grid on
% title('initial surface tension')
set(gca, 'fontsize', 14)
xlim([0 72])
ylim([0 1])
yticks(0:0.25:1)
xticks(0:20:80)
set(gca,'box','on')

% find value of Segers' curve (will be the mean init. surf. tens. of the
% individual bubbles, since it is matched to that!)
% [~, strain_0_idx] = min(abs(strain_Segers-0));
% xline(surf_tens_par_combi(strain_0_idx)*1e3, '--',...
%     'color', [0.8500 0.3250 0.0980], 'linewidth', 1.5)

% change plot size, font size etc.
[fig_init_surf_tens, ~] = Plotting_size_settings(fig_init_surf_tens, gca, size_fig);

if saving_figs_paper_combined

    % %edit height (match the surface tension median and elasticity plots' height)
    % cm_to_points = 28.3464567;
    % fig_width = [8.5 17]*cm_to_points; % single, double column
    % % fig_init_surf_tens.Position(4) = fig_size*0.8*fig_width(1);
    % fig_init_surf_tens.Position(4) = fig_elas.Position(4);

    saving_path = [saving_paper_path{data_idx} 'initial_surface_tension'];

    saveas(fig_init_surf_tens, [saving_path, '.fig'])
    exportgraphics(fig_init_surf_tens, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_init_surf_tens, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_init_surf_tens, [saving_path, '.png'])
    saveas(fig_init_surf_tens, [saving_path, '.svg'])

end


if saving_figs
    saveas(gcf, [saving_path_fig...
        '/hist initial surface tension'], 'png')
    saveas(gcf, [saving_path_fig...
        '/hist initial surface tension'], 'fig')
end

% % initial surface tension as a function of error
% figure(15)
% set(gcf, 'Position', plot_position)
% hold on
% plot(initial_surface_tension_EM_array*1e3,...
%     error_EM_array, '.', 'markersize', markersize, 'color', c(1,:))
% xlabel('initial surface tension (mN/m)')
% ylabel('final error')
% grid on
% title({'error final iteration (slope + mean value)',...
%     'of error minimization vs initial surface tension'})
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/error final iteration (slope + mean value) '...
%         'of error minimization vs initial surface tension'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/error final iteration (slope + mean value) '...
%         'of error minimization vs initial surface tension'], 'fig')
% end


%% initial elasticity histogram
[counts_init_elas, centers_init_elas] =...
    hist(initial_elasticity_EM_array{data_idx}, ...
    round(1+3.322*log(length(initial_elasticity_EM_array{data_idx}))));

% if no data, make count 0
centers_init_elas = [0 centers_init_elas(1) centers_init_elas centers_init_elas(end) 1.2];
counts_init_elas = [0 0 counts_init_elas 0 0];

size_fig = 0.6;%0.5;


fig_init_elas = figure(55);
clf
set(gcf, 'Position', plot_position)
hold on
% bar(centers_surf_tens*1e3, counts_surf_tens/max(counts_surf_tens))
% if size_fig == 0.5 | size_fig==0.4
    plot(centers_init_elas, counts_init_elas/max(counts_init_elas), 'linewidth', 1.5)%1.333)
% elseif size_fig == 1
    % plot(centers_init_elas, counts_init_elas/max(counts_init_elas), 'linewidth', 2)
% end
% text(centers_surf_tens*1e3, counts_surf_tens, num2str(counts_surf_tens'),'vert','bottom','horiz','center');
xlabel('\chi_0 (N/m)')
ylabel('count (normalized)')
grid on
% title('initial elasticity')
set(gca, 'fontsize', 14)
% xlim([0 72])
ylim([0 1])
yticks(0:0.25:1)
xticks(0:0.3:1.2)
% xlim([0.001 inf])
set(gca,'box','on')

% find value of Segers' curve
[~, strain_0_idx] = min(abs(strain_Segers-0));
xline(elasticity_Segers(strain_0_idx), '--',...
    'color', [0.8500 0.3250 0.0980], 'linewidth', 1.5)

% change plot size, font size etc.
[fig_init_elas, ~] = Plotting_size_settings(fig_init_elas, gca, size_fig);

if saving_figs_paper_combined

    % %edit height (match the surface tension median and elasticity plots' height)
    % cm_to_points = 28.3464567;
    % fig_width = [8.5 17]*cm_to_points; % single, double column
    % % fig_init_elas.Position(4) = fig_size*0.8*fig_width(1);
    % fig_init_elas.Position(4) = fig_elas.Position(4);

    saving_path = [saving_paper_path{data_idx} 'initial_elasticity'];

    saveas(fig_init_elas, [saving_path, '.fig'])
    exportgraphics(fig_init_elas, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_init_elas, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_init_elas, [saving_path, '.png'])
    saveas(fig_init_elas, [saving_path, '.svg'])

end


if saving_figs
    saveas(gcf, [saving_path_fig...
        '/hist initial elasticity'], 'png')
    saveas(gcf, [saving_path_fig...
        '/hist initial elasticity'], 'fig')
end


%% Pressure amplitude

if S.ub(4) - S.lb(4) > 1
    [counts_pres_amp, centers_pres_amp] =...
        hist(pressure_amplitude_EM_array, ...
        round(1+3.322*log(length(pressure_amplitude_EM_array))));
    
    figure(12)
    clf
    set(gcf, 'Position', plot_position)
    hold on
    bar(centers_pres_amp*1e-3, counts_pres_amp)
    % text(centers_pres_amp*1e-3, counts_pres_amp, num2str(counts_pres_amp'),'vert','bottom','horiz','center');
    xlabel('pressure amplitude (kPa)')
    ylabel('count')
    grid on
    title('pressure amplitude')
    set(gca, 'fontsize', 14)
    xlim([S.lb(4) S.ub(4)])
    xline(P.Pa*1e-3, 'r-', 'linewidth', 2)
    
    if saving_figs
        saveas(gcf, [saving_path_fig...
            '/hist pressure amplitude'], 'png')
        saveas(gcf, [saving_path_fig...
            '/hist pressure amplitude'], 'fig')
    end
end


%% RMS error
% c=colororder;
% 
% figure(10)
% set(gcf, 'Position', plot_position)
% hold on
% plot(initial_radii_EM_array*1e6, RMSE_array, '.', 'markersize', markersize, 'color', c(1,:))
% grid on
% xlabel('initial radius from error minimization (µm)')
% ylabel('RMS error')
% title('Root-mean-square error of surface tension (EM and measured curve (Segers))')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/RMS error of surface tension (EM and measured curve (Segers)) vs R_0'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/RMS error of surface tension (EM and measured curve (Segers)) vs R_0'], 'fig')
% end

%  % RMSE as a function of measurement index
%  figure(2)
%  hold on
%  set(gcf, 'Position', plot_position)
%  plot(selIdx_new, RMSE_array(selIdx), '.', 'markersize', markersize, 'color', c(1,:))
%  xlabel('measurement index')
%  ylabel('RMS error')
%  title('Root-mean-square error of surface tension (EM and measured curve (Segers))')
%  grid on
%
%  if saving_figs
%      saveas(gcf, [saving_path_fig...
%          '/RMS error of surface tension (EM and measured curve (Segers)) vs index'], 'png')
%      saveas(gcf, [saving_path_fig...
%          '/RMS error of surface tension (EM and measured curve (Segers)) vs index'], 'fig')
%  end


%     % check if the EM error is good measure of fit to measured curve (Segers)
%     % linear fit
%     [x, wfmIndx] = sort(error_EM_full(~isnan(error_EM)));
%     y = RMSE(~isnan(error_EM_full));
%     y = y(wfmIndx);
%     RMSE_error_fit = polyfit(x', y',1);

% figure(3)
% set(gcf, 'Position', plot_position)
% hold on
% plot(error_EM_array, RMSE_array, '.', 'markersize', markersize, 'color', c(1,:))
% %     plot(x, polyval(RMSE_error_fit, x), 'linewidth', 1.2)
% grid on
% xlabel('error from EM')
% ylabel('RMS error')
% %     legend('data', ['linear fit: slope ' num2str(RMSE_error_fit(1))],...
% %         'location', 'best')
% title('error from error minimization vs RMS error')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/RMS error vs EM error'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/RMS error vs EM error'], 'fig')
% end



% 
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% plot(RMSE(selIdx{data_idx}), EC_ratio(selIdx{data_idx}), '.', 'markersize', markersize, 'color', c(1,:))
% grid on
% ylabel('E/C ratio')
% xlabel('RMS error')
% title('RMSE error vs E/C ratio')
% 


%% expansion/compression ratio distribution (to see if initial sruface
% tension distribution makes sense
[counts_EC, centers_EC] =...
    hist(EC_ratio(selIdx{data_idx}), ...
    round(1+3.322*log(length(EC_ratio(selIdx{data_idx})))));

size_fig = 0.6;

fig_expans_compres = figure(54);
clf
set(gcf, 'Position', plot_position)
hold on
% bar(centers_EC, counts_EC)
% if size_fig == 0.5
    plot(centers_EC, counts_EC/max(counts_EC), 'color', 'k', 'linewidth', 1.5)
% elseif size_fig == 1
%     plot(centers_EC, counts_EC/max(counts_EC), 'color', 'k', 'linewidth', 2)
% end
% text(centers_EC, counts_EC, num2str(counts_EC'),'vert','bottom','horiz','center');
xlabel('E/C ratio')
ylabel('count (normalized)')
grid on
title('E/C ratio')
xlim([0.5 1.2])
ylim([0 1])
yticks(0:0.25:1)
xline(1, 'r', 'linewidth', 1.5)
set(gca, 'fontsize', 14)

% change plot size, font size etc.
[fig_expans_compres, ~] = Plotting_size_settings(fig_expans_compres, gca, size_fig);

if saving_figs_paper_combined
   
    saving_path = [saving_paper_path{data_idx} 'expansion_compression_ratio'];

    saveas(fig_expans_compres, [saving_path, '.fig'])
    exportgraphics(fig_expans_compres, [saving_path, '.pdf'], 'ContentType','vector')
    exportgraphics(fig_expans_compres, [saving_path, '.eps'], 'ContentType','vector')
    exportgraphics(fig_expans_compres, [saving_path, '.png'])
    saveas(fig_expans_compres, [saving_path, '.svg'])

end


if saving_figs
    saveas(gcf, [saving_path_fig...
        '/hist EC ratio'], 'png')
    saveas(gcf, [saving_path_fig...
        '/hist EC ratio'], 'fig')
end


%% error
% error as function of measurement index
%  figure(6)
%  set(gcf, 'Position', plot_position)
%  hold on
%  plot(selIdx_new, error_EM_array(selIdx_new), '.', 'markersize', markersize, 'color', c(1,:))
%  title('error EM as a function of idx')
%  xlabel('measurement index')
%  ylabel('error EM as a function of index')
%  grid on
%
%  if saving_figs
%      saveas(gcf, [saving_path_fig...
%          '/error EM as a function of index'], 'png')
%      saveas(gcf, [saving_path_fig...
%          '/error EM as a function of index'], 'fig')
%  end


% % initial radius as a function of error
% figure(9)
% set(gcf, 'Position', plot_position)
% hold on
% plot(initial_radii_EM_array*1e6, error_EM_array, '.', 'markersize', markersize, 'color', c(1,:))
% xlabel('radius (µm)')
% ylabel('final error')
% grid on
% title('error final iteration (slope + mean value) of error minimization vs initial radius')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/error final iteration (slope + mean value) of error minimization vs initial radius'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/error final iteration (slope + mean value) of error minimization vs initial radius'], 'fig')
% end
% 
% 
% % error distribution (to compare to other bubble populations)
% [counts_error, centers_error] =...
%     hist(error_EM_min(selIdx{data_idx}), ...
%     round(1+3.322*log(length(error_EM_min(selIdx{data_idx})))));
% 
% figure(91)
% set(gcf, 'Position', plot_position)
% hold on
% bar(centers_error, counts_error)
% text(centers_error, counts_error, num2str(counts_error'),'vert','bottom','horiz','center');
% xlabel('error')
% ylabel('count')
% grid on
% title('Error after error minimization routine')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/hist error'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/hist error'], 'fig')
% end

% %% pressure amplitude
% [counts, centers] =...
%     hist(pressure_amplitude_EM(selIdx{data_idx})*1e-3, ...
%     round(1+3.322*log(length(pressure_amplitude_EM(selIdx{data_idx})))));
% 
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% bar(centers, counts)
% text(centers, counts, num2str(counts'),'vert','bottom','horiz','center');
% xlabel('pressure amplitude (kPa)')
% xlim([0.8*P.Pa*1e-3, 1.2*P.Pa*1e-3])
% ylabel('count')
% grid on
% title('pressure amplitude from error minimization routine')
% 
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/hist pressure amplitude'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/hist pressure amplitude'], 'fig')
% end
% 

%% other plots

% % initial radius as a function of TIME: time after start is relTimeIndx/15
% % seconds
% figure(10)
% set(gcf, 'Position', plot_position)
% hold on
% % scatter(selIdx_new, R0EstCor(selIdx_new)*1e6, markersize,...
% %     'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', 'none',...
% %     'MarkerFaceAlpha', 0.5);
% % plot(selIdx, initial_radii_EM(selIdx)*1e6, '.', 'markersize', markersize)
% % color code based on phase regime
% plot(NaN, NaN, '.', 'color', [0 0.4470 0.7410], 'markersize', markersize)
% plot(NaN, NaN, '.', 'color', [0.8500 0.3250 0.0980], 'markersize', markersize)
%
% for i = 1: length(selIdx{data_idx})
%     if delay_EM_array(i) <= -pi
%         plot(time_idx_array(i), initial_radii_EM_array(i)*1e6, '.',...
%             'color', [0 0.4470 0.7410],...
%             'markersize', markersize)
%     else
%         plot(time_idx_array(i), initial_radii_EM_array(i)*1e6, '.',...
%             'color', [0.8500 0.3250 0.0980],...
%             'markersize', markersize)
%     end
% end
% legend('error minimization result (phase <= -\pi)',...
%     'error minimization result (phase > -\pi)', 'location', 'best')
% xlabel('time of measurement (minutes after start)')
% ylabel('initial radius (µm)')
% grid on
% title('initial radius as a function of time')
%
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/initial radius as a function of time'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/initial radius as a function of time'], 'fig')
% end

% % phase as function of time
% figure(13)
% set(gcf, 'Position', plot_position)
% hold on
% plot(time_idx_array(selIdx{data_idx}), delay_EM_array/pi, '.', 'markersize', markersize, 'color', c(1,:))
% title('phase EM as a function of time')
% xlabel('time of measurement (minutes after start)')
% ylabel('phase/pi EM')
% grid on
%
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/phase over pi EM vs time'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/phase over pi EM vs time'], 'fig')
% end



% % initial surface tension as a function of time
% figure(16)
% set(gcf, 'Position', plot_position)
% hold on
% plot(time_idx_array(selIdx{data_idx}),...
%     initial_surface_tension_EM_array*1e3, '.', 'markersize', markersize, 'color', c(1,:))
% title('initial surface tension EM as a function of time')
% xlabel('time of measurement (minutes after start)')
% ylabel('initial surface tension (mN/m)')
% grid on
%
% if saving_figs
%     saveas(gcf, [saving_path_fig...
%         '/initial surface tension EM vs time'], 'png')
%     saveas(gcf, [saving_path_fig...
%         '/initial surface tension EM vs time'], 'fig')
% end

% initial surface tension vs R0
figure(17)
clf
set(gcf, 'Position', plot_position)
hold on
plot(initial_radii_EM_array{data_idx}*1e6,...
    initial_surface_tension_EM_array{data_idx}*1e3, '.', 'markersize', markersize)%, 'color', c(1,:))
title('initial surface tension vs. R_0')
xlabel('initial radius (µm)')
ylabel('initial surface tension (mN/m)')
grid on

if saving_figs
    saveas(gcf, [saving_path_fig...
        '/initial surface tension EM vs R0'], 'png')
    saveas(gcf, [saving_path_fig...
        '/initial surface tension EM vs R0'], 'fig')
end


end

%% (small bubbles low E/C?)
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% plot(initial_radii_EM_array*1e6, EC_ratio(selIdx{data_idx}), '.', 'markersize', markersize, 'color', c(1,:))
% title('initial surface tension vs. E/C ratio')
% ylabel('expansion compression ratio')
% xlabel('initial radius (µm)')
% grid on
% 
% %% RMSE vs strainVal: if strainVal is too high, smapling not good and thus high RMSE? or error from fitting?
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% plot(strainVal(selIdx{data_idx}), RMSE(selIdx{data_idx}), '.', 'markersize', markersize, 'color', c(1,:))
% grid on
% xlabel('strainVal')
% ylabel('RMS error')
% title('RMS error vs strainVal')
% 
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% plot(strainVal(selIdx{data_idx}), error_EM_array, '.', 'markersize', markersize, 'color', c(1,:))
% grid on
% xlabel('strainVal')
% ylabel('error')
% title('error vs strainVal')

% 
% %%
% figure()
% set(gcf, 'Position', plot_position)
% hold on
% plot(EC_ratio(selIdx{data_idx}), RMSE(selIdx{data_idx}), '.', 'markersize', markersize, 'color', c(1,:))
% grid on
% xlabel('E/C ratio')
% ylabel('RMS error')
% title('RMS error vs E/C ratio')
% 


%%
% figure()
% set(gcf, 'Position', plot_position)
% plot(R0est(selIdx{data_idx})*1e6, SNR(selIdx{data_idx}), '.', 'markersize', markersize, 'color', c(1,:))
% title('initial radius vs. SNR')
% ylabel('SNR')
% xlabel('initial radius (µm)')
% grid on
% 
% %%
% figure()
% set(gcf, 'Position', plot_position)
% plot(strainVal, SNR, '.', 'markersize', markersize, 'color', c(1,:))
% title('strainVal vs. SNR')
% ylabel('SNR')
% xlabel('strainVal')
% grid on

% Apply the R0, sigma(R0), driving pressure delay and pressure amplitude found from 
% error minimization (using the function lsqcurvefit)
% and calculate f, and elastic and viscous contributions, and from that determine
% the surface tension

% Charlotte Nawijn, University of Twente, 2023


function  [strain_vec, dstrain_vec, strain_grid, dstrain_grid,...
    strain_vec_crop, surf_tension, strain_input_crop, surf_tension_crop, surf_tension_interp, surf_tension_old,...
    P_VE_int, P_VE_int_old, elastic_contrib, elastic_contrib_old] ...
    = SS_func_surface_tension(error_min_R0, error_min_sig_0, error_min_delay, ...
    error_min_Pa, input_param, wfmIndx, S, P, D, data_name)


R0 = error_min_R0;     % in m
sig_0 = error_min_sig_0;    % in N/m
delay = error_min_delay*P.Fs;
Pa = error_min_Pa;


%% Calculate viscoelastic contrib (non-dimensionalized)

% The low-frequency transducer has its own phase shift, input.phase_LF,
% which shifts the signal within the envelope! see 'scope_traces_determine_phase_shift_20231023.m'. 
Pacc_2_temp = Pa*D.envelop_hanning.*sin(2*pi*P.fUS.*D.time + input_param.phase_LF);    % incident acoustic pressure (Pa)

% Then, the bubble position may vary, which we account for by shifting the entire signal
% by a shift 'delay', within the range -2*pi to 0. A negative shift means
% the transmitted signal is shifted to arrive sooner than 'expected'
Pacc_2_shifted = zeros(size(D.Pacc));

shift = round(-delay);     % shift in number of samples
if shift > 0  % shift to the right (later)
    Pacc_2_shifted(shift+1:end) = Pacc_2_temp(1:end-shift);
elseif shift <= 0  % shift to the left (earlier)
    Pacc_2_shifted(1:end+shift) = Pacc_2_temp(-shift+1:end);
end
Pacc_2 = Pacc_2_shifted;

if S.plot_all
    figure()
    plot(D.time*1e6, D.Pacc*1e-3, 'linewidth', 2)
    hold on
    plot(D.time*1e6, Pacc_2*1e-3)
    xlabel('time (µs)')
    ylabel('acoutic pressure (kPa)')
    title('transmitted pressure')
    legend('original', 'shifted')
    grid on
end

P_VE = -(1 + input_param.strain_filtered) .* input_param.d2strain_dt2 ...
    - 3/2 * (input_param.dstrain_dt).^2 ...
    + P.T^2/(R0^2 *P.rho) * (...
    (P.P0 + 2*sig_0/R0) * (1./(input_param.strain_filtered + 1)).^(3*P.kap) .*...
    (1 - 3*P.kap/P.cw * R0/P.T * input_param.dstrain_dt)...
    - P.P0 - Pacc_2);
P_VE_old = -(1 + input_param.strain_filtered) .* input_param.d2strain_dt2 - ...
    3/2 * (input_param.dstrain_dt).^2 + ...
    P.T^2/(P.R0^2 *P.rho) * (...
    (P.P0 + 2*P.sig_0/P.R0) *...
    (1./(input_param.strain_filtered + 1)).^(3*P.kap) .*...
    (1 - 3*P.kap/P.cw * P.R0/P.T * input_param.dstrain_dt) - ...
    P.P0 - D.Pacc);



if S.plot_all
    figure()
    plot3(input_param.strain_filtered, input_param.dstrain_dt, P_VE, '.', 'MarkerSize', 4)
    xlabel('strain dR/R0')
    ylabel('strain rate d(dR/R0)/P.dt')
    zlabel('viscoelastic contrib (non-dimensional)');
    title({'viscoelastic contrib'...
        , ['measurement index: ' num2str(wfmIndx)]})
    grid on

    figure()
    plot3(input_param.strain_filtered, input_param.dstrain_dt, P_VE_old, '.', 'MarkerSize', 4)
    xlabel('strain dR/R0')
    ylabel('strain rate d(dR/R0)/P.dt')
    zlabel('viscoelastic contrib (non-dimensional)');
    title({'original viscoelastic contrib'...
        , ['measurement index: ' num2str(wfmIndx)]})
    grid on
end

if S.plot_all
    figure()
    sgtitle({'viscoelastic contrib'...
        , ['measurement index: ' num2str(wfmIndx)]})
    subplot(1,2,1)
    plot(input_param.strain_filtered, P_VE, '.-', 'markersize', 3)
    xlabel('strain dR/R0')
    ylabel('viscoelastic contrib (non-dimensional)')
    grid on
    subplot(1,2,2)
    plot(input_param.dstrain_dt, P_VE, '.-', 'markersize', 3)
    xlabel('strain rate d(dR/R0)/P.dt')
    ylabel('viscoelastic contrib (non-dimensional)')
    grid on

    figure()
    sgtitle({'original viscoelastic contrib'...
        , ['measurement index: ' num2str(wfmIndx)]})
    subplot(1,2,1)
    plot(input_param.strain_filtered, P_VE_old, '.-', 'markersize', 3)
    xlabel('strain dR/R0')
    ylabel('viscoelastic contrib (non-dimensional)')
    grid on
    subplot(1,2,2)
    plot(input_param.dstrain_dt, P_VE_old, '.-', 'markersize', 3)
    xlabel('strain rate d(dR/R0)/P.dt')
    ylabel('viscoelastic contrib (non-dimensional)')
    grid on
end


%% put P_VE on grid of strain and strain rate
if S.mask_cutoff > 0 
    strain_vec = linspace(min(input_param.strain_filtered)...
        , max(input_param.strain_filtered), 200);        
    dstrain_vec = linspace(-max(abs(input_param.dstrain_dt))...
        , max(abs(input_param.dstrain_dt)), 200);       
else
    strain_vec = linspace(median(mink(input_param.strain_filtered, S.range_vec))...
        , median(maxk(input_param.strain_filtered, S.range_vec)), 200);       
    dstrain_vec = linspace(median(mink(input_param.dstrain_dt, S.range_vec))...
        , median(maxk(input_param.dstrain_dt, S.range_vec)), 200);       
end

[strain_grid, dstrain_grid] = meshgrid(strain_vec, dstrain_vec);


% interpolate scattered data onto grid
P_VE_int = griddata(input_param.strain_filtered * 1e2,...
    input_param.dstrain_dt, P_VE,...
    strain_grid * 1e2, dstrain_grid, 'linear');

P_VE_int_old =  griddata(input_param.strain_filtered * 1e2,...
    input_param.dstrain_dt, P_VE_old,...
    strain_grid * 1e2, dstrain_grid, 'linear');


if S.plot_all
    % 3D
    figure()
    surf(strain_grid, dstrain_grid, P_VE_int,'edgecolor','none')
    xlabel('radial strain (dR/R0)')
    ylabel({'radial strain rate','(d(dR/R0)/dt)'})
    zlabel({'viscoelastic contribution', '(nondimensional)'});
    title('total viscoelastic contribution')
    % title({'corrected viscoelastic contrib on grid'...
    %     , ['measurement index: ' num2str(wfmIndx)]})
    grid on
    set(gca, 'fontsize', 14)
    xlim([-0.15 0.15])

    if S.saving_figs
        saving_folder = [cd '/error minimization figures/' data_name '/for presentations'];

        if ~exist(saving_folder, 'dir')
            mkdir(saving_folder)
        end

        saveas(gcf, [saving_folder '/bubble_'...
            num2str(wfmIndx) '_viscoelastic_contrib'], 'png')
        saveas(gcf, [saving_folder '/bubble_'...
            num2str(wfmIndx) '_viscoelastic_contrib'], 'fig')
    end

    figure()
    subplot(1,2,1)
    surf(strain_grid, dstrain_grid, P_VE_int_old,'edgecolor','none')
    xlabel('strain (R-R_0)/R_0')
    ylabel('strain rate d((R-R_0)/R_0)/P.dt')
    zlabel({'function f','(non-dimensionalized)'});
    grid on
    view(360,0)
    set(gca, 'fontsize', 14)

    subplot(1,2,2)
    surf(strain_grid, dstrain_grid, P_VE_int_old,'edgecolor','none')
    xlabel('radial strain (R-R_0)/R_0')
    ylabel('radial strain rate d((R-R_0)/R_0)/dt')
    sgtitle({'total viscoelastic contribution'})% interpolated'})
    grid on
    view(90,0)
    set(gca, 'fontsize', 14)

    figure()
    surf(strain_grid, dstrain_grid, P_VE_int_old,'edgecolor','none')
    xlabel('strain dR/R0')
    ylabel('strain rate d(dR/R0)/P.dt')
    zlabel('viscoelastic contrib (non-dimensionalized)');
    title({'original viscoelastic contrib on grid'...
        , ['measurement index: ' num2str(wfmIndx)]})
    grid on
    set(gca, 'fontsize', 14)

    % side and topview
    figure()
    sgtitle({'corrected viscoelastic contrib on grid'...
        , ['measurement index: ' num2str(wfmIndx)]})
    subplot(1,2,1)
    surf(strain_grid, dstrain_grid, P_VE_int,'edgecolor','none')
    xlabel('strain (dR/R0)')
    ylabel('strain rate (d(dR/R0)/P.dt)')
    zlabel('viscoelastic contrib (non-dimensionalized)');
    grid on
    view(360,0)
    set(gca, 'fontsize', 14)
    subplot(1,2,2)
    surf(strain_grid, dstrain_grid, P_VE_int,'edgecolor','none')
    xlabel('strain (dR/R0)')
    ylabel('strain rate (d(dR/R0)/P.dt)')
    zlabel('viscoelastic contrib (non-dimensionalized)');
    grid on
    view(90,0)
    set(gca, 'fontsize', 14)
end




%% Elastic contribution from full viscoelastic contribution

if S.elastic_from_full_viscoelastic

    % find point where derivative of strain is closest to 0
    [~, PosV0] = min(abs(dstrain_vec));

    % Average around Rdot = 0 to find the elastic contribution
    elastic_contrib = mean(P_VE_int(PosV0-5:PosV0+5,:), 1);
    elastic_contrib_old = mean(P_VE_int_old(PosV0-5:PosV0+5,:), 1);

    % plot the total elastic contribution
    if S.plot_all
        figure()
        plot(strain_vec, elastic_contrib_old, 'linewidth', 1.25)
        hold on
        plot(strain_vec, elastic_contrib, 'linewidth', 1.25)
        grid on
        xlabel('strain dR/R0')
        ylabel('elastic contribution (N.D.)')
        title({'elastic contribution',...
            ['measurement index: ' num2str(wfmIndx)]})

        legend('data', 'data (corrected R0, sigma(R0), delay, pressure amplitude)',...
            'measured curve (Segers)')
    end


else
    %% Elastic contribution from envelope of strain
    [~, LOCS_lower, ~] = findpeaks(-input_param.strain_filtered);
    [~, LOCS_upper, ~] = findpeaks(input_param.strain_filtered);

    zero_strain_rate_idx = [LOCS_lower, LOCS_upper];

    elastic_contrib_temp = -(1 + input_param.strain_filtered(zero_strain_rate_idx)) .*...
        input_param.d2strain_dt2(zero_strain_rate_idx) ...
        + P.T^2/(R0^2 *P.rho) * (...
        (P.P0 + 2*sig_0/R0) * (1./(input_param.strain_filtered(zero_strain_rate_idx) + 1)).^(3*P.kap) ...
        - P.P0 - Pacc_2(zero_strain_rate_idx));

    [strain_zero_strain_rate_sorted, sort_idx] = ...
        sort(input_param.strain_filtered(zero_strain_rate_idx));

    % put on regular vector: strain_vec
    elastic_contrib = lininterp1(strain_zero_strain_rate_sorted, ...
        elastic_contrib_temp(sort_idx), strain_vec);

    % plot the total elastic contribution
    if S.plot_all
        figure()
        plot(strain_vec, elastic_contrib, 'linewidth', 1.25)
        grid on
        xlabel('strain dR/R0')
        ylabel('elastic contribution (N.D.)')
        title({'elastic contribution'...
            , ['measurement index: ' num2str(wfmIndx)]})
        legend('data (corrected R0, sigma(R0), delay and pressure amplitude)')
    end

    elastic_contrib_old_temp = -(1 + input_param.strain_filtered(zero_strain_rate_idx)) .*...
        input_param.d2strain_dt2(zero_strain_rate_idx) ...
        + P.T^2/(P.R0^2 *P.rho) * (...
        (P.P0 + 2*P.sig_0/P.R0) * (1./(input_param.strain_filtered(zero_strain_rate_idx) + 1)).^(3*P.kap) ...
        - P.P0 - D.Pacc(zero_strain_rate_idx));
    elastic_contrib_old = lininterp1(strain_zero_strain_rate_sorted, ...
        elastic_contrib_old_temp(sort_idx), strain_vec);
    


end

%% Calculate surface tension
% from non-dimensionalization: elastic contribution is
% P.T^2 / (R_0^2 * P.rho) * 2 *surf_tension_ND(R) / ((dR/R0+1))
surf_tension_ND = elastic_contrib .* (strain_vec + 1) / 2;
surf_tension_ND_old = elastic_contrib_old .* (strain_vec + 1) / 2;

% dimensionalize: surf_tension = nondimensionalized surface tension
% times P.rho*R_0^3 / P.T^2
surf_tension_temp = surf_tension_ND * ((P.rho * R0^3)/P.T^2);
surf_tension_old = surf_tension_ND_old * ((P.rho * R0^3)/P.T^2);


surf_tension= surf_tension_temp;
if S.mask_cutoff > 0
    surf_tension(1:round(S.mask_cutoff*length(surf_tension_temp))) = NaN;
    surf_tension(round((1-S.mask_cutoff)*length(surf_tension_temp)):end) = NaN;
end

% determine surface tension 
NaNindices = isnan(surf_tension);
strain_vec_crop = strain_vec(~NaNindices);

strain_input = min(input_param.strain_filtered): S.strain_input_step: max(input_param.strain_filtered);
if S.mask_cutoff > 0
    strain_input_crop = strain_input(round(S.mask_cutoff*length(strain_input)):...
        round((1-S.mask_cutoff)*length(strain_input)));
else
    strain_input_crop = strain_input;
end

surf_tension_interp = lininterp1(strain_vec_crop, surf_tension(~NaNindices), strain_input_crop);

surf_tension_crop = surf_tension(~NaNindices);


if S.plot_all
    figure()
    plot(strain_vec, elastic_contrib_old, 'linewidth', 1.25)
    hold on
    plot(strain_vec, elastic_contrib, 'linewidth', 1.25)
    grid on
    xlabel('strain dR/R0')
    ylabel('elastic contribution (N.D.)')
    title({'elastic contribution'...
        , ['measurement index: ' num2str(wfmIndx)]})
    legend('data', 'data (corrected R0, sigma(R0), delay, pressure amplitude)',...
        'data spline smoothed', 'measured curve (Segers)')
end



end


% This function calculates the surface tension by first determining
% viscoelastic contribution (elastic+viscous contribution of  the bubble shell)
% as a stress-strain analysis for a single bubble shell
% characterization, fitting the initial radius (R0), delay and amplitude of the
% acoustic pressure, and initial surface tension (sig_R0). Then the surface
% tension is determined, and the error compared to a parametrization is
% calculated, which is of DPPC:DPPE-PEG5k, 9:1 mol% bubbles, from
% Segers et al Soft Matter 2018 "High-precision acoustic measurements
% of the nonlinear dilatational elasticity of phospholipid coated
% monodisperse microbubbles"

% Charlotte Nawijn, University of Twente, 2023

function xout = SS_func_error_lsqnonlin(x, input_param, wfmIndx, S, P, D, F)

R0 = x(1)/input_param.weights(1);
sig_0 = x(2)/input_param.weights(2);           
delay = x(3)*P.Fs/input_param.weights(3);       
Pa = x(4)/input_param.weights(4); 

%% Obtain direct input parametrization surface tension (Segers)
if S.mask_cutoff > 0
    strain_input = min(input_param.strain_filtered):...
        S.strain_input_step:...
        max(input_param.strain_filtered);
elseif S.mask_cutoff == 0
    strain_input = median(mink(input_param.strain_filtered, S.range_vec)):...
        S.strain_input_step:...
        median(maxk(input_param.strain_filtered, S.range_vec));
end

R_input = (strain_input*R0 + R0)*1e6; % in µm
[~, sigma_par_input] = A0cor3_parallel(R_input, F.fit_surf_tens, sig_0, R0);

% shift input surface tension sigma_par_input to overlap with initial
% surface tension (caused by difference in bubble radius, but shape
% parametrization is the same for all sizes (see mail Tim 20/7/2021)
shift_sigma_um = R0*1e6 - R_input(find(sigma_par_input > sig_0, 1, 'first') -1);
shift_sigma = round(shift_sigma_um / diff(R_input(1:2)));
if shift_sigma >= 0
    sigma_par_input_shifted = zeros(size(sigma_par_input));
    sigma_par_input_shifted(shift_sigma+1: end)...
        = sigma_par_input(1: end-shift_sigma);      % shift to the right
elseif shift_sigma < 0
    sigma_par_input_shifted = zeros(size(sigma_par_input));
    sigma_par_input_shifted(1: shift_sigma+end)...
        = sigma_par_input(-shift_sigma+1: end);      % shift to the right
    sigma_par_input_shifted(end+shift_sigma:end) = sigma_par_input_shifted(end+shift_sigma-1);
else
    sigma_par_input_shifted = NaN(size(sigma_par_input));
end

if S.plot_all
    figure()
    plot(R_input, sigma_par_input*1e3)
    hold on
    plot(R_input, sigma_par_input_shifted*1e3)
    plot(R0*1e6, sig_0*1e3, 'r.', 'MarkerSize', 30)
    legend('initial surface tension input', 'shifted',...
        'location', 'northwest')
    grid on
end


%% Calculate viscoelastic contribution (non-dimensionalized)

% The low-frequency transducer has its own phase shift, input.phase_LF,
% which shifts the signal within the envelope! see 'scope_traces_determine_phase_shift_20231023.m'. 
Pacc_2_temp = Pa*D.envelop_hanning.*sin(2*pi*P.fUS.*D.time + input_param.phase_LF);    % incident acoustic pressure (Pa)

% Then, the bubble position may vary, which we account for by shifting the entire signal
% by a shift 'delay', within the range -2*pi to 0. A negative shift means
% the transmitted signal is shifted to arrive sooner than 'expected'
Pacc_2_shifted = zeros(size(D.Pacc));

shift = round(-delay);     % shift in number of samples
if shift > 0   % shift to the right (later)
    Pacc_2_shifted(shift+1:end) = Pacc_2_temp(1:end-shift);
elseif shift <= 0   % shift to the left (earlier)
    Pacc_2_shifted(1:end+shift) = Pacc_2_temp(-shift+1:end);
end
Pacc_2 = Pacc_2_shifted;

if S.plot_all
    figure()
    plot(D.time*1e6, Pacc_2_temp*1e-3, 'linewidth', 1.5)
    hold on
    plot(D.time*1e6, Pacc_2*1e-3, 'linewidth', 1)
    xlabel('time (µs)')
    ylabel('acoutic pressure (kPa)')
    title('(interpolated) transmitted pressure')
    % legend('original', 'shifted')
    grid on
end

% Viscoelastic contribution P_VE
P_VE = -(1 + input_param.strain_filtered) .* input_param.d2strain_dt2 ...
    - 3/2 * (input_param.dstrain_dt).^2 ...
    + P.T^2/(R0^2 *P.rho) * (...
    (P.P0 + 2*sig_0/R0) * (1./(input_param.strain_filtered + 1)).^(3*P.kap) .*...
    (1 - 3*P.kap/P.cw * R0/P.T * input_param.dstrain_dt)...
    - P.P0 - Pacc_2);

%% put viscoelastic contribution on grid of strain and strain rate
% define the vectors of strain and strain rate
if S.mask_cutoff > 0
    strain_vec = linspace(min(input_param.strain_filtered)...
        , max(input_param.strain_filtered), 200);           
    dstrain_vec = linspace(-max(abs(input_param.dstrain_dt))...
        , max(abs(input_param.dstrain_dt)), 200);          
elseif S.mask_cutoff == 0
    strain_vec = linspace(median(mink(input_param.strain_filtered, S.range_vec))...
        , median(maxk(input_param.strain_filtered, S.range_vec)), 200);        
    dstrain_vec = linspace(median(mink(input_param.dstrain_dt, S.range_vec))...
        , median(maxk(input_param.dstrain_dt, S.range_vec)), 200);
end

[strain_grid, dstrain_grid] = meshgrid(strain_vec, dstrain_vec);

P_VE_int = griddata(input_param.strain_filtered * 1e2, ...
    input_param.dstrain_dt, P_VE, ...
    strain_grid * 1e2, dstrain_grid, 'linear');

if S.plot_all
    figure()
    surf(strain_grid, dstrain_grid, P_VE_int,'edgecolor','none')
    xlabel('strain (dR/R0)')
    ylabel('strain rate (d(dR/R0)/dt)')
    zlabel('viscoelastic contribution (non-dimensionalized)');
    title({'viscoelastic contribution on grid from linear interpolation'...
        , ['measurement index: ' num2str(wfmIndx)]})
    grid on

    if S.saving_figs
        folder = ['error minimization figures/intermediate plots/bubble ' num2str(wfmIndx)];
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        saveas(gcf, [folder '/bubble_' num2str(wfmIndx) '_viscelas_original_interp.png'])
        saveas(gcf, [folder '/bubble_' num2str(wfmIndx) '_viscelas_original_interp.fig'])
    end

end



%% Elastic contribution from full viscoelastic contribution

if S.elastic_from_full_viscoelastic
    % find point where the derivative of strain is closest to 0
    [~, PosV0] = min(abs(dstrain_vec));
    
    % Average around Rdot = 0 to find the elastic contribution
    elastic_contrib = mean(P_VE_int(PosV0-5:PosV0+5,:), 1);
    
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
end

%% Calculate surface tension
% from non-dimensionalization: elastic contribution is
% T^2 / (R_0^2 * rho) * 2 *surf_tension_ND(R) / ((dR/R0+1))
surf_tension_ND = elastic_contrib .* (strain_vec + 1) / 2;

% dimensionalize: surf_tension = nondimensionalized surface tension
% times rho*R_0^3 / T^2
surf_tension_temp = surf_tension_ND * ((P.rho * R0^3)/P.T^2);

%10/11/21: take away outer parts surface tension
surf_tension= surf_tension_temp;
if S.mask_cutoff > 0
    surf_tension(1:round(S.mask_cutoff*length(surf_tension_temp))) = NaN;
    surf_tension(round((1-S.mask_cutoff)*length(surf_tension_temp)):end) = NaN;
end


% define surface tension to use in error minimization
NaNindices = isnan(surf_tension);
strain_crop = strain_vec(~NaNindices);


% expected surface tension (from parametrization Tim)
surf_tension_Segers = sigma_par_input_shifted;
if S.mask_cutoff > 0
    surf_tension_Segers_crop = surf_tension_Segers(round(S.mask_cutoff*length(surf_tension_Segers)):...
        round((1-S.mask_cutoff)*length(surf_tension_Segers)));

    strain_input_crop = strain_input(round(S.mask_cutoff*length(strain_input)):...
        round((1-S.mask_cutoff)*length(strain_input)));
else
    surf_tension_Segers_crop = surf_tension_Segers;
    strain_input_crop = strain_input;
end

surf_tension_interp = lininterp1(strain_crop, surf_tension(~NaNindices), strain_input_crop);
surf_tension_smooth = smooth(surf_tension_interp, S.smoothingN);
surf_tension_smooth = surf_tension_smooth';

if S.plot_all
    figure()
    plot(strain_vec, surf_tension*1e3)
    hold on
    plot(strain_input_crop, surf_tension_interp*1e3)
    plot(strain_input_crop, surf_tension_Segers_crop*1e3, 'k--')
    plot(0, sig_0*1e3, 'r.', 'MarkerSize', 30)
    xlabel('strain')
    ylabel('surface tension (mN/m)')
    grid on
    legend('data', 'interpolated', 'Segers', 'initial')

end



%% Determine error in slopes buckled and elastic regime

% strain range_EM from buckled end (as index from the input strain,
% used for the curve from Segers)
slope1_index_start = find(strain_input_crop >= strain_crop(1), 1, 'first'); %1;  % from the start of the cropped strain
slope1_index_end = min(slope1_index_start + round(S.range_EM/S.strain_input_step) -1, ...
    find(surf_tension_Segers_crop >= 0 ...
    & surf_tension_Segers_crop <= 0.0001, 1, 'last'));

% strain range_EM from ruptured end (as index from the input strain,
% used for the curve from Segers)
slope2_index_end = find(strain_input_crop <= strain_crop(end), 1, 'last');

slope2_index_start = max(slope2_index_end - round(S.range_EM/S.strain_input_step) + 1,...
    find(surf_tension_Segers_crop ==...
    max(surf_tension_Segers_crop), 1, 'first'));

if S.plot_all
    figure()
    plot(strain_input_crop, surf_tension_smooth*1e3)
    hold on
    plot(strain_input_crop, surf_tension_Segers_crop*1e3)
    xline(strain_input_crop(slope1_index_start),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope1_index_end),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope2_index_start),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope2_index_end),...
        'HandleVisibility', 'off')
    xlabel('strain')
    ylabel('surface tension (mN/m)')
    grid on
    title('regions over which error in surface tension is determined')
end

% define regimes over which to minimize: 
% 1 = buckled,
% 2= ruptured
regime1_y = surf_tension_smooth(slope1_index_start:slope1_index_end);   % in N/m
regime2_y = surf_tension_smooth(slope2_index_start:slope2_index_end);   % in N/m


if (length(regime1_y) <= 1) || (length(regime2_y) <= 1)
    xout = 1e5 * ones(1, 2*round(S.range_EM/S.strain_input_step));
else


len_desired =  round(S.range_EM/S.strain_input_step);  % desired length for fitting regimes

% if length(regime1_y) < round(S.range_EM/S.strain_input_step)  %round(radius_range_elastic/6)
%     regime1_y = [1e5 .* ones(1, length(input_param.regime1)-length(regime1_y)), regime1_y];
% end
% alternative: interpolate to get the required number of points
if length(regime1_y) < len_desired
    regime1_x_old = strain_input_crop(slope1_index_start:slope1_index_end);

    regime1_x_interp = linspace(regime1_x_old(1), regime1_x_old(end), round(S.range_EM/S.strain_input_step) );
    regime1_y = lininterp1  (regime1_x_old, regime1_y, regime1_x_interp);
end

% if length(regime2_y) < round(S.range_EM/S.strain_input_step) %round(radius_range_elastic/6)
%     regime2_y = [1e5 .* ones(1, length(input_param.regime2)-length(regime2_y)), regime2_y];
% end
% alternative: interpolate to get the required number of points
if length(regime2_y) < len_desired
    regime2_x_old = strain_input_crop(slope2_index_start:slope2_index_end);

    regime2_x_interp = linspace(regime2_x_old(1), regime2_x_old(end), round(S.range_EM/S.strain_input_step) );
    regime2_y = lininterp1(regime2_x_old, regime2_y, regime2_x_interp);
end


%% define output
xout = [regime1_y - zeros(size(regime1_y)), regime2_y - P.sig*ones(size(regime2_y))];


%% plot result
if S.attempts == 1 && S.plot_all
    color = get(gca, 'ColorOrder');
    
    figure(100)
    clf
    set(gca, 'fontsize', 12)
    xlim([-0.15 0.15])
    ylim([-20 100])
    hold on
    plot(strain_vec, surf_tension*1e3, 'linewidth', 2, 'color', color(2,:))
    plot(strain_input_crop, surf_tension_Segers_crop*1e3, 'k--', 'linewidth',2 )
    plot(0, sig_0*1e3, 'r.', 'markersize', 30)
    plot(strain_input_crop, surf_tension_interp*1e3, 'linewidth', 2, 'color', color(5,:))
    plot(strain_input_crop, surf_tension_smooth*1e3, 'linewidth', 2, 'color', color(5,:))
    hold off
    grid on
    xlabel('strain dR/R0')
    ylabel('surface tension (mN/m)')
    title('surface tension from fit')
    legend('data (after fitting)', 'measured curve (Segers)', 'initial surface tension', 'spline smoothed', 'location', 'northwest')
    xline(strain_input_crop(slope1_index_start),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope1_index_end),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope2_index_start),...
        'HandleVisibility', 'off')
    xline(strain_input_crop(slope2_index_end),...
        'HandleVisibility', 'off')
    pause(0.1)
end

end

end


%% lininterp1 function
% Version 1.3.0.0 (4.82 KB) by Jeffrey Wu 
% linear interpolation, given set of X and V values, and an x query
% assumes X values are in strictly increasing order
%
% Differences from matlab built-in :
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%       

function v = lininterp1(X, V, x)

if length(X) ~= length(V)
    error('X and V sizes do not match'); 
end

for n = 1: length(x)
    pindex = find((x(n) >= X), 1, 'last');
    index = find((x(n) <= X), 1, 'first');
    if isempty(pindex)
        % warning('interpolating before beginning');
        pindex = index;
        slope = 0;
    elseif isempty(index)
        % warning('interpolating after end');
        index = pindex;
        slope = 0;
    elseif pindex == index
        slope = 0;
    else
        Xp = X(pindex);
        slope = (x(n) - Xp) / (X(index) - Xp);
    end
    v(n) = V(pindex) * (1 - slope) + V(index) * slope;
end
end



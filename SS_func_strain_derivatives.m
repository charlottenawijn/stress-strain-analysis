% Determine the two time derivatives of the strain

% Charlotte Nawijn, University of Twente, 2023

function [strain_filtered, dstrain_dt, d2strain_dt2] ...
    = SS_func_strain_derivatives(strain_full, strain_filtered_temp,...
    f_fft_appended, wfmIndx, S, P, D)


k = 2*pi*f_fft_appended;
N = length(strain_full(wfmIndx,:));

% determine Fourier transform
FT_strain_filtered = fft(strain_filtered_temp);
FT_strain_filtered(round(N/2+1):end) = 0;        %make single-sided

if S.plot_all
    f_fft = linspace(0, 1/P.dt, length(strain_full(wfmIndx, :)));   %0: P.Fs/N: P.Fs/2;

    figure()
    plot(f_fft/1e6, abs(fft(strain_full(wfmIndx, :))))
    hold on
    % plot(f_fft_appended/1e6, abs(FT_strain_appended))
    plot(f_fft_appended/1e6, abs(FT_strain_filtered))
    xlabel('frequency (MHz)')
    xlim([0 5])
    legend('original radial strain', 'filtered radial strain' )
    % legend('original radial strain', 'appended radial strain'...
    %     , 'filtered radial strain' )
    grid on
end

% Rbar is the strain: dR/R0=(R-R0 / R0) in the .pdf file detailing the non-dimensionalization
dstrain_dt_temp = 2.*real(ifft(1i*k.*FT_strain_filtered));            % k is angular freq
dstrain_dt_temp = P.T * dstrain_dt_temp;    % non-dimensionalize (dRbar/dtbar)

d2strain_dt2_temp = 2.*real(ifft(-k.^2.*FT_strain_filtered));
d2strain_dt2_temp = P.T^2 * d2strain_dt2_temp;    % non-dimensionalize (dRbar/dtbar)

% cut off appended sections (to prevent edge effects)
strain_filtered_old = strain_filtered_temp(S.ApLg+1:end-S.ApLg);
strain_filtered_full = strain_filtered_old;

dstrain_dt_old = dstrain_dt_temp(S.ApLg+1:end-S.ApLg);
dstrain_dt_full = dstrain_dt_old;

d2strain_dt2_full = d2strain_dt2_temp(S.ApLg+1:end-S.ApLg);


% plot results
if S.plot_all
    figure()
    subplot(3,1,1)
    plot(D.time_full*1e6, strain_filtered_full)
    title('strain (filtered): dR/R0')
    ylabel('amplitude')
    xlim([D.time(1)*1e6 D.time(end)*1e6])
    hold on

    subplot(3,1,2)
    plot(D.time_full*1e6, dstrain_dt_full)
    title('first derivative')
    ylabel('amplitude')
    xlim([D.time(1)*1e6 D.time(end)*1e6])

    subplot(3,1,3)
    plot(D.time_full*1e6, d2strain_dt2_full)
    title('second derivative')
    ylabel('amplitude')
    xlabel('time (µs)')
    xlim([D.time(1)*1e6 D.time(end)*1e6])
end


%% Cropping the filtered strain and its derivatives
% envelope
strain_env = envelope(strain_filtered_full, 30, 'peak');

if S.plot_all
    figure()
    plot(D.time_full*1e6, strain_filtered_full)
    hold on
    plot(D.time_full*1e6, strain_env, 'linewidth', 1.25)
    grid on
    xlabel('time (µs)')
    ylabel('strain (filtered): dR/R0')
    title('filtered strain and envelope')
    legend('strain filtered', 'envelope')
end


% crop the strain, its derivatives, the acoustic pressure and time
% arrays
strain_filtered = strain_filtered_full;%(start_idx: end_idx);
dstrain_dt = dstrain_dt_full;%(start_idx: end_idx);
d2strain_dt2 = d2strain_dt2_full;%(start_idx: end_idx);



%% shift wrt driving pulse
env_strain = envelope(strain_filtered, round(P.T*P.Fs), 'peak');
env_strain = smooth(env_strain, 1000);
% 
% figure(3)
% clf
% plot(D.envelop_hanning_full)
% hold on
% plot(env_strain/max(env_strain))
% 
% [c, lags] = xcorr(D.envelop_hanning_full, env_strain, 'normalized');
% [val, shift_idx] = max(c);
% 
% % figure(2)
% % clf
% % hold on
% % stem(lags,c)
% % plot(lags(shift_idx), val, 'r', 'markersize', 20)
% % grid on
% 
% shift_samples = lags(shift_idx);
% shift_us = lags(shift_idx)/P.Fs*1e6;


%% plot
if ~S.calculate_on_server
    figure()
    p1 = plot(D.time_full*1e6, strain_filtered, 'linewidth', 1.25);
    p1.Color(4) = 1;
    hold on
    p2 = plot(D.time*1e6, D.Pacc/max(D.Pacc)*max(strain_filtered));
    p2.Color(4) = 0.5;
    legend('strain filtered cropped',  'driving pressure (a.u.)','location', 'northwest')
    grid on
    xlabel('time (µs)')
    ylabel('\epsilon')
    title({['measurement index: ' num2str(wfmIndx)], ...
        'filtered & cropped strain and driving pressure'})

    % plot(D.time_full*1e6, env_strain, 'linewidth', 2)


    if S.saving_figs
        saving_folder = ['error minimization figures/'];

        if S.section_analysis == 1
            saving_folder = [saving_folder 'updown/'];
        elseif S.section_analysis == 2
            saving_folder = [saving_folder 'up/'];
        elseif S.section_analysis ==3
            saving_folder = [saving_folder 'down/'];
        end

        saving_folder = [saving_folder 'range_EM_' num2str(S.range_EM) '_mask_cutoff_'...
            num2str(S.mask_cutoff) '/' S.data_name];

        if ~exist(saving_folder, 'dir')
            mkdir(saving_folder)
        end

        saveas(gcf, [saving_folder '/bubble_'...
            num2str(wfmIndx) '_strain'], 'png')
        saveas(gcf, [saving_folder '/bubble_'...
            num2str(wfmIndx) '_strain'], 'fig')
    end
end


end
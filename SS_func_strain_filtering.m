% Filtering strain around the fundamental and harmonic

% Charlotte Nawijn, University of Twente, 2023


function [strain_filtered_temp, f_fft_appended] = SS_func_strain_filtering(...
    strain_full, wfmIndx, S, P, D)


%% Frequency spectrum of dR/R0 (strain) -------------------------------
N = length(strain_full(wfmIndx,:));
FT_strain = fft(strain_full(wfmIndx,:));
FT_strain_temp = 2*abs(FT_strain ./ N);
FT_strain_temp(round(N/2+1):end) = 0;       % make single sided
FTlog_strain = 20.*log10(FT_strain_temp./max(FT_strain_temp));      % in dB

% frequency vector
f_fft = linspace(0, 1/P.dt, N);   %0: P.Fs/N: P.Fs/2;

if S.plot_all
    figure()
    plot(f_fft*1e-6, FTlog_strain)
    hold on
    title({'Frequency spectrum of full radial strain (dR/R0)'...
        , ['measurement index: ' num2str(wfmIndx)]})
    xlabel('frequency (MHz)')
    ylabel('amplitude (dB)')
    xlim([0 10])
    grid on
    xl_d = xline(P.fUS/1e6, 'r--', 'driving');
    xl_h = xline(P.fUS*2/1e6, 'r--', 'harmonic');
    xl_d.Alpha = 0.3;
    xl_h.Alpha = 0.3;
    ylim([-100 -0])
end


%% Filter
% append signal to beginning and end (to prevent edge effects)
strain1 = 2*strain_full(wfmIndx,1) - fliplr(strain_full(wfmIndx, 2:S.ApLg+1));
strain2 = 2*strain_full(wfmIndx, end) - fliplr(strain_full(wfmIndx,end-S.ApLg:end-1));
strain_appended = [strain1, strain_full(wfmIndx, :) strain2];

if S.plot_all
    figure()
    hold on
    plot([zeros(1,S.ApLg) strain_full(wfmIndx, :), zeros(1,S.ApLg)])
    plot(strain_appended)
    title('appended signal for filtering')
    xlabel('points')
    ylabel('radial strain dR/R_0')
    grid on
    xlim([0 2*S.ApLg])
    legend('original strain', 'appended strain')
end

% find Fourier transform
N_appended = length(strain_appended);
FT_strain_appended = fft(strain_appended);
FT_strain_appended_temp = 2*abs(FT_strain_appended ./ N_appended);
FT_strain_appended_temp(round(N_appended/2+1):end) = 0;       % make single sided
FTlog_strain_appended = 20.*log10(FT_strain_appended_temp./...
    max(FT_strain_appended_temp));      % in dB

% frequency vector
f_fft_appended = linspace(0, 1/P.dt, N_appended);   %0: P.Fs/N: P.Fs/2;

% find the indices of the driving and 2nd harmonic frequencies
p_f_driving = find(f_fft_appended <= P.fUS, 1, 'last');      % index of driving freq
p_f_harmonic = find(f_fft_appended <= 2*P.fUS, 1, 'last');   % index of harmonic freq

mask = zeros(size(f_fft_appended));

for h_idx = 1:2 % 7
    % make the mask around the i+1'th harmonic (i = 1 is fundamental)
    p_f_low(h_idx) = find(f_fft_appended <= (h_idx)*P.fUS - S.width_mask/2, 1, 'last');
    p_f_high(h_idx) = find(f_fft_appended <= (h_idx)*P.fUS + S.width_mask/2, 1, 'last');

    % create the mask
    mask(p_f_low(h_idx):p_f_high(h_idx)) = 1;

end
if ~S.filtering
    % if not filtering: AND DIVIDE STRAIN_FILTERED_TEMP BY 2! since mask is
    % not for single-sided spectrum now, or use strain_appended
    mask = ones(1, length(f_fft_appended));
end

% apply the mask
FT_strain_filtered = FT_strain_appended .* mask;
FTlog_strain_filtered = FTlog_strain_appended .* mask;

% plot unfiltered and filtered spectrum of the strain to show the mask
if S.plot_all
    figure()
    plot(f_fft_appended*1e-6, FTlog_strain_appended)
    hold on
    plot(f_fft_appended*1e-6, FTlog_strain_filtered)
    xlim([0 10])
    xlabel('frequency (MHz)')
    ylabel('amplitude (dB)')
    title({'Frequency spectrum of strain (dR/R0)'...
        , ['measurement index: ' num2str(wfmIndx)]})
    legend('unfiltered', 'filtered')
    grid on
    ylim([-100 -0])
end

%% Return to the time domain and compare the (un)filtered strain ------
strain_filtered_temp = 2.*real(ifft(FT_strain_filtered));

if ~S.filtering
    strain_filtered_temp = [zeros(size(strain1)), strain_full(wfmIndx, :), zeros(size(strain2))];
end

if S.plot_all
    figure()
    plot(D.time_full*1e6, strain_full(wfmIndx,:)) %plot(P.T, Strain(wfmIndx,:))
    hold on
    plot(D.time_full*1e6, strain_filtered_temp(S.ApLg+1: end-S.ApLg))
    legend('original radial strain', 'filtered radial strain')
    xlabel('time (µs)')
    ylabel('strain (dR/R0)')
    set(gca,'FontSize',14)
    title({'(filtered) radial strain'...
        , ['measurement index: ' num2str(wfmIndx)]})
    grid on


    figure()
    sgtitle(['measurement index: ' num2str(wfmIndx)])
    subplot(2,1,1)
    plot(D.time*1e6, D.Pacc*1e-3)
    hold on
    ylabel('P_{ac} (kPa)')
    set(gca,'FontSize', 14)
    grid on
    title('driving pressure', 'Fontsize', 12)
    xlim([D.time_full(1)*1e6 D.time_full(end)*1e6])
    subplot(2,1,2)      % strain dr/R0
    plot(D.time_full*1e6, strain_filtered_temp(S.ApLg+1: end-S.ApLg))%plot(t1,Strain(wfmIndx,coordsVec))
    hold on
    ylabel('dR/R0')
    set(gca,'FontSize', 14)
    grid on
    title('radial strain filtered', 'Fontsize', 12)
    xlabel('time (µs)')
    xlim([D.time_full(1)*1e6 D.time_full(end)*1e6])
end

end
% The surface tension from a parametrization is
% calculated, which is of DPPC:DPPE-PEG5k, 9:1 mol% bubbles, from
% Segers et al Soft Matter 2018 "High-precision acoustic measurements
% of the nonlinear dilatational elasticity of phospholipid coated
% monodisperse microbubbles"

% Charlotte Nawijn, University of Twente, 2023

% function [F.strain_input_crop, F.surf_tension_Segers_crop] ...
function [F] ...
    = SS_func_surface_tension_Segers(error_min_R0, error_min_sig_0, input_param, S, F)

R0 = error_min_R0;     % in m
sig_0 = error_min_sig_0;    % in N/m
% phase = error_min_phase;    

strain_input = min(input_param.strain_filtered): S.strain_input_step: max(input_param.strain_filtered);

R_input = (strain_input*R0 + R0)*1e6; % in µm
[~, sigma_par_input] = A0cor3_parallel(R_input, F.fit_surf_tens, sig_0, R0);


% shift input surface tension sigma_par_input to overlap with initial
% surface tension (caused by difference in bubble radius, but shape
% parametrization is the same for all sizes (see mail Tim 20/7/2021)
% idx_sig_0 = find(sigma_par_input > sig_0, 1, 'first');

% if isempty(idx_sig_0)
%     continue
% end
% if idx_sig_0 <= 1 || isnan(idx_sig_0)
%     continue
% end

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
end

surf_tension_Segers = sigma_par_input_shifted;

if S.mask_cutoff > 0
    F.surf_tension_Segers_crop = surf_tension_Segers(round(S.mask_cutoff*length(surf_tension_Segers)):...
        round((1-S.mask_cutoff)*length(surf_tension_Segers)));

    F.strain_input_crop = strain_input(round(S.mask_cutoff*length(strain_input)):...
        round((1-S.mask_cutoff)*length(strain_input)));
else
    F.surf_tension_Segers_crop = surf_tension_Segers;
    F.strain_input_crop = strain_input;
end

end
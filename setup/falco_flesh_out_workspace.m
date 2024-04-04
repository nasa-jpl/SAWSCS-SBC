% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Function to flesh out the workspace prior to wavefront estimation and control.

function [mp, out] = falco_flesh_out_workspace(mp)

mp = falco_set_optional_variables(mp); % Optional/hidden boolean flags and variables

mp = falco_set_spectral_properties(mp);
mp = falco_set_jacobian_weights(mp); % Zernike Modes and Subband Weighting of Control Jacobian

%--Pupil Masks
mp = falco_gen_chosen_pupil(mp);
mp = falco_gen_chosen_apodizer(mp);
mp = falco_gen_chosen_lyot_stop(mp);

%falco_plot_superposed_pupil_masks(mp); %--Visually inspect relative pupil mask alignment

% ******************* by erkin ******************
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_512pix.fits');
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_10mm_512pix_sharp_v2.fits');
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_10mm_512pix_soft.fits');
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_cir95_512pix_soft.fits');
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_cir98_512pix_soft.fits');
%lyo  = fitsread('/home/esidick/2022habex/dat/lyot_cir99_512pix_soft.fits');

%lyot = fitsread('/home/esidick/2022habex/dat_macos/lyo256_ratio098.fits');
% Should have used the above
%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_cir96_512pix_soft.fits'); %used this in dat_bb_256pix/

lyot = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/lyo_256pix.fits');


%load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/lyo_512pix_098D lyo
   %lyot = lyo;

%lyot = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/lyo512_ogap3_igap1_v2.fits');   
%lyot = fitsread('/home/esidick/2022habex/dat_macos/lyo256_ratio098.fits');

mp.P4.compact.mask = lyot; %(load your file here)
mp.P4.full.mask    = lyot; %(load your file here)

%--Focal planes
mp = falco_gen_FPM(mp); %% Generate FPM if necessary. If HLC, uses DM8 and DM9.
mp = falco_get_FPM_coordinates(mp); 
mp = falco_get_Fend_resolution(mp); 
mp = falco_configure_dark_hole_region(mp); %% Software Mask for Correction (corr) and Scoring (score)
mp = falco_set_spatial_weights(mp); %--Spatial weighting for control Jacobian. 
% mp.Fend.mask = ones(mp.Fend.Neta,mp.Fend.Nxi); %% Field Stop at Fend (as a software mask) (NOT INCLUDED YET)



%%--Wavefront sensor
mp = falco_configure_WFS(mp);

%--DM1 and DM2
mp = falco_configure_dm1_and_dm2(mp); %% Flesh out the dm1 and dm2 structures
mp = falco_gen_DM_stops(mp);
mp = falco_set_dm_surface_padding(mp); %% DM Surface Array Sizes for Angular Spectrum Propagation with FFTs

before_set_initial_Efields = 1

mp = falco_set_initial_Efields(mp);

before_get_PSF_norm_factor = 1

mp = falco_get_PSF_norm_factor(mp);
% mp = falco_gen_contrast_over_NI_map(mp); %--Contrast to Normalized Intensity Map Calculation (NOT INCLUDED YET)

after_get_PSF_norm_factor = 1


%figure(1), clf, imagesc(mp.P1.full.mask), axis xy image, colormap(gray), colorbar
%figure(2), clf, imagesc(mp.P3.full.mask), axis xy image, colormap(gray), colorbar
%figure(3), clf, imagesc(mp.P4.full.mask), axis xy image, colormap(gray), colorbar
disp('paused-workspace'), %pause




out = falco_init_storage_arrays(mp); %% Initialize Arrays to Store Performance History

%--Save the config file
fn_config = [mp.path.config mp.runLabel,'_config.mat'];
save(fn_config,'mp')
fprintf('Saved the config file: \t%s\n',fn_config)

fprintf('\nBeginning Trial %d of Series %d.\n',mp.TrialNum,mp.SeriesNum);

end %--END OF FUNCTION
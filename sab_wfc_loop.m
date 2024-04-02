% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

model_mode = 'wfc';
pr2_hex_vvc

load([falcodir 'macos/hex_dat1/efc3_trial3_opd1_gap8_256pix_bw20_27nm_2to12lamD_chg6_4lam_itr410'],'dm1','dm2','tx','cbx','fom','Im','bvalo')
    dm1_save = dm1;
    dm2_save = dm2;

load([falcodir 'macos/IFhex_opd/target_opd_amp_gap08_opd1_27nm_bw20_1500nm_post_efc'],'opd_target','amp_target','ce_target')
    mp.wfc.opd_target = opd_target;
    mp.wfc.amp_target = amp_target;
    mp.wfc.ce_target  = ce_target;

pupil_save = fitsread([falcodir 'macos/hex_maps/pupil_8mm_gap_8pix_256_matched.fits']);

load hex_maps/dat_dave_trial3_opdnm_3d opdnm_3d
opdnm0 = opdnm_3d(:,:,1);

pupil = pupil_save .* exp(1i*2*pi*opdnm0*1e-9/mp.lambda0);
mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;
[xm, ym] = meshgrid(xv);
rm = abs(xm + 1i * ym);
maskc = (rm > 2) & (rm <= 12);

if ~exist('G')
    load IFhex_opd/G_dm1_dm2_lim1e3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx
    mp.wfc.indx = indx;
    mp.wfc.km1  = km1;
    mp.wfc.km2  = km2;
    mp.wfc.G    = G;
end

% Load pre-calculated OPDs from Redding
trial_name = 'trial8'; trial = read_trial(trial_name); % possible options for trial_name: trial4, trial6
for jj = 1:size(trial,2)
  delta_opd = trial(jj).delta; % relevant variable name is trial3.delta(1:10).dopd_drift and .dopd_final
  ns = size(delta_opd,2);

  for ii = 1:ns

    opdnm = opdnm0 + pad(zernike_remove(1e3*delta_opd(ii).dopd_final,[1:3]),length(pupil));
    figure(3),show_opd(opdnm,sprintf('OPD drift case %0.2d',ii),'nm'),drawnow

    pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    dm1 = dm1_save; % + dm1x_3d(:,:,ii);
    dm2 = dm2_save; % + dm2x_3d(:,:,ii);

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

    mp.wfc.fn = sprintf('%smacos/IFhex_opd/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_%s_case%i_%0.2i', ...
	                     mp.falcodir, trial_name, jj, ii);

% *******************************
    flag = 6
    mp.wfc.model_id = flag; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp, 8 = target_amp_opd_eval
% *******************************
    Im = falco_get_summed_image(mp);
    %disp('paused'), pause
    
    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    cx = [-12 -8];
    fs = 14;


    figure(2), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        title(sprintf('post-wfc-case%i: Cb = %0.3g ',  ii, cb0));
        axis(12*[-1 1 -1 1])
        drawnow
   %print('-dpng', sprintf('../Figs/fig_psf_case%i',ii))

    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')
    rmsx = rms_opd;
    x = [1:size(rmsx,1)]-1;
    x = x';
    figure(1), clf, semilogy(x, rmsx(:,1)*1e3, 'ro-', x, dmrms(:,1)*1e3, 'bs-'), grid
    legend(sprintf('wfe-rms: last = %0.1fpm', rmsx(end,1)*1e3), 'dm1-rms [pm]', 'Location', 'Northeast')
    drawnow
  end % ii - loop
end % jj - loop

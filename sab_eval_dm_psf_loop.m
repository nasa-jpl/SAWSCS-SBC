% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

model_mode = 'eval';
pr2_hex_vvc

ifPlot = 0; ifPrint = 0; ifSave = 1; ifSavePSF = 0;

load([falcodir 'macos/hex_dat1/efc3_trial3_opd1_gap8_256pix_bw20_27nm_2to12lamD_chg6_4lam_itr410'],'dm1','dm2','tx','cbx','fom','Im','bvalo')
dm1_save = dm1;
dm2_save = dm2;

load([falcodir 'macos/IFhex_opd/target_opd_amp_gap08_opd1_27nm_bw20_1500nm_post_efc'],'opd_target','amp_target','ce_target')
    mp.wfc.opd_target = opd_target;
    mp.wfc.amp_target = amp_target;
    mp.wfc.ce_target  = ce_target;

pupil_save = fitsread([falcodir 'macos/hex_maps/pupil_8mm_gap_8pix_256_matched.fits']);

% Below just uses an example opd map to create the pupil mask
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
tightenAxis = abs(xv) <= 12;  xvt = xv(tightenAxis);

idd = 3000

%load hex_maps/trial3_opd1_final opdx %SAB2024
trial_name = 'trial13'; trial = read_trial(trial_name);

nt = size(trial,2);
ns = size(trial(1).delta,2);
cbxrb = zeros(ns, nt); cbxi = zeros(ns, nt); cbxo = zeros(ns, nt);
wfxrb = zeros(ns, nt); dwfxrb = zeros(ns, nt); wfxi = zeros(ns, nt); dwfxi = zeros(ns, nt);
offset = 0;
if ifSavePSF, offset = 5; end
clear errFac cbrbmean cbimean cbomean wfrbmean dwfrbmean wfimean dwfimean psfrbtot psfitot psfotot
before_met = 0; before_met_str = 'dopd_drift'; after_met = 0; after_met_str = 'dopd_final';
for jj = 1+offset:nt %6 %1:nt
  delta_opd = trial(jj).delta; % relevant variable name is trial3.delta(1:10).dopd_drift and .dopd_final
  if isfield(trial(1).delta,  after_met_str),  after_met = 1; end
  if isfield(trial(1).delta, before_met_str), before_met = 1; end
  if isfield(trial(jj),'errFac')
    errFac(jj) = trial(jj).errFac;
  elseif isfield(trial(jj),'sig_dz')
    errFac(jj) = trial(jj).sig_dz;
  end

  flag_movie = 0;
  if flag_movie
      video = VideoWriter('dat_dm_err1/movie_psf_1000pm.avi'); %create the video object
      video.FrameRate = 1;
      open(video); %open the file for writing
  end

  for ii = 1:ns %7 %1:ns
    if before_met
      % Compute pre-WFC Contrast
      opdrb = 1e3*zernike_remove(delta_opd(ii).dopd_drift,[1:3]);
      opdnm = opdnm0 + pad(opdrb,length(pupil));
      wfxrb(ii,jj) = rms2(nonzeros(opdnm));
      dwfxrb(ii,jj) = rms2(nonzeros(opdrb));
      pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

      if ifPlot, smallerLim = min(max(opdrb(:)), -min(opdrb(:))); clim2 = [-smallerLim smallerLim];
        figure(11),show_opd(opdrb,sprintf('Drift WFE, ErrFac = %d',errFac(jj)),'nm',clim2)
        if ifPrint,print('-dpng',sprintf('../Figs/DriftWFE%0.2i',ii)), end, end
      mp.P1.compact.mask = pupil; %(load your file here)
      mp.P1.full.mask    = pupil; %(load your file here)

      mp.dm1.V = dm1_save;
      mp.dm2.V = dm2_save;

      % *******************************
      flag = 1
      mp.wfc.model_id = flag; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp, 8 = target_amp_opd_eval
      % *******************************
      Im = falco_get_summed_image(mp);
      psfrb = Im .* maskc;
      if ifSavePSF, psfrbtot(:,:,ii,jj-offset) = psfrb; end
      cbrb = mean(nonzeros(psfrb(:)))
      cbxrb(ii,jj) = cbrb;

      if ifPlot, figure(1), clf % SAB We added these lines because having trouble using sub_plot_rb_psf
	clim = [-11.4 -6];
	imagesc(xvt,xvt,log10(pad(psfrb,sum(tightenAxis))),clim); axis xy image, colormap(jet); colorbar
        title([sprintf('pre-wfc NI(%0.2i): Cb = %0.3g',  ii, cbrb) ]), drawnow
        if ifPrint, print('-dpng',sprintf('../Figs/preNI%0.2i',ii)), end
      end
    end % if before_met
% ===========================================================================================================
    if after_met
      % Compute post-wfc Contrast
      opdcor = 1e3*zernike_remove(delta_opd(ii).dopd_final,[1:3]);
      opdnm = opdnm0 + pad(opdcor,length(pupil));
      wfxi(ii,jj) = rms2(nonzeros(opdnm));
      dwfxi(ii,jj) = rms2(nonzeros(opdcor));
      pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

      if ifPlot, figure(12),show_opd(opdcor,sprintf('Metrology-corrected WFE, ErrFac = %d',errFac(jj)),'nm',clim2)
        if ifPrint,print('-dpng',sprintf('../Figs/MetCorWFE%0.2i',ii)), end, end
      mp.P1.compact.mask = pupil; %(load your file here)
      mp.P1.full.mask    = pupil; %(load your file here)

      Im = falco_get_summed_image(mp);
      psfi = Im .* maskc; 
      if ifSavePSF, psfitot(:,:,ii,jj) = psfi; end
      cbi = mean(nonzeros(psfi(:)))
      cbxi(ii,jj) = cbi;

      if ifPlot, figure(2), clf % SAB We added these lines because having trouble using sub_plot_rb_psf
        imagesc(xvt,xvt,log10(pad(psfi,sum(tightenAxis))),clim); axis xy image, colormap(jet); colorbar
        title([sprintf('post-wfc NI(%0.2i): Cb = %0.3g',  ii, cbi) ]), drawnow
	if ifPrint, print('-dpng',sprintf('../Figs/postNI%0.2i',ii)), end
      end
    end % if after_met
% ===========================================================================================================
    % -------------------------------------------------
    mp.wfc.fn = sprintf('%smacos/IFhex_opd/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_%s_case%i_%0.2i', ...
	                     mp.falcodir, trial_name, jj, ii);
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdi', 'ampi', 'opdx', 'ampx', 'rms_opd', 'rms_amp','dmrms')

    mp.dm1.V = dm1x; %(:,:,end); % dm1_save + dm1x_3d(:,:,ii);
    mp.dm2.V = dm2x; %(:,:,end);

    Im = falco_get_summed_image(mp);
    psfo = Im .* maskc; 
    if ifSavePSF, psfotot(:,:,ii,jj) = psfo; end
    cbo = mean(nonzeros(psfo(:)))
    cbxo(ii,jj) = cbo;

    if ifPlot, figure(3), clf % SAB We added these lines because having trouble using sub_plot_rb_psf
        imagesc(xvt,xvt,log10(pad(psfo,sum(tightenAxis))),clim); axis xy image, colormap(jet); colorbar
        title([sprintf('dm1&dm2 control NI(%0.2i): Cb = %0.3g',  ii, cbo) ]), drawnow
	if ifPrint, print('-dpng',sprintf('../Figs/afterNI%0.2i',ii)), end
    end
    % ------------------------------
        %sub_plot_rb_psf, drawnow
    % ------------------------------

    if flag_movie
        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(video,frame);
    end
    if 0
        figure(2), clf
        imagesc(xv,xv,log10(psfi)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        title([sprintf('case-%i: Cb = %0.3g',  ii, cb0) ]), drawnow
    end
  end % ii - loop
  if flag_movie, close(video); end

  cbrbmean(jj) = mean(cbxrb(:,jj)); cbimean(jj) = mean(cbxi(:,jj)); cbomean(jj) = mean(cbxo(:,jj));
  wfrbmean(jj) = mean(wfxrb(:,jj)); dwfrbmean(jj) = mean(dwfxrb(:,jj));
  wfimean(jj) = mean(wfxi(:,jj)); dwfimean(jj) = mean(dwfxi(:,jj));

  x = [1:ns];
  figure(4), clf, semilogy(x, cbxrb(:,jj), 'ro-', x, cbxi(:,jj), 'bo-', x, cbxo(:,jj), 'gs-'), grid
  legend('pre-wfc','post-wfc','dm1 + dm2 wfc'), axis([1 10 8e-11 5e-7]) % axis([1 10 9e-11 9.6e-11])
  title(sprintf('wfc w/ %s = 1500nm, trial %d','\lambda',jj))
  ylabel('mean broadband contrast, Cb'),xlabel('wfc case #')

end % jj - loop

figure(5)
yyaxis left, loglog(errFac,dwfrbmean,'k*-',errFac,dwfimean,'k*-.')
ylabel('WFE (nm)')

yyaxis right, loglog(errFac,cbrbmean,'ro-',errFac,cbimean,'bo-',errFac,cbomean,'go-'),grid on
title('wfc w/ \lambda = 1500nm')
ylabel('mean broadband contrast, C_b'),xlabel('Metrology error factor')

legend('Drift WFE','Metrology-corrected WFE','pre-wfc C_b','post-wfc C_b','dm1 + dm2 wfc C_b')

if ifSave
  if ifSavePSF
    save(sprintf('IFhex_opd/%s_contrast_psf',trial_name),'nt','ns','errFac','cbxrb','cbxi','cbxo','cbrbmean','cbimean', ...
        'cbomean','wfxrb','dwfxrb','wfxi','dwfxi','wfrbmean','dwfrbmean','wfimean','dwfimean','opdnm0', ...
	'psfrbtot','psfitot','psfotot')
  else
    save(sprintf('IFhex_opd/%s_contrast',trial_name),'nt','ns','errFac','cbxrb','cbxi','cbxo','cbrbmean','cbimean', ...
	'cbomean','wfxrb','dwfxrb','wfxi','dwfxi','wfrbmean','dwfrbmean','wfimean','dwfimean','opdnm0')
  end
end

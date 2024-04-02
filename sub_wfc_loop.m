% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

    idd = 1
    %disp('paused'), pause

model_mode = 'wfc';
pr2_hex_vvc

if 1


%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
load([falcodir 'macos/hex_dat1/efc3_trial3_opd1_gap8_256pix_bw20_27nm_2to12lamD_chg6_4lam_itr410'],'dm1','dm2','tx','cbx','fom','Im','bvalo')
    dm1_save = dm1;
    dm2_save = dm2;

    tc = fom(end)*10;

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
load([falcodir 'macos/IFhex_opd/target_opd_amp_gap08_opd1_27nm_bw20_1500nm_post_efc'],'opd_target','amp_target','ce_target')
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc_no_pm opd_target amp_target ce_target
    mp.wfc.opd_target = opd_target;
    mp.wfc.amp_target = amp_target;
    mp.wfc.ce_target  = ce_target;
%disp('paused'), pause

%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
%opdrb = opd_rb_drift_100pm;
%rmsrb = stat2d(opdrb);

if 0 % Make RB drift error SAB2024
   %load dat_rb_wfe/opdrb_100pm opdrb_nm
   load dat_rb_wfe/opdrb_1000pm opdrb_nm
   ns = size(opdrb_nm, 3);
   opdrb = pad(opdrb_nm(:,:,1), 256);
   rmsrb = stat2d(opdrb);

   % ********************************
   fac_rb = 300;  % 1, 3, 10, 30, 100 300, 1000, 3000 pm
   % ********************************
   opdrb = opdrb * fac_rb / rmsrb(1); %[nm]
    %figure(1), clf, imagesc(opdrb); axis xy image, colormap(jet), colorbar, return
end
%pupil_save = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
pupil_save = fitsread([falcodir 'macos/hex_maps/pupil_8mm_gap_8pix_256_matched.fits']);
%pupil_save = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');
%load dat_dm_err1/opd10_erkin dopd_c_nm
%load dat_dm_err1/opd30_erkin dopd_c_nm

%load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
load hex_maps/trial3_opd1_final opdx %SAB2024
opdx_fromDave = opdx;
ns = size(opdx_fromDave, 3);
opdnm0 = 1e-3*opdx_fromDave(:,:,1);

%load ~esidick/2022habex/dat_macos/opdnom_256pix_macos opdnm

%dopd_c_nm = opdnm;

%    opdnm = dopd_c_nm + 0 * opd_rb_drift_100pm * 1e-3;
    opdnm = opdnm0; % + opdrb * 1e-3;

    pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

%%

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);


if 0 % ====================================================================

    flag = 3
    mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, 4 = dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

    %disp('paused'), pause
    id1 = 1
    %ii = 1;
    %mp.wfc.fn = sprintf('~esidick/Afalco/falco20200916/macos/dat_rb_wfe/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_opdrb_100pm_run%i', ii);
    %load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')

id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

Im = falco_get_summed_image(mp);
    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));
return
end

if 0

    fs = 14;
    cx = [-12 -8];
    cx1 = [-7 -4];

    figure(1), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title(['pre-wfc: 30nm-20%-12\lambda, rb-drift = 100pm ']);
        xlabel([sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%']);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psf1
   %return
    

    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp); %<***** thput ***

%end
    
figure(2), clf
    imagesc(xv,xv,log10(ImSimOffaxis)); axis xy image, grid
    %hold on, plot(xc, yc, 'y', xc1, yc1, 'r', 'Linewidth', 2), hold off
    parms = fun_newaxes(fs, 0, 0);  
    title([sprintf('input-30nm: with flattening, Tc = %0.2f', thput*100) '%']), 
print -dpng ../Figs/fig_psf3

return
end % ====================================================================
%%

if 1
    a= 20
    %load IFopd/G_dm1_lim5_m1 G km1
    %load IFopd/G_dm2_lim5_m1 G km2
    %load IFopd/G_dm1_dm2_lim5_m2_1064nm G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm km1 km2 indx
    %load IFopd/dwddm_dm1_nm_over_nm_full indx maskopd
    %load IFopd/G_dm1_dm2_lim5_m2_355nm_30nm_bw20_flat G km1 km2 indx
    %load IFopd/G_dm1_lim3_m2_355nm_30nm_bw20_gap08 G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim3_m3_355nm_30nm_bw20_gap08 G km1 km2 indx
    %load IFopd/G_dm1_lim3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim3_m3_1500nm_30nm_bw20_gap08_no_ota G km1 km2 indx
    load IFhex_opd/G_dm1_dm2_lim1e3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx

%else

    mp.wfc.indx = indx;
    mp.wfc.km1  = km1;
    mp.wfc.km2  = km2;
    mp.wfc.G    = G;
end

%return

if 0    
    aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-opd: rms = %0.2f, pv = %0.1fnm', rms0));
    % figure(13),show_opd(opd_target,'opd: evaluated with \lambda = 1500nm','nm',cx)

    aa = amp_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-amp: rms = %0.2f, pv = %0.1f', rms0));
    % figure(14),show_opd(amp_target,'amp: evaluated with \lambda = 1500nm','')
end

id2 = 100

if 0
    ma = 64;
    del = 1; % nm

    for jj = 1:10
        dm1x_3d(:,:,jj) = del * randn(ma,ma);
        dm2x_3d(:,:,jj) = del * randn(ma,ma);
    end
    %dm1x = 2*(0.5 - rand(ma,ma));
    %dm2x = 2*(0.5 - rand(ma,ma));
    %dm1x(find(dm1x >  del)) = del;
    %dm1x(find(dm1x < -del)) =-del;
    %dm2x(find(dm2x >  del)) = del;
    %dm2x(find(dm2x < -del)) =-del;
    %save dat_dm_err1/dm_err_1nm dm1x dm2x
    save dat_dm_err/dm_err_100pm dm1x_3d dm2x_3d
%else
    %load dat_dm_err1/dm_err_1nm dm1x dm2x
    %load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
    load dat_dm_err2/dm_err_1000pm dm1x_3d dm2x_3d
end

%dmxx = [dm1x_3d(:,:,1) dm2x_3d(:,:,1)];
%rmso = stat2d(dmxx);

%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_150pm_itr20'];
    idd = 10
    %disp('paused'), pause

% ******************************************
%ff = 1 / rmso(1)
% ******************************************
idd = 10

for ii = 1:ns

    %opdrb = pad(opdrb_nm(:,:,ii), 256);

    opdnm = 1e-3*opdx_fromDave(:,:,ii); % + opdrb;

    pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    dm1 = dm1_save; % + dm1x_3d(:,:,ii);
    dm2 = dm2_save; % + dm2x_3d(:,:,ii);

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

    %mp.wfc.fn = sprintf('~esidick/Afalco/falco20200916/macos/dat_rb_wfe/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_opdrb_1000pm_run%i', ii);
    %mp.wfc.fn = sprintf('~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_1000pm_run%i', ii);
    mp.wfc.fn = sprintf('%smacos/IFhex_opd/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_1000pm_case%0.2i', mp.falcodir,ii);

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
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title(sprintf('pre-efc-10nm: \\sigma = 1nm, Cb = %0.3g ',  cb0));
        title(sprintf('post-wfc-case%i: Cb = %0.3g ',  ii, cb0));
        axis(12*[-1 1 -1 1])
        drawnow
   %print('-dpng', sprintf('../Figs/fig_psf_case%i',ii))

    %load(mp.wfc.fn, 'dm1x', 'opdx', 'rmsx', 'dmrms')
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')
    rmsx = rms_opd;
    x = [1:size(rmsx,1)]-1;
    x = x';
    figure(1), clf, semilogy(x, rmsx(:,1)*1e3, 'ro-', x, dmrms(:,1)*1e3, 'bs-'), grid
    legend(sprintf('wfe-rms: last = %0.1fpm', rmsx(end,1)*1e3), 'dm1-rms [pm]', 'Location', 'Northeast')
    drawnow

end % ii - loop

return
end % -----------------------------------------------
%%








    %idd = 500
    %disp('paused'), pause

if 1


%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_10nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_20nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw05_30nm_2to12lamD_chg6_4lam_itr380 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw10_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw15_30nm_2to12lamD_chg6_4lam_itr370 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 

%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_00nm_2to12lamD_chg6_4lam_itr210 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw05_30nm_2to12lamD_chg6_3lam_itr330 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw10_30nm_2to12lamD_chg6_4lam_itr330 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw15_30nm_2to12lamD_chg6_4lam_itr330 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr300 dm1 dm2 tx cbx fom Im beta_value 
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 

    dm1_save = dm1;
    dm2_save = dm2;

    tc = fom(end)*10;

load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
opdrb = opd_rb_drift_100pm;
rmsrb = stat2d(opdrb);

%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_1500nm_m3_opdrb_100pm_itr20'];

% ********************************
fac_rb = 3000;  % 1, 3, 10, 30, 100 300, 1000, 3000 pm
% ********************************
opdrb = opdrb * fac_rb / rmsrb(1); %[nm]

    %figure(1), clf, imagesc(opdrb); axis xy image, colormap(jet), colorbar, return

%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');
%load dat_dm_err1/opd10_erkin dopd_c_nm
%load dat_dm_err1/opd30_erkin dopd_c_nm

load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
opdnm0 = opdnm;

%load ~esidick/2022habex/dat_macos/opdnom_256pix_macos opdnm

%dopd_c_nm = opdnm;

%    opdnm = dopd_c_nm + 0 * opd_rb_drift_100pm * 1e-3;
%    opdnm = opdnm0 + opdrb * 1e-3;

    pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

%%

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

flag = 4;

if 0 % ====================================================================

    mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, 4 = dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
    disp('paused'), pause

id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

Im = falco_get_summed_image(mp);
    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

%end

    fs = 14;
    cx = [-12 -8];
    cx1 = [-7 -4];

    figure(1), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title(['post-30nm-20%-12\lambda, ' sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%']);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psf1
   return
    

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
    a= 2
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
    load IFopd/G_dm1_dm2_lim3_m3_1500nm_30nm_bw20_gap08 G km1 km2 indx

%else

    mp.wfc.indx = indx;
    mp.wfc.km1  = km1;
    mp.wfc.km2  = km2;
    mp.wfc.G    = G;
end

if 0    
    aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-opd: rms = %0.2f, pv = %0.1fnm', rms0));

    aa = amp_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-amp: rms = %0.2f, pv = %0.1f', rms0));
end


if 1
    ma = 64;
    del = 0.01; % nm

    for jj = 1:10
        
    dm1x_3d(:,:,jj) = del * randn(ma,ma);
    dm2x_3d(:,:,jj) = del * randn(ma,ma);
    %dm1x = 2*(0.5 - rand(ma,ma));
    %dm2x = 2*(0.5 - rand(ma,ma));
    %dm1x(find(dm1x >  del)) = del;
    %dm1x(find(dm1x < -del)) =-del;
    %dm2x(find(dm2x >  del)) = del;
    %dm2x(find(dm2x < -del)) =-del;
    %save dat_dm_err1/dm_err_1nm dm1x dm2x
    save dat_dm_err/dm_err_10pm dm1x_3d dm2x_3d

    end
else
    %load dat_dm_err1/dm_err_1nm dm1x dm2x
    load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
end

dmxx = [dm1x dm2x];
rmso = stat2d(dmxx);

%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_150pm_itr20'];
    idd = 300
    disp('paused'), pause

% ******************************************
ff = 0 / rmso(1)
% ******************************************

    dm1 = dm1_save + dm1x * ff;
    dm2 = dm2_save + dm2x * ff;

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1550nm';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/eval_30nm_amp_opd_1064nm_at_500nm';
%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_bb_30nm_bw20_355nm_flat_m2'];
%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_sig_01nm_itr20'];
%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_opdrb_600pm_itr40'];
mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_opdrb_3000pm_itr20'];
%mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_1500nm_m2_1500pm_itr20'];

%load dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1064nm dm1x dm2x opdx ampx
%    mp.wfc.dm1x = dm1x;
%    mp.wfc.dm2x = dm2x;


% *******************************
    %mp.wfc.model_id = 7; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
% *******************************
    %Im = falco_get_summed_image(mp);

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm_500nm_distorted opd_nm_distorted amp_distorted
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm_1064nm opd_target amp_target 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm opd_target amp_target 

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw05_1064nm opd_target amp_target 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_355nm_flat opd_target amp_target ce_target

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target 
load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target

    mp.wfc.opd_target = opd_target;
    mp.wfc.amp_target = amp_target;
    mp.wfc.ce_target  = ce_target;

    %opd_nm_target = opd_target;
%save ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_4brandon opd_nm_target amp_target opd_nm_distorted amp_distorted
%return
    
%end % if - main
%% 

% *******************************
flag = 6
    mp.wfc.model_id = flag; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp, 8 = target_amp_opd_eval
% *******************************
    Im = falco_get_summed_image(mp);
    
    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

%end

    cx = [-7 -3];
    %cx = [-12 -8];
    fs = 14;


    figure(2), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title(sprintf('pre-efc-10nm: \\sigma = 1nm, Cb = %0.3g ',  cb0));
        title(sprintf('pre-efc-30nm: Cb = %0.3g ',  cb0));
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfi3

   end

    %load(mp.wfc.fn, 'dm1x', 'opdx', 'rmsx', 'dmrms')
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')
    rmsx = rms_opd;
    x = [1:size(rmsx,1)]-1;
    x = x';
    figure(1), clf, semilogy(x, rmsx(:,1)*1e3, 'ro-', x, dmrms(:,1)*1e3, 'bs-'), grid
    legend(sprintf('wfe-rms: last = %0.2fpm', rmsx(end,1)*1e3), 'dm1-rms [pm]', 'Location', 'Northeast')

return
%end % -----------------------------------------------
%%







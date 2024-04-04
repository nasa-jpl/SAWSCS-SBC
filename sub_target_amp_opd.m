flag = 1
disp('paused'), pause

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case2_20nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case3_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw05_10nm_2to12lamD_chg6_3lam_itr255 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw10_10nm_2to12lamD_chg6_4lam_itr350 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw15_10nm_2to12lamD_chg6_4lam_itr250 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw20_10nm_2to12lamD_chg6_4lam_itr150 dm1 dm2 tx cbx fom Im beta_value 

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw05_30nm_2to12lamD_chg6_4lam_itr380 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw10_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw15_30nm_2to12lamD_chg6_4lam_itr370 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 

load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr300 dm1 dm2 tx cbx fom Im beta_value 

pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
load dat_dm_err1/opd30_erkin dopd_c_nm

opdnm = dopd_c_nm ;

pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)


    mp.wfc.model_id = 7; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

Im = falco_get_summed_image(mp);

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    fs = 14;
    cx = [-12 -8];
    %cx = [-10 -6];

    %close all

    figure(1), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title([sprintf('post-wfc-30nm: Cb = %0.3g ',  cb0)]);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfi3

return
%end % -----------------------------------------------
%








flag = 10
disp('paused'), pause

if 0

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_10nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_20nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  

%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw05_30nm_2to12lamD_chg6_4lam_itr380 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw10_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw15_30nm_2to12lamD_chg6_4lam_itr370 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr340 dm1 dm2 tx cbx fom Im beta_value 
load /home/esidick/Afalco/falco20200916/macos/hex_dat1/efc1_gap8_256pix_bw20_00nm_2to12lamD_chg6_4lam_itr120 dm1 dm2 tx cbx fom Im bvalo 

    dm1_save = dm1;
    dm2_save = dm2;

load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm

%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
%pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/hex_maps/pupil_8mm_gap_8pix_256_matched.fits');

%load dat_dm_err1/opd10_erkin dopd_c_nm
%load dat_dm_err1/opd10_erkin dopd_c_nm
%load dat_dm_err1/opd30_erkin dopd_c_nm

%load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
load hex_maps/dat_opdnm opdnm %opdh mask    
dopd_c_nm = opdnm;

    opdnm = dopd_c_nm + 0 * opd_rb_drift_100pm * 1e-3;

    pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

%%

if 1
    %load IFopd/G_dm1_lim5_m1 G km1
    %load IFopd/G_dm2_lim5_m1 G km2
    %load IFopd/G_dm1_dm2_lim5_m2_1064nm G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm G km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm km1 km2 indx
    %load IFopd/G_dm1_dm2_lim5_m1_1064nm km1 km2 indx
    %load IFopd/dwddm_dm1_nm_over_nm_full indx maskopd
    load IFhex_opd/G_dm1_lim5_m0_500nm_00nm_bw20_gap08 G km1 km2 indx

    mp.wfc.indx = indx;
    mp.wfc.km1 = km1;
    mp.wfc.km2 = km2;
    mp.wfc.G = G;
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


if 0
    ma = 64;
    del = 1; % nm
    %dm1x = del * randn(ma,ma);
    %dm2x = del * randn(ma,ma);
    dm1x = 2*(0.5 - rand(ma,ma));
    dm2x = 2*(0.5 - rand(ma,ma));
    %dm1x(find(dm1x >  del)) = del;
    %dm1x(find(dm1x < -del)) =-del;
    %dm2x(find(dm2x >  del)) = del;
    %dm2x(find(dm2x < -del)) =-del;
    %save dat_dm_err1/dm_err_1nm dm1x dm2x
    save dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
else
    load dat_dm_err1/dm_err_1nm dm1x dm2x
    %load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
end

ff = 1;

    dm1 = dm1_save + dm1x * ff;
    dm2 = dm2_save + dm2x * ff;

    mp.dm1.V = dm1 * 0;
    mp.dm2.V = dm2 * 0;

%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1550nm';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/eval_30nm_amp_opd_1064nm_at_500nm';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_bb_30nm_bw05_1064nm_v2';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_30nm_500nm';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_20nm_500nm';
%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_10nm_500nm';
mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/IFhex_opd/wfc_dm1_flat_30nm_500nm_gray_gap08_m1';

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
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gray_gap08_wfe00_bw00_500nm opd_target amp_target ce_target 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gray_gap08_wfe00_bw00_500nm opd_target amp_target ce_target 
load ~esidick/Afalco/falco20200916/macos/IFhex_opd/target_opd_amp_gap08_00nm_bw20_500nm_post_efc opd_target amp_target ce_target
    mp.wfc.opd_target = opd_target * 0;
    mp.wfc.amp_target = amp_target * 0;
    mp.wfc.ce_target  = ce_target;

    %opd_nm_target = opd_target;
%save ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_4brandon opd_nm_target amp_target opd_nm_distorted amp_distorted
%return
    

% *******************************
    mp.wfc.model_id = 4; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp, 8 = target_amp_opd_eval
% *******************************
    Im = falco_get_summed_image(mp);

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    
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

return
end % -----------------------------------------------
%%
%load ~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_30nm_500nm_gray_gap08 dm1x opdx
%load ~esidick/Afalco/falco20200916/macos/IFhex_opd/wfc_dm1_flat_30nm_500nm_gray_gap08 dm1x  opdx rms_opd dmrms
load ~esidick/Afalco/falco20200916/macos/IFhex_opd/wfc_dm1_flat_30nm_500nm_gray_gap08_m1 dm1x  opdx rms_opd dmrms

%load hex_maps/dat_opdnm opdnm %opdh mask  
%dm1_flat = dm1x(:,:,end);
%save hex_maps/flatten_opd1_dm1 dm1_flat opdnm
%return

ns = size(opdx,3);
figure(1), clf
mx = [];
rmsx = [];

for ii = 1:20
    opd = opdx(:,:,ii);
    rmsi = stat2d(opd);
    cx = rmsi(1) * 3 * [-1 1];

    rmsx = [rmsx; rmsi];
    mx   = [mx; ii];

    figure(1), subplot(4,5,ii), imagesc(opd, cx), axis xy image, colormap(jet), niceaxes, 
        title(sprintf('%i: %0.1f',  ii, rmsi(1))); drawnow
end
   print -dpng ../Figs/fig_opd

    figure(2), clf
        plot(mx, rmsx(:,1), 'ro-'), grid
        parms = fun_newaxes(fs, 2, 5);  
        xlabel('WFC Iterations');
        ylabel('Residual RMS [nm]');
   print -dpng ../Figs/fig_rms

    figure(3), clf
        imagesc(opd, cx), axis xy image, colormap(jet), niceaxes, 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        title(sprintf('residual of flattened wfe'));
        xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rmsi));
   print -dpng ../Figs/fig_opd1

    aa = dm1x(:,:,end); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(4), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('flattening-dm1'));
        xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1
   
return
%end % -----------------------------------------------
%%






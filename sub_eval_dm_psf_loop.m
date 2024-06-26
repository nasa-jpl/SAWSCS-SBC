% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

    idd = 1
    %disp('paused'), pause

model_mode = 'eval';
pr2_hex_vvc

if 1


load([falcodir 'macos/hex_dat1/efc3_trial3_opd1_gap8_256pix_bw20_27nm_2to12lamD_chg6_4lam_itr410'],'dm1','dm2','tx','cbx','fom','Im','bvalo')
    dm1_save = dm1;
    dm2_save = dm2;

    tc = fom(end)*10;

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
load([falcodir 'macos/IFhex_opd/target_opd_amp_gap08_opd1_27nm_bw20_1500nm_post_efc'],'opd_target','amp_target','ce_target')
    mp.wfc.opd_target = opd_target;
    mp.wfc.amp_target = amp_target;
    mp.wfc.ce_target  = ce_target;

%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
%opdrb = opd_rb_drift_100pm;
%rmsrb = stat2d(opdrb);

if 0 % SAB2024
    %load dat_rb_wfe/opdrb_100pm opdrb_nm
    load dat_rb_wfe/opdrb_1000pm opdrb_nm
    ns = size(opdrb_nm, 3);
    ii = 2;
    opdrb = pad(opdrb_nm(:,:,ii), 256);
    rmsrb = stat2d(opdrb);

    % ********************************
    fac_rb = 10;  % 1, 3, 10, 30, 100 300, 1000, 3000 pm
    % ********************************
    %opdrb = opdrb * fac_rb / rmsrb(1); %[nm]
    %figure(1), clf, imagesc(opdrb); axis xy image, colormap(jet), colorbar, return
end
%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
pupil_save = fitsread([falcodir 'macos/hex_maps/pupil_8mm_gap_8pix_256_matched.fits']);
%load dat_dm_err1/opd10_erkin dopd_c_nm
%load dat_dm_err1/opd30_erkin dopd_c_nm

%load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
load hex_maps/trial3_opd1_final opdx %SAB2024
opdx_fromDave = opdx;
ns = size(opdx_fromDave, 3);
opdnm = 1e-3*opdx_fromDave(:,:,1);

%load ~esidick/2022habex/dat_macos/opdnom_256pix_macos opdnm

%dopd_c_nm = opdnm;

%    opdnm = dopd_c_nm + 0 * opd_rb_drift_100pm * 1e-3;
%    opdnm = opdnm0 + 0*opdrb * 1e0;

    pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

%%

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

flag = 1;

if 0 % ====================================================================

    mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, 4 = dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
    ii = 1
    %disp('paused'), pause
    mp.wfc.fn = sprintf('%smacos/IFhex_opd/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_1000pm_case%i', mp.falcodir, ii);
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')

id = 1;
    mp.dm1.V = dm1_save * id;
    mp.dm2.V = dm2_save * id;

Im = falco_get_summed_image(mp);
    psfn0 = Im .* maskc; 
    cb00 = mean(nonzeros(psfn0(:)));

    id1 = 2

    fs = 14;
    cx = [-12 -8];
    cx1 = [-7 -4];

    figure(1), clf
        imagesc(xv,xv,log10(psfn0)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title(['pre-wfc: 30nm-20%-12\lambda, dm-drift = 0pm ']);
        %xlabel([sprintf('Cb = %0.3g, Tc = %0.2f',  cb00, tc) '%']);
        xlabel([sprintf('Cb = %0.3g',  cb00)]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psf1
   return
    
%{
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp); %<***** thput ***

%end
    
figure(2), clf
    imagesc(xv,xv,log10(ImSimOffaxis)); axis xy image, grid
    %hold on, plot(xc, yc, 'y', xc1, yc1, 'r', 'Linewidth', 2), hold off
    parms = fun_newaxes(fs, 0, 0);  
    title([sprintf('input-30nm: with flattening, Tc = %0.2f', thput*100) '%']), 
print -dpng ../Figs/fig_psf3
%}

return
end % ====================================================================
%%



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
%else
    %load dat_dm_err1/dm_err_1nm dm1x dm2x
    %load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
    load dat_dm_err2/dm_err_3000pm dm1x_3d dm2x_3d
end

idd = 3000

%dmxx = [dm1x dm2x];
%rmso = stat2d(dmxx);
if 0 % SAB:  doesn't seem to do anything, probably another study's parameters
   % ******************************************
   ff = 0; %0 / rmso(1)
   % ******************************************
   load dat_rb_wfe/cbx_post_wfc_100pm cbxx psfn0 cb00
   psfn1 = psfn0;
   cb1 = cb00;
end

flag_movie = 0;

    if flag_movie
        video = VideoWriter('dat_dm_err1/movie_psf_1000pm.avi'); %create the video object
        video.FrameRate = 1;
        open(video); %open the file for writing
    end

cbxi = [];
cbxo = [];

for ii = 1:ns

    %opdrb = pad(opdrb_nm(:,:,ii), 256);
    %rmsrbi = stat2d(opdrb)*1e3; % pm

    opdnm = 1e-3*opdx_fromDave(:,:,ii); %+ opdrb;

    pupil = pupil_save .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    %mp.wfc.fn = sprintf('~esidick/Afalco/falco20200916/macos/dat_rb_wfe/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_opdrb_1000pm_run%i', ii);
    mp.wfc.fn = sprintf('%smacos/IFhex_opd/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_1000pm_case%0.2i', mp.falcodir,ii);
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdi', 'ampi', 'opdx', 'ampx', 'rms_opd', 'rms_amp','dmrms')

    %rmsrbi = stat2d(opdi)*1e3; % pm

    mp.dm1.V = dm1_save * 1;
    mp.dm2.V = dm2_save * 1;

    %dm1 = dm1x(:,:,1); % dm1_save + dm1x_3d(:,:,ii);
    %dm2 = dm2x(:,:,1); % dm2_save + dm2x_3d(:,:,ii);

    rmsrbo = stat2d(opdx)*1e3;

% *******************************
    flag = 1
    mp.wfc.model_id = flag; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp, 8 = target_amp_opd_eval
% *******************************
    Im = falco_get_summed_image(mp);
    psfi = Im .* maskc; 
    cbi = mean(nonzeros(psfi(:)))
    cbxi = [cbxi cbi];

    figure(1), clf % SAB We added these lines because having trouble using sub_plot_rb_psf
        imagesc(xv,xv,log10(psfi)); axis xy image, colormap(jet); colorbar
        title([sprintf('pre-wfc NI(%0.2i): Cb = %0.3g',  ii, cbi) ]), drawnow
	print('-dpng',sprintf('../Figs/preNI%0.2i',ii))

    % -------------------------------------------------
    mp.dm1.V = dm1x; %(:,:,end); % dm1_save + dm1x_3d(:,:,ii);
    mp.dm2.V = dm2x; %(:,:,end);

    Im = falco_get_summed_image(mp);
    psfo = Im .* maskc; 
    cbo = mean(nonzeros(psfo(:)))
    cbxo = [cbxo cbo];

    figure(2), clf % SAB We added these lines because having trouble using sub_plot_rb_psf
        imagesc(xv,xv,log10(psfo)); axis xy image, colormap(jet); colorbar
        title([sprintf('post-wfc NI(%0.2i): Cb = %0.3g',  ii, cbo) ]), drawnow
	print('-dpng',sprintf('../Figs/postNI%0.2i',ii))

    cx = [-11.5 -8.5];
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

%save dat_rb_wfe/cbx_post_wfc_100pm cbxx psfn0 cb00
%save dat_rb_wfe/cbx_pre_pos_wfc_1000pm cbxi cbxo  
%save dat_dm_err/cbx_pre_pos_wfc_3000pm cbxi cbxo   % SAB cbxi mean contrast with drift, cbxo mean contrast w/o drift

cmean = mean(cbxo)

x = [1:10];
figure(2), clf, semilogy(x, cbxi, 'ro-', x, cbxo, 'bs-'), grid
legend('pre-wfc NI','post-wfc NI'), axis([1 10 7.7e-10 8.2e-10])
%axis([1 10 1e-10 1e-9])

return
%end % -----------------------------------------------
%%

    %load(mp.wfc.fn, 'dm1x', 'opdx', 'rmsx', 'dmrms')
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms')
    rmsx = rms_opd;
    x = [1:size(rmsx,1)]-1;
    x = x';
    figure(1), clf, semilogy(x, rmsx(:,1)*1e3, 'ro-', x, dmrms(:,1)*1e3, 'bs-'), grid
    legend(sprintf('wfe-rms: last = %0.1fpm', rmsx(end,1)*1e3), 'dm1-rms [pm]', 'Location', 'Northeast')
    drawnow
return
end % -----------------------------------------------
%%
load dat_dm_err/cbx_pre_pos_wfc_10pm cbxi cbxo  
    ya1 = cbxi; yb1 = cbxo;
load dat_dm_err/cbx_pre_pos_wfc_30pm cbxi cbxo  
    ya2 = cbxi; yb2 = cbxo;
load dat_dm_err/cbx_pre_pos_wfc_100pm cbxi cbxo  
    ya3 = cbxi; yb3 = cbxo;
load dat_dm_err/cbx_pre_pos_wfc_300pm cbxi cbxo  
    ya4 = cbxi; yb4 = cbxo;
load dat_dm_err/cbx_pre_pos_wfc_1000pm cbxi cbxo  
    ya5 = cbxi; yb5 = cbxo;
load dat_dm_err/cbx_pre_pos_wfc_3000pm cbxi cbxo  
    ya6 = cbxi; yb6 = cbxo;

    y1 = [mean(ya1) mean(ya2) mean(ya3) mean(ya4) mean(ya5) mean(ya6)];
    y2 = [mean(yb1) mean(yb2) mean(yb3) mean(yb4) mean(yb5) mean(yb6)];

legx{1} = 'dm-drift = 10pm';
legx{2} = 'dm-drift = 30pm';
legx{3} = 'dm-drift = 100pm';
legx{4} = 'dm-drift = 300pm';
legx{5} = 'dm-drift = 1000pm';
legx{6} = 'dm-drift = 3000pm';

x = [1:10];

figure(1), clf, semilogy(x,ya1, 'ro-', x,ya2, 'bs-', x,ya3, 'gd-',x,ya4, 'm^-', x,ya5, 'c>-', x,ya6, 'ko-', x,yb1, 'ro--', x,yb2, 'bs--', x,yb3, 'gd--', x,yb4, 'm^--', x,yb5, 'c--', x,yb6, 'ko--', x, x*0+7.8e-11, 'g--'), grid
    parms = fun_newaxes(12, 2, 5);  
    xlabel('random error case number')
    ylabel('broadband mean contrast, Cb')
    title('solid: pre-wfc. dashed: post-wfc')
    axis([1 10 5e-11 1e-5])
    legend(legx,'Location', 'North')
print -dpng ../Figs/fig_cb    

x1 = [10 30 100 300 1000 3000];

figure(2), clf, loglog(x1,y1, 'ro-', x1,y2, 'bs-', x1, x1*0+7.8e-11, 'm--'), grid
    parms = fun_newaxes(12, 2, 5);  
    xlabel('drift rms [pm]')
    ylabel('broadband mean contrast, Cb')
    title('average values of 10 cases')
    axis([10 3000 5e-11 1e-5])
    legend('pre-wfc', 'post-wfc', 'nominal','Location', 'Northwest')
print -dpng ../Figs/fig_cb1   

return
%end % -----------------------------------------------
%%




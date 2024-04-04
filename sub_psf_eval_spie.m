id_update = 20
disp('paused'), pause

fs = 14;


flag_rb = 0;

% ********************************
    rms_rb = 3000;  % 1, 3, 10, 30, 100 pm
% ********************************

%if 1

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
    dm1_save = dm1;
    dm2_save = dm2;

load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
    dmxx = [dm1x dm2x];
    rmso = stat2d(dmxx);

% ******************************************
    ff = 0 / rmso(1);
    ff = 0.3;
% ******************************************

    dm1a = dm1 + dm1x * ff;
    dm2a = dm2 + dm2x * ff;

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;
    ms = length(find(abs(xv)<13));

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);


load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
    opdrb = opd_rb_drift_100pm;
    rmsrb = stat2d(opdrb);
    opdrb = opdrb * rms_rb / rmsrb(1); %[pm]

%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
%load dat_dm_err1/opd30_erkin dopd_c_nm

pupil_map = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
opdnm0 = opdnm;

    %opdnm = dopd_c_nm + fac * opd_rb_drift_100pm * 1e-3;
    if flag_rb, opdnm = opdnm0 + opdrb * 1e-3; end

    pupil = pupil_map .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

%%

% *******************************
    mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
% *******************************
    %mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20_new'];
    mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20'];
    load(mp.wfc.fn, 'psfx', 'cbx');

    mp.wfc.fn = ['~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_03nm_itr20'];
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp', 'dmrms');
    
    rmsx = rms_opd;
    
    ns = size(opdx, 3);
    %ns = 5;

    

if 1    

    mp.dm1.V = dm1_save;
    mp.dm2.V = dm2_save;
    %mp.dm1.V = dm1x(:,:,1);
    %mp.dm2.V = dm2x; %(:,:,1);
    mp.dm1.V = dm1a;
    mp.dm2.V = dm2a;

    pupil = pupil_map .* exp(1i*2*pi*opdnm0*1e-9/mp.lambda0);
    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    Im = falco_get_summed_image(mp);

    pupil = pupil_map .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    psfn0 = Im .* maskc; 
    cb0 = mean(nonzeros(psfn0(:)));

    fs = 14;
    cx = [-11 -8];
    cx = [-10 -7];

    figure(1), clf
        imagesc(xv,xv,log10(psfn0), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-dm-50pm-12\\lambda: Cb = %0.3g ',  cb0)]);
        title([sprintf('pre-wfc: drift = 180pm, Cb = %0.3g ',  cb0)]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfi
   return
%end


%load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
pmax = max(amp_target(:))

flag_movie = 0;

if flag_movie
    video = VideoWriter('dat_img/erkin_movie_dm_drift_1500pm_1500nm_dm1_dm2_wfc.avi'); %create the video object
    video.FrameRate = 2;
    open(video); %open the file for writing
end

   %mx  = [0]; cbx = [cb0];
   mx = []; cbx = []; 
   rms_opd = [];
   rms_amp = [];

   clear psfx
   figure(1), clf

   for ii = 21 %1:21

        mp.dm1.V = dm1x(:,:,ii);
        %mp.dm2.V = dm2x(:,:,ii); % for dm1/dm2 wfc
        %mp.dm2.V = dm2a;  % for dm1-wfc only
        mp.dm2.V = dm2a;  % for dm1-wfc only

        Im = falco_get_summed_image(mp);
        psfn = Im .* maskc; 
        cbi = mean(nonzeros(psfn(:)))

        psfx(:,:,ii) = psfn;

        ampi = ampx(:,:,ii) / pmax;
        opdi = opdx(:,:,ii) * 1e3;

        mx = [mx ii];
        cbx = [cbx cbi];

        aa = opdi; rmsa = stat2d(aa);     rms_opd = [rms_opd; rmsa];     cx_opd = rmsa(1)*3*[-1 1]; 
        bb = ampi; rmsb = stat2d(bb);     rms_amp = [rms_amp; rmsb];     cx_amp = rmsb(1)*2*[-1 1]; 

        if 1

    fp = figure(1); colormap jet; fp.Position = [24   547   812   251];

    subplot(1,3,1), imagesc(bb,cx_amp), axis xy image, colormap(jet), colorbar, niceaxes
        title(sprintf('itr%i: \\Deltaamp', ii)), 
        xlabel(sprintf('rms = %0.3g', rmsb(1))), 

    subplot(1,3,2), imagesc(aa,cx_opd), axis xy image, colormap(jet), colorbar, niceaxes
        title(sprintf('itr%i: \\Deltaopd [pm]', ii)), 
        xlabel(sprintf('rms = %0.1fpm', rmsa(1))), 

    subplot(1,3,3), imagesc(log10(pad(psfn,ms)), [-11 -8]), axis xy image, colormap(jet), colorbar, niceaxes
        title(sprintf('itr%i: contrast, Cb = %0.3g', ii, cbi)), 

        end

if flag_movie

        sub_plot_spie_drift

        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(video,frame);
end

end % for ii - loop
      
    if flag_movie, close(video); end %close the file

%end % main - if

figure(3), clf, semilogy(mx, cbx, 'ro-', mx, mx*0+cb0,'b--'), grid, title('Cb'), drawnow
figure(4), clf, semilogy(mx', rms_opd(:,1)*1e3, 'ro-'), grid, title('opd'), drawnow
figure(5), clf, semilogy(mx', rms_amp(:,1), 'ro-'), grid, title('amp'), drawnow

    %save dat_dm_err1/post_wfc_case6_30nm mx cbx rms_opd rms_amp psfx
    %save dat_dm_err1/post_wfc_30nm_bw20_1064nm_flat mx cbx rms_opd rms_amp psfx
    %save dat_dm_err1/post_wfc_gap08_30nm_bw20_500nm_12lam_sig09 mx cbx rms_opd rms_amp psfx psfn0
    %save dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20_new mx cbx rms_opd rms_amp psfx psfn0
    %save dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20 mx cbx rms_opd rms_amp psfx psfn0
    %save dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_200pm_itr40 mx cbx rms_opd rms_amp psfx psfn0
    %save dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_3000pm_itr20_1500nm mx cbx rms_opd rms_amp psfx psfn0 cb0
    %save dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_1500pm_itr20_1500nm mx cbx rms_opd rms_amp psfx psfn0

return
%end % main - if

        mp.dm1.V = dm1a;
        mp.dm2.V = dm2a;

        Im1 = falco_get_summed_image(mp);
        psfn = Im1 .* maskc; 
        cbi = mean(nonzeros(psfn(:)));

    figure(3), clf
        imagesc(xv,xv,log10(psfn), [-10 -7]); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('pre-wfc-\\gamma09-12\\lambda: Cb = %0.3g ',  cbi)]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfi2

return
%end % -----------------------------------------------
%%
%load dat_dm_err1/post_wfc_30nm_bw20_1064nm_flat mx cbx rms_opd rms_amp psfx

    figure(4), clf, semilogy(cbx, 'ro-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('WFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('input rms = 30nm, bw = 20%, 12\labda, 500nm, \beta=-2')
        %legend(s1, 'Location', 'Northeast')
    print -dpng ../Figs/fig_wfc

return
%end % -----------------------------------------------
%%

psfn = psfx(:,:,1); %Im .* maskc; 
    cb = mean(nonzeros(psfn(:)));

    cx = [-11 -8];

    figure(1), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('pre-dm-50pm-12\\lambda: Cb = %0.3g ',  cb)]);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfi_dm_50

psfn = psfx(:,:,end); %Im .* maskc; 
    cb = mean(nonzeros(psfn(:)));

    cx = [-11 -8];

    figure(2), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('post-wfc: drift = 180pm, Cb = %0.3g ',  cbx(end))]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfo

return
end % -----------------------------------------------
%%
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm           opd_target amp_target ce_target
%load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
a = max(amp_target(:))

aa = ampx(:,:,end)/a; rms0 = stat2d(aa); cx = rms0(1)*2*[-1 1];

figure(1), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    %hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('post-wfc \\Deltaamp: rms = %0.3g', rms0(1))), 
print -dpng ../Figs/fig_ampo

aa = opdx(:,:,end)*1e3; rms0 = stat2d(aa); cx = rms0(1)*2*[-1 1];

figure(2), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('post-wfc \\Deltaphase: rms = %0.1fpm', rms0(1))), 
print -dpng ../Figs/fig_opdo

return
%end % -----------------------------------------------
%%



% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;
addpath(genpath('/proj/exospin/users/joegreen/'));

%fn = '/proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin'
%return

if 0

    %load dat/dat_ep2_ce_no_dm_20230209 ce
    %load dat/dat_ep2_ce_with_dm_20230209 ce
    load dat/dat_ep3_ce_at_apod_no_dm_20230209 ce
    ce1 = ce;
    load dat/dat_ep3_ce_at_apod_with_dm2_poke_20230209 ce

np = 512;
amp  = abs(ce1);
amx = max(amp(:));
mask = amp/amx > 0.1;
amp = (abs(ce) - amp) / amx;

lam0 = 500;

opd = mask .* real((angle(ce) - angle(ce1))* lam0 / 2 / pi);
[opdf, zk] = zernike_remove(opd, [1]);
rmsa = stat2d(opdf); cx = rmsa(1)*3*[-1 1];

nx = [-np/2:np/2-1];
[xc,yc] = fun_circle(nx, np/2);

figure(1), clf, imagesc(nx,nx,amp);  axis xy image; colormap(jet), niceaxes
    hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('\Deltaamp at apod: falco, dm2-act(32,32) = 5nm'), 
    xlabel(sprintf('rms = %0.3g, pv = %0.3g', stat2d(amp))), 
print -dpng ../Figs/fig_amp1

figure(2), clf, imagesc(nx,nx,opdf);  axis xy image; colormap(jet), niceaxes
    hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    %title('opd at apod: falco, no dm'), 
    title('\Deltaopd at apod: falco, dm2-act(32,32) = 5nm'), 
    xlabel(sprintf('rms = %0.2f, pv = %0.1f nm', rmsa)), 
print -dpng ../Figs/fig_opd1

return
%end % -----------------------------------------------
%%



load maps/opd_dm_post_maint_4falco 

load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 

    aa = [dm1 dm1_new_nm_4falco]; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

    figure(1), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title('dm1: pre-drift                      post-drift')
        xlabel(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = [dm2 dm2_new_nm_4falco]; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

    figure(2), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title('dm2: pre-drift                      post-drift')
        xlabel(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2

    aa = opd_dm1_new_nm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

    figure(3), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('new opd: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_opd
   

return
%end % -----------------------------------------------
%%

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_00mm_bw20_30nm_2to12lamD_chg6_4lam_itr400 dm1 dm2 tx cbx fom Im bvalo 
    x1 = [0:length(cbx)-1]; y1 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
    x2 = [0:length(cbx)-1]; y2 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_16mm_bw20_30nm_2to12lamD_chg6_4lam_itr420 dm1 dm2 tx cbx fom Im bvalo 
    x3 = [0:length(cbx)-1]; y3 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_24mm_bw20_30nm_2to12lamD_chg6_4lam_itr440 dm1 dm2 tx cbx fom Im bvalo 
    x4 = [0:length(cbx)-1]; y4 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/sharp_gap_24mm_bw20_30nm_2to12lamD_chg6_4lam_itr360 dm1 dm2 tx cbx fom Im bvalo 
    x5 = [0:length(cbx)-1]; y5 = cbx;

s1 = ['gap = 00mm: ' sprintf('final Cb = %0.3g', y1(end))];
s2 = ['gap = 08mm: ' sprintf('final Cb = %0.3g', y2(end))];
s3 = ['gap = 16mm: ' sprintf('final Cb = %0.3g', y3(end))];
s4 = ['gap = 24mm: ' sprintf('final Cb = %0.3g', y4(end))];
s5 = ['gap = 24mm-sharp: ' sprintf('final Cb = %0.3g', y5(end))];

end

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
    x1 = [0:length(cbx)-1]; y1 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_flat_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
    x2 = [0:length(cbx)-1]; y2 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_opdhi_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
    x3 = [0:length(cbx)-1]; y3 = cbx;

    
s1 = ['nominal:   ' sprintf('final Cb = %0.3g', y1(end))];
s2 = ['flattened: ' sprintf('final Cb = %0.3g', y2(end))];
s3 = ['filtered:    ' sprintf('final Cb = %0.3g', y3(end))];


fs = 14;

    %figure(1), clf, semilogy(x1, y1, 'r-',x2, y2, 'b-',x3, y3, 'g-',x4, y4, 'm-', x5,y5,'k-'), grid, 
    figure(1), clf, semilogy(x1, y1, 'r-',x2, y2, 'b-',x3, y3, 'g-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Cb vs Itr: gap = 8mm, 4-\lambda')
        legend(s1, s2,s3,'Location', 'Northeast')
        axis([0 420 1e-10 1e-4])
    print -dpng ../Figs/fig_efc

return
%end % -----------------------------------------------
%%
load dat_dm_err/flatten_init_wfc_dmx_opdx dmx opdx

    aa = real(dmx(:,:,end)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = real(opdx(:,:,end)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(15), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-flattened opd'));
        xlabel(sprintf('rms = %0.3f, pv = %0.2fnm', rms0));
   print -dpng ../Figs/fig_opd
   
   %save maps/du_dm1_dm2_nm_20221219 dm1 dm2
   %save maps/du_dm1_dm2_nm_20221228 dm1 dm2 opd_dm1_new_nm

return
%end % -----------------------------------------------
%%
flag = 1

%load IFopd/wfc_dmx_opdx dmx opdx
%load dat_dm_err/flatten_init_wfc_dmx_opdx dmx opdx
%load dat_dm_err1/wfc_dm1_dm2_opdx_10nm dm1x dm2x opdx ampx
%load dat_dm_err/flatten_init_wfc_dmx_opdx dmx opdx

%load dat_dm_err1/wfc_dm1_dm2_case1_10nm dm1x dm2x opdx ampx
%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm dm1x dm2x opdx ampx
%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_355nm dm1x dm2x opdx ampx
%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1064nm dm1x dm2x opdx ampx

%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_bb_30nm_bw05_1064nm_v2 dm1x dm2x opdx ampx
load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_bb_30nm_bw20_355nm_flat_m2 dm1x dm2x opdx ampx
%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_30nm_500nm dm1x opdx

ns = size(opdx,3);
%aa = opdnm; rmsi = stat2d(aa); cx = rmsi(1)*3*[-1 1];
aa = opdx(:,:,1); %amp_target; 
rmsi = stat2d(aa); cx = rmsi(1)*3*[-1 1];

nrow = 4;
ncol = 4;

figure(2), clf, jj = 0;
    figure(2), subplot(nrow,ncol,jj+1), imagesc(aa), axis xy image, colormap(jet), niceaxes, 
    %figure(2), clf, imagesc(aa), axis xy image, colormap(jet), niceaxes,  colorbar
        title(sprintf('%i', jj)), drawnow
        xlabel(sprintf('rms=%0.3f', rmsi(1))), drawnow
%return

rmsx = [];
rmsa = [];

for jj = 1:ns

    opdi = real(opdx(:,:,jj));
    %opdf = real(ampx(:,:,jj));
    [opdf, zk] = zernike_remove(opdi, [1]);

    if 0
        m = 256; 
        x = [1:m*m];
        bb = opd(:);
        k1 = find(bb>0); ap = bb(k1); meanp = mean(ap)*100;
        k2 = find(bb<0); am = bb(k2); meanm = mean(am)*100;
        k1 = find(bb>meanp); bb(k1) = 0; 
        k2 = find(bb<meanm); bb(k2) = 0; 
        %opdf = reshape(bb,m,m);
        %mask = abs(opd) ~= 0;
        %[dw,indx1] = m2v(mask);
    end

    %dw = m2v(opd,indx1);
    %opdf = v2m(dw,indx1);

    %da = m2v(ampx(:,:,jj),indx1);
    %ampi = v2m(da,indx1);
    
    rmsi = stat2d(opdf); 
    cx = rmsi(1)*3*[-1 1];
    rmsx = [rmsx; rmsi];

    ampi = ampx(:,:,jj);
    rmsii = stat2d(ampi); 
    rmsa = [rmsa; rmsii];

   %figure(4), clf, plot(x, sort(opd(:)),'ro-'), grid
   %hold on, plot(x,x*0+meanp,'b--'), hold off
   %hold on, plot(x,x*0+meanm,'g--'), hold off
   %title(sprintf('jj = %i',jj))
   %disp('paused'), pause

    %figure(2), subplot(nrow,ncol,jj), imagesc(opdf, cx), axis xy image, colormap(jet), niceaxes, 
     %   title(sprintf('itr=%i', jj-1)), 
     %   xlabel(sprintf('%0.2f', rmsi(1))), drawnow
end
   %print -dpng ../Figs/fig_opdo

   mx = [1:size(rmsx,1)]-1;

figure(10), clf, semilogy(mx, rmsx(:,1)*1e3, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('pm-rms of \Deltawfe at apodizer')
    title('dm1/dm2-wfc with 355nm, \beta = -2')
    %title('wf-flattening with dm1: \lambda = 500nm, \beta = -2')
    %axis([0 10 0.02 0.03])
    legend(sprintf('last rms = %0.1fpm', rmsx(end,1)*1e3), 'Location', 'Northwest')
print -dpng ../Figs/fig_rms1

%end

figure(11), clf, semilogy(mx, rmsa(:,1), 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp at apodizer')
    title('dm1/dm2-wfc with 355nm, \beta = -2')
    legend(sprintf('last rms = %0.3g', rmsa(end,1)*1e0), 'Location', 'Northwest')
print -dpng ../Figs/fig_rms2

%return

opd = opdf;

m = 256; 
bb = opd(:);
kk = find(bb>0); ap = bb(kk); mean1 = mean(ap);
kk = find(bb<0); am = bb(kk); mean2 = mean(am);

kk = find(bb > 100); bb(kk) = 0; 
kk = find(bb <-100); bb(kk) = 0; opd = reshape(bb,m,m);

opd = opdf * 1e3;

    rmsi = stat2d(opd); 
    cx = rmsi(1)*3*[-1 1];

    figure(1), clf, imagesc(opd, cx), axis xy image, colormap(jet), niceaxes, 
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-wfc: \\beta = -2: rms = %0.2fpm', rmsi(1))), drawnow
   print -dpng ../Figs/fig_opd1

   figure(4), clf, plot(sort(opdf(:)),'ro-'), grid


%save dat_dm_err1/dm1_dm2_opd_rms_30nm_bw10_1064nm rmsx rmsa
return
%end % -----------------------------------------------
%%
load dat_dm_err1/dm1_opd_rms_10nm rmsx
    y1 = rmsx(:,1); x1 = [1:length(y1)]'-1;
load dat_dm_err1/dm2_opd_rms_10nm rmsx
    y2 = rmsx(:,1); x2 = [1:length(y2)]'-1;
load dat_dm_err1/dm1_dm2_opd_rms_10nm rmsx
    y3 = rmsx(:,1); x3 = [1:length(y3)]'-1;

figure(10), clf, plot(x1,y1, 'ro-', x2,y2, 'bs-', x3,y3, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('nm-rms of \Deltawfe at apodizer')
    %ylabel('rms of \Deltaamp at apodizer')
    title('\\sigma = 1.0 nm')
    %axis([0 10 0.02 0.03])
    legend('dm1-wfc', 'dm2-wfc', 'dm1/dm2-wfc', 'Location', 'Northeast')
print -dpng ../Figs/fig_rms

return
%end % -----------------------------------------------
%%

%load /home/esidick/Afalco/falco20200916/dat_macos_20230125/sixMST_cir97_2to12lamD_chg6_1lam_run1_cbx_itr50 cbx 
load /home/esidick/Afalco/falco20200916/dat_macos_20230125/d200_cir97_2to12lamD_chg6_1lam_run1_cbx_itr65 cbx 
figure(1), clf, semilogy(cbx, 'ro-'), grid
cb0 = cbx(end);
x1 = [0];
cbxx = [cb0];

Colx = ['ro- '; 'bs- '; 'gd- '; 'm^- '; 'c>- ';'ks- '];
%load dat_dm_err/opd_rms_01nm rmsx
load dat_dm_err/cbx_10nm cbx

xx = [1 2:2:10];
%fn = 'dat_dm_err/opd_rms_%02inm.mat';
fn = 'dat_dm_err/cbx_%02inm.mat';
ns = length(xx);
%mx = [1:size(cbx,1)]-1;
mx = [1:size(cbx,2)]-1;
figure(10), clf
clear legx

for jj = 1:ns
    fni = sprintf(fn,xx(jj))
    %load(sprintf(fn,xx(jj)), 'rmsx');
    load(sprintf(fn,xx(jj)), 'cbx');
    figure(10), semilogy(mx, cbx, Colx(jj,:)), hold on
    legx{jj} = sprintf('\\sigma = %0.1fnm: Cb = %0.3g', xx(jj)/10, cbx(end));

    x1 = [x1 xx(jj)/10];
    cbxx = [cbxx cbx(end)];
end

figure(10), semilogy(mx, mx*0+cb0, 'g--'), hold off, grid
    legx{jj+1} = sprintf('\\sigma = 0.0nm: Cb = %0.3g', cb0);
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    %ylabel('nm-rms of \Deltawfe at apodizaer')
    ylabel('Mean Contrast, Cb')
    legend(legx, 'Location', 'Northeast')
    axis([0 20 5e-11 1e-7])
print -dpng ../Figs/fig_cb
    
figure(9), semilogy(x1, cbxx, 'ro-'), grid
    legx{jj+1} = 'nominal';
    parms = fun_newaxes(14, 2, 5);  
    xlabel('dm error std, \sigma [nm]')
    %ylabel('nm-rms of \Deltawfe at apodizaer')
    ylabel('Mean Contrast, Cb')
    %legend(legx, 'Location', 'Northeast')
    axis([0 1 5e-11 6e-10])
print -dpng ../Figs/fig_cb1

return
%end % -----------------------------------------------
%%
%load IFopd/dm1_opdtarget dm1 opd_target
%load IFopd/dat_post_wfc_opd opd
%load IFopd/target_opd_amp_bb_30nm_bw20_500nm_post_efc opd_target amp_target ce_target
load IFopd/target_opd_amp_bb_30nm_bw20_1064nm_flat opd_target amp_target ce_target

%[opdf, zk] = zernike_remove(opd, [1]);


    aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('target: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('Phase: RMS = %0.1fnm', rms0(1)));
   print -dpng ../Figs/fig_opd1

    aa = amp_target; aa = aa/max(aa(:)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(15), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('itr1 \\Deltaopd: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('Amplitude: RMS = %0.2f', rms0(1)));
   print -dpng ../Figs/fig_amp1

return
%end % -----------------------------------------------
%%
%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10 dopd_c_nm opd_10nm opd_c_nm opd_nom_nm
%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd20 dopd_c_nm opd_20nm opd_c_nm opd_nom_nm
%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd30 dopd_c_nm opd_30nm opd_c_nm opd_nom_nm

load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opdc30 dopd_c_nm opd_30nm opd_c_nm opd_nom_nm
%aa =  '/proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/'

rms0 = stat2d(dopd_c_nm);
opd1 = dopd_c_nm * 2 / rms0(1);

load ~esidick/2023_6mst/6mst/dat/pds_opd_nm psdnm

rms0 = stat2d(psdnm);
opd2 = psdnm / rms0(1);

%aa = dopd_c_nm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
%aa = opd1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
aa = opd2+opd1; rms1 = stat2d(aa); cx = rms0(1)*3*[-1 1];
aa = (opd2+opd1)*30/rms1(1); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

dopd_c_nm = aa;
save dat_dm_err1/opd30_erkin dopd_c_nm

    figure(2), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('scaled sum opd: rms = %0.1f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_opd1

return
%end % -----------------------------------------------
%%
%id = 2, disp('paused'), pause

%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_3pm opd_rb_drift_10pm opd_rb_drift_30pm opd_rb_drift_100pm
%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10_0323
load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd20_0323
%fn = '/proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/'
%whos, return

%aa = opd_rb_drift_100pm; rms0 = stat2d(aa); cx = rms0(1) * 3 * [-1 1];
aa = dopd_c_20nm; rms0 = stat2d(aa);
aa = aa * 2 / rms0(1);
rms0 = stat2d(aa); cx = rms0(1) * 3 * [-1 1];

    figure(1), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('residual opd-20nm: rms = %0.1f', rms0(1)));
   print -dpng ../Figs/fig_opd1

%return
%end % -----------------------------------------------
%%

rms_input = 20;

rms0 = stat2d(aa);
%opd1 = dopd_c_nm * 2 / rms0(1);
opd1 = aa * 2 / rms0(1);

%load ~esidick/2023_6mst/6mst/dat/pds_opd_nm psdnm
load ~esidick/2023_6mst/6mst/dat/pds_opd_nm_case3 psdnm

rms0 = stat2d(psdnm);
opd2 = psdnm / rms0(1);

    figure(2), clf, imagesc(opd2), axis xy image, colormap(jet), niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('scaled psd-opd: rms = 1'))
    print -dpng ../Figs/fig_psdm

%aa = dopd_c_nm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
%aa = opd1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
aa = opd2+opd1; rms1 = stat2d(aa); cx = rms0(1)*3*[-1 1];

aa = (opd2+opd1) * rms_input /rms1(1); 
rms0 = stat2d(aa); 
cx = rms0(1)*3*[-1 1];

dopd_c_nm = aa;
%save dat_dm_err1/opd30_erkin dopd_c_nm
%save dat_dm_err1/opd20_erkin dopd_c_nm
%save dat_dm_err1/opd10_erkin dopd_c_nm

    figure(3), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('scaled/combined opd: rms = %0.1f nm', rms0(1)));
   print -dpng ../Figs/fig_opd2

return
%end % -----------------------------------------------
%%
%load ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_ce_at_apod_with_dm_rb_errors amp opd ce
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_1064nm_flat opd_target amp_target

%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_00nm_bw20_500nm opd_target amp_target ce_target
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_1064nm opd_target amp_target 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_1064nm_flat opd_target amp_target 

load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_355nm_flat opd_target amp_target ce_target
load ~esidick/Afalco/falco20200916/macos/IFopd/distorted_opd_amp_bb_30nm_bw20_355nm_flat opd_target amp_target ce_distorted
    
    %aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    aa = (amp_target>0).*angle(ce_distorted./(ce_target+eps))*355/2/pi; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('\\Deltaopd: post-efc, input = 30nm, flattened'));
        xlabel(sprintf('rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_opd

    aa = amp_target; aa = aa / max(aa(:)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('amp: post-efc, input = 30nm, flattened'));
        xlabel(sprintf('rms = %0.2f, pv = %0.2f', rms0));
   print -dpng ../Figs/fig_amp

return
%end % -----------------------------------------------
%%
id = 3

clear legx
%load dat_dm_err1/post_wfc_case4_10nm mx cbx rms_opd rms_amp %psfx
load dat_dm_err1/post_wfc_case6_30nm_355nm mx cbx rms_opd rms_amp %psfx
    ns = length(cbx); kk = 1:12;
    x1 = mx(kk); y1 = cbx; 
    %x1 = mx(kk); y1 = rms_opd(:,1); legx{1} = 'input-rms = 10nm';
%load dat_dm_err1/post_wfc_case5_20nm mx cbx rms_opd rms_amp 
load dat_dm_err1/post_wfc_case6_30nm_500nm mx cbx rms_opd rms_amp %psfx
    x2 = mx(kk); y2 = cbx; y2(1) = 2.75e-9; 
    %x2 = mx(kk); y2 = rms_opd(:,1)*1e3; legx{2} = 'input-rms = 20nm';
    %x2 = mx(kk); y2 = rms_amp(:,1); legx{2} = 'input-rms = 20nm';
%load dat_dm_err1/post_wfc_case6_30nm mx cbx rms_opd rms_amp %psfx
load dat_dm_err1/post_wfc_case6_30nm_1064nm mx cbx rms_opd rms_amp %psfx
    x3 = mx(kk); y3 = cbx; 
    %x3 = mx(kk); y3 = rms_opd(:,1); legx{3} = 'input-rms = 30nm';
load dat_dm_err1/post_wfc_case6_30nm_1550nm mx cbx rms_opd rms_amp %psfx
    x4 = mx(kk); y4 = cbx; 
    %x4 = mx(kk); y4 = rms_opd(:,1); legx{4} = 'input-rms = 30nm';

    x = [1:11]';

legx{1} = 'wfc-\lambda = 355nm';
legx{2} = 'wfc-\lambda = 500nm';
legx{3} = 'wfc-\lambda = 1064nm';
legx{4} = 'wfc-\lambda = 1550nm';
legx{5} = 'nominal: wfc-\lambda = 1550nm';

figure(1), clf, semilogy(x1,y1, 'ro-', x2,y2, 'bs-', x3,y3, 'gd-', x4,y4, 'm^-', x1, y1*0+1.74e-11, 'c--'), grid, 
%figure(3), clf, semilogy(x,y1, 'ro-', x,y2, 'bs-', x,y3, 'gd-', x,y4, 'm^-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('mean contrast, Cb')
    %ylabel('residual opd-rms [pm]')
    %ylabel('\Deltaamp-rms')
    title('Cb vs WFC-Iterations')
    %title('\Deltaopd vs WFC-Iterations')
    %axis([0 12 1e-4 1e-2])
    legend(legx, 'Location', 'Northeast')
print -dpng ../Figs/fig_cb

return
%end % -----------------------------------------------
%%
load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_10nm_v3_2to12lamD_chg6_4lam_itr320 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_1dm_noflat_bw20_30nm_2to12lamD_chg6_4lam_itr10 dm1 dm2 tx cbx fom Im beta_value 
figure(2), clf, imagesc(dm2), axis xy image, colorbar
    x1 = [1:length(cbx)]-1; y1 = cbx; legx{1} = 'input-rms = 10nm';
load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_20nm_v3_2to12lamD_chg6_4lam_itr320 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_1dm_dm2_bw20_30nm_2to12lamD_chg6_4lam_itr160 dm1 dm2 tx cbx fom Im beta_value 
    x2 = [1:length(cbx)]-1; y2 = cbx; legx{2} = 'input-rms = 20nm';
load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr320 dm1 dm2 tx cbx fom Im beta_value 
    x3 = [1:length(cbx)]-1; y3 = cbx; legx{3} = 'input-rms = 30nm';

%figure(1), clf, semilogy(x1,y1, 'r-', x2,y2, 'b-', x3,y3, 'g-'), grid, 
figure(1), clf, semilogy(x1,y1, 'ro-', x2,y2,'bs-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('efc iterations')
    ylabel('mean contrast, Cb')
    %title('Cb vs EFC-Iterations: \\sigma = 1.0 nm')
    title('Cb vs EFC-Iterations')
    legend(legx, 'Location', 'Northeast')
    %axis([0 100 1e-11 1e-3])
print -dpng ../Figs/fig_efc

return
%end % -----------------------------------------------
%%
load dat_dm_err1/dm_err_1nm dm1x dm2x
y1 = [dm1x(:); dm2x(:)];

x = [-1:0.1:1];
y = hist(y1, x);
y = 100 * y / length(y1);

figure(1), clf, plot(x,y, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('dm drift error [nm]')
    ylabel('histogram [%]')
    title('truncated normal dist with \sigma  =1nm')
print -dpng ../Figs/fig_hist


return
%end % -----------------------------------------------
%%

load dat_dm_err1/post_wfc_case6_30nm mx cbx rms_opd rms_amp psfx

    video = VideoWriter('../Figs/my_movie.avi'); %create the video object
    video.FrameRate = 2;
    open(video); %open the file for writing

for ii = 1:size(psfx,3)

    figure(1), imagesc(log10(psfx(:,:,ii)), [-12 -7]), axis xy image, colormap(jet), colorbar, niceaxes
        title(sprintf('itr = %i: Cb = %0.3g (input-rms = 30nm)', ii, cbx(ii+1))), drawnow

        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(video,frame);

end % for ii - loop
      
    close(video); %close the file

return
%end % -----------------------------------------------
%%
x = [0 5 10 15 20];
%cb = [1.74 3.45 3.59 5.09 10.6] * 1e-11;
cb = [9.33e-2 1.06 1.4 1.05 0.995] * 1e-10;
tc = [28.69 13.24 13.03 7.68 6.91];

px = polyfit(x,cb,2);
y = polyval(px, x);

figure(2), clf, plot(x,cb, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('Bandwidth [%]')
    ylabel('Mean Contrast, Cb')
    title('EFC with Flattening: Error = 30nm')
    %legend('Data', '2nd-Order Fit', 'Location', 'Northwest')
    axis([0 20 0 2e-10])
print -dpng ../Figs/fig_cb

figure(3), clf, plot(x,tc, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('Bandwidth [%]')
    ylabel('Core Throughput, Tc')
    title('EFC with Flattening: Error = 30nm')
    %legend('Data', '2nd-Order Fit', 'Location', 'Northwest')
    %axis([0 20 0 2e-10])
print -dpng ../Figs/fig_tc

return
%end % -----------------------------------------------
%%
x = [0 10 20 30];
cb = [581] * 1e-10;
tc = [27.02 16.94 6.31 7.18];
ydm1 = [34.6 44.3 59.2 58.4];
ydm2 = [36.7 46.2 61.3 56.2];

x = [0 8 16 24];
cb = [5.81 2.93 29.9 102] * 1e-12;
tc = [34.91 34.87 32.44 25.62];
ydm1 = [4.9 8.35 17.09 34.1];
ydm2 = [5.45 9.15 17.9 36.4];
ydm = sqrt(ydm1.^2 + ydm2.^2);

x = [0 8 16 24 25];
cb = [0.763 0.78 1.21 1.46 1.06] * 1e-10;
tc = [8.22 7.7 7.82 9.56 9.98];
ydm1 = [61.44 78.9 78.2 75.6 55.8];
ydm2 = [68 84.4 83.6 80.8 58.1];
ydm = sqrt(ydm1.^2 + ydm2.^2);


tx = 'Segment Gap Width [mm]';

figure(1), clf, semilogy(x,cb, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel(tx)
    ylabel('Mean Contrast, Cb')
    title('Cb: EFC After Flattening')
    %axis([0 20 0 2e-10])
print -dpng ../Figs/fig_cb

figure(2), clf, plot(x,tc, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel(tx)
    ylabel('Core Throughput, Tc')
    title('Tc: EFC After Flattening')
print -dpng ../Figs/fig_tc

px = polyfit(x,ydm,2);
y2 = polyval(px, x);

%figure(3), clf, plot(x,ydm, 'ro-', x, y2, 'b--'), grid, 
figure(3), clf, plot(x,ydm, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel(tx)
    ylabel('RSS of DM Command RMS')
    title('DM-RMS: EFC After Flattening')
    %legend('Data', sprintf('Fit: y = %0.2fx^2 + (%0.1f)x + %0.1f',px), 'Location', 'Northwest')
print -dpng ../Figs/fig_rms

ns = length(ydm);
kk = 1:ns-1;
[xs, hh] = sort(ydm(kk));
tc1 = tc(kk);
y1 = tc1(hh);
px = polyfit(xs,y1,1);
y2 = polyval(px, xs);

%figure(4), clf, plot(xs,y1, 'ro-', xs, y2, 'b--'), grid, 
figure(4), clf, plot(xs,y1, 'ro-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('RSS of DM Command RMS, Sorted')
    ylabel('Core Throughput, Tc')
    %ylabel('Mean Contrast, Cb')
    title('Tc vs DM Command RMS')
    %legend('Data', sprintf('Linear-Fit: y = %0.2fx + %0.1f',px), 'Location', 'Southwest')
print -dpng ../Figs/fig_fom


return
%end % -----------------------------------------------
%%
load dat_dm_err1/opd10_erkin dopd_c_nm
maskz = dopd_c_nm ~= 0;
zk = [0 0 0 randn(1,12)];
opdz = zernike_compose(maskz, zk);
rmsz = stat2d(opdz);
fac = rmsz(1);
rmsz = rmsz / fac;
opdz = opdz / fac;

    aa = opdz; rms0 = rmsz; cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('random z4-z15: rms = %0.2f, pv = %0.1f', rmsz));
   print -dpng ../Figs/fig_dm1

return
%end % -----------------------------------------------
%%

pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
mask1 = opdnm ~= 0;


figure(1), imagesc(2*mask1-pupil), axis xy image, colormap(gray), %colorbar, %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('white = opd-mask, gray = pupil')
print -dpng ../Figs/fig_psf

figure(2), imagesc(opdnm), axis xy image, colormap(jet), colorbar, %niceaxes
return
%end % -----------------------------------------------
%%


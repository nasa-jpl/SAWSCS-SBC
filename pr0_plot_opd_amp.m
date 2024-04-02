% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
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

load dat_dm_err1/target_opd_amp_gap08_bb_30nm_bw20_355nm_unwrapped lam_um mask unwrapped_nm unwrapped_rad
%whos, return


load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x
dm = [dm1x dm2x];
rms0 = stat2d(dm);
a1 = rms0(1) * 0.3;
a2 = rms0(1) * 0.6;
a3 = rms0(1) * 0.9;
%pp = [a1 a2 a3]

idd = 700

if 0

    %save ~esidick/Afalco/falco20200916/macos/maps_psd/post_efc_opd_amp_30nm_bw20_1500nm_gray08_owa12_big opd_target amp_target ce_target
    %load ~esidick/Afalco/falco20200916/macos/maps_psd/post_efc_opd_amp_30nm_bw20_1500nm_gray08_owa12 opd_target amp_target ce_target   
    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_500nm opd_target amp_target ce_target
    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_500nm_post_efc opd_target amp_target ce_target
    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target
    load IFopd/target_opd_amp_gap08_bb_30nm_bw20_500nm_post_efc_no_pm opd_target amp_target ce_target
    amp = amp_target;
    %opd = opd_target;
    opd = unwrapped_nm .* (amp~=0);

%opd = mask .* real((angle(ce) - angle(ce1))* lam0 / 2 / pi);
%[opdf, zk] = zernike_remove(opd, [1]);
%rmsa = stat2d(opdf); cx = rmsa(1)*3*[-1 1];

%nx = [-np/2:np/2-1];
%[xc,yc] = fun_circle(nx, np/2);

aa = amp; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(1), clf, imagesc(aa);  axis xy image; colormap(jet), niceaxes
    %hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('amp: no pm, \lambda = 1500nm, \Delta\lambda/\lambda = 0.2'), 
    xlabel(sprintf('rms = %0.3g, pv = %0.1f', rms0)), 
print -dpng ../Figs/fig_amp

aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(2), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('unwrapped phase: \lambda = 355nm, \Delta\lambda/\lambda = 0.2'), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
print -dpng ../Figs/fig_opd

return
%end % -----------------------------------------------
%%
idd = 1

    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_500nm_post_efc_no_pm opd_target amp_target ce_target
    %load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc_no_pm opd_target amp_target ce_target
    load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target
    amp = amp_target;
    opd0 = opd_target;

    mask = amp_target ~= 0;

    lam = 1500;

    re = real(ce_target);
    im = imag(ce_target); %figure(4), clf, imagesc([re im]), axis image, colorbar
    facx = [0.1:0.1:1];
    ns = length(facx);

    rmsx = [];
    opdx = [];

for ii = 1:ns
    fac = facx(ii);

    opd = mask .* (atan2(im*fac, re) * lam / 2 / pi);
    rmsa = stat2d(opd); 
    rmsx = [rmsx; rmsa];
    opdx = [opdx opd];
end

%figure(3), clf, plot(facx, rmsx(:,1)'/rmsx(end,1), 'ro-'), grid
figure(3), hold on, plot(facx, rmsx(:,1)'/rmsx(end,1), 'c^-'), hold off

figure(4), clf, imagesc(opdx), axis image, colorbar
return

%nx = [-np/2:np/2-1];
%[xc,yc] = fun_circle(nx, np/2);

aa = amp; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(1), clf, imagesc(aa);  axis xy image; colormap(jet), niceaxes
    %hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('amp: no pm, \lambda = 1500nm, \Delta\lambda/\lambda = 0.2'), 
    xlabel(sprintf('rms = %0.3g, pv = %0.1f', rms0)), 
%print -dpng ../Figs/fig_amp

aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(2), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('unwrapped phase: \lambda = 355nm, \Delta\lambda/\lambda = 0.2'), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
%print -dpng ../Figs/fig_opd

return
%end % -----------------------------------------------
%%

load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x

aa = dm1x; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(1), clf, imagesc(aa);  axis xy image; colormap(jet), %niceaxes
    %hold on, plot(xc, yc, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('dm1 drift'), 
    xlabel(sprintf('rms = %0.3g, pv = %0.1fnm', rms0)), 
print -dpng ../Figs/fig_dm1


aa = dm2x; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(2), clf, imagesc(aa);  axis xy image; colormap(jet), %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('dm2 drift'), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
print -dpng ../Figs/fig_dm2

aa = [dm1x dm2x]; rms0 = stat2d(aa); 
mx = [0.1 0.2 0.3 0.6 0.9] * rms0(1)


return
%end % -----------------------------------------------
%%
%load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm           opd_target amp_target ce_target
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
a = max(amp_target(:))
ns = 21;

load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_sig_01nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y1 = 1e3*rmsx(:,1)'; x1 = [1:length(y1)]-1; z1 = fun_map_rms(ampx/a, ns); 
load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_sig_02nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y2 = 1e3*rmsx(:,1)'; x2 = [1:length(y2)]-1; z2 = fun_map_rms(ampx/a, ns);
load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_sig_03nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y3 = 1e3*rmsx(1:21,1)'; x3 = [1:length(y3)]-1; z3 = fun_map_rms(ampx/a, ns);

load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_01nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y1a = 1e3*rmsx(:,1)'; z1a = fun_map_rms(ampx/a, ns);
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_02nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y2a = 1e3*rmsx(:,1)'; z2a = fun_map_rms(ampx/a, ns);
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_03nm_itr20 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y3a = 1e3*rmsx(:,1)'; z3a = fun_map_rms(ampx/a, ns);

s1 = '\gamma = 0.1';
s2 = '\gamma = 0.2';
s3 = '\gamma = 0.3';

s1 = sprintf('dm drift rms = %ipm', aa1);
s2 = sprintf('dm drift rms = %ipm', aa2);
s3 = sprintf('dm drift rms = %ipm', aa3);

figure(1), clf, semilogy(x1,y1a, 'ro-', x2,y2a, 'bs-', x3,y3a, 'gd-'), grid, 
hold on, semilogy(x1,y1, 'ro--', x2,y2, 'bs--', x3,y3, 'gd--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of residual wfe [pm]')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3, 'Location', 'Northeast')
print -dpng ../Figs/fig_opd2

figure(2), clf, semilogy(x1,z1a, 'ro-', x2,z2a, 'bs-', x3,z3a, 'gd-'), grid, 
hold on, semilogy(x1,z1, 'ro--', x2,z2, 'bs--', x3,z3, 'gd--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3, 'Location', 'Northeast')
print -dpng ../Figs/fig_amp2

return
%end % -----------------------------------------------
%%
%load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
a = max(amp_target(:))
ns = 21;

%load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_opdrb_200pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_1500nm_m3_050pm_itr20 rms_opd rms_amp
    y1 = 1e3*rms_opd(:,1)'; x1 = [1:length(y1)]-1; z1 = rms_amp(:,1)/a; 
%load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_opdrb_400pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_1500nm_m3_100pm_itr20 rms_opd rms_amp
    y2 = 1e3*rms_opd(:,1)'; x2 = [1:length(y2)]-1; z2 = rms_amp(:,1)/a; 
%load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_355nm_m2_opdrb_600pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_gap08_bb_30nm_bw20_1500nm_m3_150pm_itr20 rms_opd rms_amp
    y3 = 1e3*rms_opd(:,1)'; x3 = [1:length(y3)]-1; z3 = rms_amp(:,1)/a; 

%load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_opdrb_200pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_050pm_itr20 rms_opd rms_amp
    y1a = 1e3*rms_opd(:,1)'; z1a = rms_amp(:,1)/a; 
%load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_opdrb_400pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_100pm_itr20 rms_opd rms_amp
    y2a = 1e3*rms_opd(:,1)'; z2a = rms_amp(:,1)/a; 
%load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_opdrb_600pm_itr40 rms_opd rms_amp
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_1500nm_m3_150pm_itr20 rms_opd rms_amp
    y3a = 1e3*rms_opd(:,1)'; z3a = rms_amp(:,1)/a; 

s1 = sprintf('dm drift rms = 50pm');
s2 = sprintf('dm drift rms = 100pm');
s3 = sprintf('dm drift rms = 150pm');

%s1 = sprintf('segment rb drift rms = 200pm');
%s2 = sprintf('segment rb drift rms = 400pm');
%s3 = sprintf('segment rb drift rms = 600pm');


figure(1), clf, semilogy(x1,y1a, 'ro-', x2,y2a, 'bs-', x3,y3a, 'gd-'), grid, 
hold on, semilogy(x1,y1, 'ro--', x2,y2, 'bs--', x3,y3, 'gd--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of residual wfe [pm]')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3, 'Location', 'Northeast')
print -dpng ../Figs/fig_opd

figure(2), clf, semilogy(x1,z1a, 'ro-', x2,z2a, 'bs-', x3,z3a, 'gd-'), grid, 
hold on, semilogy(x1,z1, 'ro--', x2,z2, 'bs--', x3,z3, 'gd--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3, 'Location', 'Northeast')
    axis([0 20 1e-6 1e-3])
print -dpng ../Figs/fig_amp

return
%end % -----------------------------------------------
%%
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_355nm opd_target amp_target ce_target
a = max(amp_target(:))
ns = 41;

load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_03nm_itr40 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y1 = 1e3*rmsx(:,1)'; x1 = [1:length(y1)]-1; z1 = fun_map_rms(ampx/a, ns); 
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_06nm_itr40 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y2 = 1e3*rmsx(:,1)'; x2 = [1:length(y2)]-1; z2 = fun_map_rms(ampx/a, ns);
load dat_dm_err1/wfc_dm1_dm2_gap08_bb_30nm_bw20_355nm_m3_sig_09nm_itr40 rmsx ampx %dm1x dm2x opdx ampx rmsx dmrms
    y3 = 1e3*rmsx(:,1)'; x3 = [1:length(y3)]-1; z3 = fun_map_rms(ampx/a, ns);

s1 = '\gamma = 0.3';
s2 = '\gamma = 0.6';
s3 = '\gamma = 0.9';

aa1 = 180; aa2 = aa1*6/3; aa3 = aa1*9/3;
s1 = sprintf('dm drift rms = %ipm', aa1);
s2 = sprintf('dm drift rms = %ipm', aa2);
s3 = sprintf('dm drift rms = %ipm', aa3);

figure(1), clf, semilogy(x1,y1, 'ro-', x2,y2, 'bs-', x3,y3, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of residual wfe [pm]')
    title('\Deltaphase vs wfc iterations')
    legend(s1,s2,s3, 'Location', 'Northeast')
print -dpng ../Figs/fig_opd1

figure(2), clf, semilogy(x1,z1, 'ro-', x2,z2, 'bs-', x3,z3, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp')
    title('\Deltaamp vs wfc iterations')
    legend(s1,s2,s3, 'Location', 'Northeast')
print -dpng ../Figs/fig_amp1

return
%end % -----------------------------------------------
%%
    %save dat_dm_err1/post_wfc_gap08_30nm_bw20_500nm_12lam_sig09 mx cbx rms_opd rms_amp psfx psfn0

load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_050pm_itr20_1500nm mx cbx , cb3 = [cbx(end)];
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_100pm_itr20_1500nm mx cbx , cb3 = [cb3 cbx(end)];
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_150pm_itr20_1500nm mx cbx , cb3 = [cb3 cbx(end)];
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_500pm_itr20_1500nm mx cbx , cb3 = [cb3 cbx(end)];
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_1000pm_itr20_1500nm mx cbx , cb3 = [cb3 cbx(end)];
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_1500pm_itr20_1500nm mx cbx , cb3 = [cb3 cbx(end)];

load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_050pm_itr20_1500nm mx cbx , cb1 = [cbx(1)];     cb2 = [cbx(end)];
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_100pm_itr20_1500nm mx cbx , cb1 = [cb1 cbx(1)]; cb2 = [cb2 cbx(end)];
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_150pm_itr20_1500nm mx cbx , cb1 = [cb1 cbx(1)]; cb2 = [cb2 cbx(end)];

%load dat_dm_err1/post_wfc_gap08_30nm_bw20_500nm_12lam_sig03 mx cbx %rms_opd rms_amp psfx psfn0
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_500pm_itr20_1500nm mx cbx rms_opd rms_amp,  
    y1 = [cbx]; x1 = [1:length(y1)]-1; o1 = rms_opd(:,1); a1 = rms_amp(:,1); cb1 = [cb1 cbx(1)];     cb2 = [cb2 cbx(end)];
%load dat_dm_err1/post_wfc_gap08_30nm_bw20_500nm_12lam_sig06 mx cbx %rms_opd rms_amp psfx psfn0
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_1000pm_itr20_1500nm mx cbx rms_opd rms_amp, 
    y2 = [cbx]; x2 = [1:length(y2)]-1; o2 = rms_opd(:,1); a2 = rms_amp(:,1); cb1 = [cb1 cbx(1)];     cb2 = [cb2 cbx(end)];
%load dat_dm_err1/post_wfc_gap08_30nm_bw20_500nm_12lam_sig09 mx cbx %rms_opd rms_amp psfx psfn0
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_1500pm_itr20_1500nm mx cbx rms_opd rms_amp, 
    y3 = [cbx]; x3 = [1:length(y3)]-1; o3 = rms_opd(:,1); a3 = rms_amp(:,1); cb1 = [cb1 cbx(1)];     cb2 = [cb2 cbx(end)];

aa1 = 0.5; aa2 = 1; aa3 = 1.5;
s1 = sprintf('dm drift rms = %0.1fnm', aa1);
s2 = sprintf('dm drift rms = %0.1fnm', aa2);
s3 = sprintf('dm drift rms = %0.1fnm', aa3);

figure(3), clf, semilogy(x1,y1, 'ro-', x2,y2, 'bs-', x3,y3, 'gd-', x3, x3*0+7.8e-11, 'm--'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('mean broadband contrast, Cb')
    title('contrast: wfc with 1500nm')
    legend(s1,s2,s3, 'nominal', 'Location', 'Northeast')
    %axis([0 20 5e-11 5e-7])
%print -dpng ../Figs/fig_cb3

figure(1), clf, semilogy(x1,o1, 'ro-', x2,o2, 'bs-', x3,o3, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of residual wfe [pm]')
    title('\Deltaphase: wfc with 1500nm')
    legend(s1,s2,s3, 'Location', 'Northeast')
%print -dpng ../Figs/fig_opd3

figure(2), clf, semilogy(x1,a1, 'ro-', x2,a2, 'bs-', x3,a3, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp')
    title('\Deltaamp: wfc with 1500nm')
    legend(s1,s2,s3, 'Location', 'Northeast')
%print -dpng ../Figs/fig_amp3

x = [50 100 150 500 1000 1500];
x3 = x(1:3);

figure(4), clf, loglog(x,cb1, 'ro-', x, cb3, 'bs-', x, cb2, 'gd-', x, x*0+7.8e-11, 'm--'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('dm-drift rms [pm]')
    ylabel('mean broadband contrast, Cb')
    title('wfc with \lambda = 1500nm')
    legend('pre-wfc', 'dm1 wfc', 'dm1+dm2 wfc', 'nominal', 'Location', 'Northwest')
    axis([50 1500 5e-11 2e-6])
print -dpng ../Figs/fig_cb3a


return
end % -----------------------------------------------
%%
idd = 1
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_sig01_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_200pm_itr40 mx cbx
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_050pm_itr20_1500nm mx cbx
    y1 = [1.62e-9 cbx]; x1 = [1:length(y1)]-1; 
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_sig02_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_400pm_itr40 mx cbx
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_100pm_itr20_1500nm mx cbx
    y2 = [6.24e-9 cbx]; x2 = [1:length(y2)]-1; 
load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_600pm_itr40 mx cbx
%load dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_dm_150pm_itr20_1500nm mx cbx
    y3 = [1.39e-8 cbx]; x3 = [1:length(y3)]-1; 

%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_sig01_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_200pm_itr40 mx cbx
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_050pm_itr20_1500nm mx cbx
    y1a = [1.62e-9 cbx]; x1 = [1:length(y1)]-1; 
%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_sig02_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_400pm_itr40 mx cbx
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_100pm_itr20_1500nm mx cbx
    y2a = [6.24e-9 cbx]; x1 = [1:length(y1)]-1; 
%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_sig03_itr20 mx cbx %rms_opd rms_amp psfx psfn0
%load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_600pm_itr40 mx cbx
load dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_dm_150pm_itr20_1500nm mx cbx
    y3a = [1.39e-8 cbx]; x1 = [1:length(y1)]-1; 

aa1 = 50; aa2 = 100; aa3 = 150;

s1 = sprintf('dm drift rms = %ipm', aa1);
s2 = sprintf('dm drift rms = %ipm', aa2);
s3 = sprintf('dm drift rms = %ipm', aa3);

%s1 = sprintf('segment rb drift rms = 200pm');
%s2 = sprintf('segment rb drift rms = 400pm');
%s3 = sprintf('segment rb drift rms = 600pm');

figure(3), clf, semilogy(x1,y1a, 'ro-', x2,y2a, 'bs-', x3,y3a, 'gd-', x3, x3*0+7.8e-11, 'm-'), grid, 
    hold on, semilogy(x1,y1, 'ro--', x2,y2, 'bs--', x3,y3, 'gd--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('mean broadband contrast, Cb')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3, 'nominal', 'Location', 'Northeast')
    axis([0 20 5e-11 5e-8])
%print -dpng ../Figs/fig_cbx

return
%end % -----------------------------------------------
%%
xx = [0.1 1 3 10 30 100 300 1000 3000];
ns = length(xx);
mx = [0];
cb0 = 7.8e-11;
cb1 = cbi;
cb2 = cbi;

fn1{1} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_001pm_itr20_1500nm';
fn1{2} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_003pm_itr20_1500nm';
fn1{3} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_010pm_itr20_1500nm';
fn1{4} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_030pm_itr20_1500nm';
fn1{5} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_100pm_itr20_1500nm';
fn1{6} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_300pm_itr20_1500nm';
fn1{7} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_1000pm_itr20_1500nm';
fn1{8} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_3000pm_itr20_1500nm';

fn2{1} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_001pm_itr20';
fn2{2} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_003pm_itr20';
fn2{3} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_010pm_itr20';
fn2{4} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_030pm_itr20';
fn2{5} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_100pm_itr20';
fn2{6} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_300pm_itr20_1500nm';
fn2{7} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_1000pm_itr20_1500nm';
fn2{8} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_3000pm_itr20_1500nm';

load(fn2{1},'mx', 'cbx'); y1a = [cbx cbx(end)]; x1a = [1:length(y1a)]-1;  cbb2 = [cb0 cbx(end)];      cbi = [cb0 cbx(1)];
load(fn2{2},'mx', 'cbx'); y2a = [cbx cbx(end)]; x2a = [1:length(y2a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{3},'mx', 'cbx'); y3a = [cbx cbx(end)]; x3a = [1:length(y3a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{4},'mx', 'cbx'); y4a = [cbx cbx(end)]; x4a = [1:length(y4a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{5},'mx', 'cbx'); y5a = [cbx cbx(end)]; x5a = [1:length(y5a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{6},'mx', 'cbx'); y6a = [cbx cbx(end)]; x6a = [1:length(y6a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{7},'mx', 'cbx'); y7a = [cbx cbx(end)]; x7a = [1:length(y7a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];
load(fn2{8},'mx', 'cbx'); y8a = [cbx cbx(end)]; x8a = [1:length(y8a)]-1;  cbb2 = [cbb2 cbx(end)]; cbi = [cbi cbx(1)];

load(fn1{1},'mx', 'cbx'); y1 = [y1a(1) cbx]; x1 = [1:length(y1)]-1; cbb1 = [cb0 cbx(end)];
load(fn1{2},'mx', 'cbx'); y2 = [y2a(1) cbx]; x2 = [1:length(y2)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{3},'mx', 'cbx'); y3 = [y3a(1) cbx]; x3 = [1:length(y3)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{4},'mx', 'cbx'); y4 = [y4a(1) cbx]; x4 = [1:length(y4)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{5},'mx', 'cbx'); y5 = [y5a(1) cbx]; x5 = [1:length(y5)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{6},'mx', 'cbx'); y6 = [y6a(1) cbx]; x6 = [1:length(y6)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{7},'mx', 'cbx'); y7 = [y7a(1) cbx]; x5 = [1:length(y7)]-1; cbb1 = [cbb1 cbx(end)];
load(fn1{8},'mx', 'cbx'); y8 = [y8a(1) cbx]; x5 = [1:length(y8)]-1; cbb1 = [cbb1 cbx(end)];

if 0
aa1 = 1; aa2 = 3; aa3 = 10; aa4 = 30; aa5 = 100;

s1 = sprintf('segment rb drift rms = %ipm', aa1);
s2 = sprintf('segment rb drift rms = %ipm', aa2);
s3 = sprintf('segment rb drift rms = %ipm', aa3);
s4 = sprintf('segment rb drift rms = %ipm', aa4);
s5 = sprintf('segment rb drift rms = %ipm', aa5);

figure(3), clf, semilogy(x1a,y1a, 'ro-', x2a,y2a, 'bs-', x3a,y3a, 'gd-', x4a,y4a, 'm^-', x5a,y5a, 'c>-', x3a, x3a*0+7.8e-11, 'k-'), grid, 
    hold on, semilogy(x1,y1, 'ro--', x2,y2, 'bs--', x3,y3, 'gd--', x4,y4, 'm^--', x5,y5, 'c>--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('mean broadband contrast, Cb')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3,s4,s5, 'nominal', 'Location', 'Northeast')
    axis([0 20 5e-11 1e-9])
print -dpng ../Figs/fig_cbx1

x = [1 3 10 30 100];

end

figure(2), clf, loglog(xx,cbi, 'ro-', xx, cbb1, 'bs-', xx, cbb2, 'gd-'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('segment rb-drift rms [pm]')
    ylabel('mean broadband contrast, Cb')
    title('wfc with \lambda = 1500nm')
    legend('pre-wfc', 'dm1 wfc', 'dm1+dm2 wfc',  'Location', 'Northwest')
    axis([0.1 3000 5e-11 1e-6])
print -dpng ../Figs/fig_cbx2b

return
%end % -----------------------------------------------
%%
load IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc opd_target amp_target ce_target
a = max(amp_target(:))
ns = 21;

fn1{1} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_001pm_itr20_1500nm';
fn1{2} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_003pm_itr20_1500nm';
fn1{3} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_010pm_itr20_1500nm';
fn1{4} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_030pm_itr20_1500nm';
fn1{5} = 'dat_dm_err1/post_dm1_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_100pm_itr20_1500nm';

fn2{1} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_001pm_itr20';
fn2{2} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_003pm_itr20';
fn2{3} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_010pm_itr20';
fn2{4} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_030pm_itr20';
fn2{5} = 'dat_dm_err1/post_dm1_dm2_wfc_gap08_30nm_bw20_500nm_12lam_opdrb_100pm_itr20';

load(fn1{1},'rms_opd', 'rms_amp'); y1 = 1e3*rms_opd(:,1)'; x1 = [1:length(y1)]-1; z1 = rms_amp(:,1)/a;
load(fn1{2},'rms_opd', 'rms_amp'); y2 = 1e3*rms_opd(:,1)'; x2 = [1:length(y1)]-1; z2 = rms_amp(:,1)/a;
load(fn1{3},'rms_opd', 'rms_amp'); y3 = 1e3*rms_opd(:,1)'; x3 = [1:length(y1)]-1; z3 = rms_amp(:,1)/a;
load(fn1{4},'rms_opd', 'rms_amp'); y4 = 1e3*rms_opd(:,1)'; x4 = [1:length(y1)]-1; z4 = rms_amp(:,1)/a;
load(fn1{5},'rms_opd', 'rms_amp'); y5 = 1e3*rms_opd(:,1)'; x5 = [1:length(y1)]-1; z5 = rms_amp(:,1)/a;

load(fn2{1},'rms_opd', 'rms_amp'); y1a = 1e3*rms_opd(:,1)'; x1a = [1:length(y1a)]-1; z1a = rms_amp(:,1)/a;
load(fn2{2},'rms_opd', 'rms_amp'); y2a = 1e3*rms_opd(:,1)'; x2a = [1:length(y2a)]-1; z2a = rms_amp(:,1)/a;
load(fn2{3},'rms_opd', 'rms_amp'); y3a = 1e3*rms_opd(:,1)'; x3a = [1:length(y3a)]-1; z3a = rms_amp(:,1)/a;
load(fn2{4},'rms_opd', 'rms_amp'); y4a = 1e3*rms_opd(:,1)'; x4a = [1:length(y4a)]-1; z4a = rms_amp(:,1)/a;
load(fn2{5},'rms_opd', 'rms_amp'); y5a = 1e3*rms_opd(:,1)'; x5a = [1:length(y5a)]-1; z5a = rms_amp(:,1)/a;

aa1 = 1; aa2 = 3; aa3 = 10; aa4 = 30; aa5 = 100;

s1 = sprintf('segment rb drift rms = %ipm', aa1);
s2 = sprintf('segment rb drift rms = %ipm', aa2);
s3 = sprintf('segment rb drift rms = %ipm', aa3);
s4 = sprintf('segment rb drift rms = %ipm', aa4);
s5 = sprintf('segment rb drift rms = %ipm', aa5);

figure(1), clf, semilogy(x1a,y1a, 'ro-', x2a,y2a, 'bs-', x3a,y3a, 'gd-',x4a,y4a, 'm^-', x5a,y5a, 'c>-'), grid, 
hold on, semilogy(x1,y1, 'ro--', x2,y2, 'bs--', x3,y3, 'gd--',x4,y4, 'm^--', x5,y5, 'c>--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of residual wfe [pm]')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3,s4,s5, 'Location', 'Northeast')
print -dpng ../Figs/fig_opd1

figure(2), clf, semilogy(x1a,z1a, 'ro-', x2a,z2a, 'bs-', x3a,z3a, 'gd-', x4a,z4a, 'm^-', x5a,z5a, 'c>-'), grid, 
hold on, semilogy(x1,z1, 'ro--', x2,z2, 'bs--', x3,z3, 'gd--', x4,z4, 'm^--', x5,z5, 'c>--'), hold off, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    ylabel('rms of \Deltaamp')
    title('solid: dm1+dm2 wfc. dashed: dm1 wfc')
    legend(s1,s2,s3,s4,s5, 'Location', 'Northeast')
    %axis([0 20 1e-6 1e-3])
print -dpng ../Figs/fig_amp1

return
%end % -----------------------------------------------
%%
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
    opd = opdnm;

aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(2), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('wavefront phase: rms = %0.1fnm', rms0(1))), 
print -dpng ../Figs/fig_opd

return
%end % -----------------------------------------------
%%

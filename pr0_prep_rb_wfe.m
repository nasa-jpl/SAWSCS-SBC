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

id1 = 2

if 0

%load /proj/jwst2/6MST/data/jzlou/dwdx_6mst_urad_um_10042023.mat OPDMask_g dwdx indx wnom
%load /proj/jwst2/6MST/data/jzlou/dwdx_6mst_noSegGap_urad_um_10122023.mat
%load IFopd/G_dm1_lim3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx
load /proj/jwst2/6MST/data/jzlou/dwdx_6mst_noSegGap_urad_um_10132023.mat
%whos, return


[nr, nact] = size(dwdx);

%wx = dwdx(:,7);
drb = randn(nact,1);
wx  = dwdx * drb;
opd = v2m(wx,indx);
rms0 = stat2d(opd);
opd = 0.1 * opd / rms0(1);
rms0 = stat2d(opd);

[nr, nc] = size(opd);
nh = (nr-1)/2;
nx = [-nh:nh];
[xc, yc] = fun_circle(nx, nh);

figure(1), clf, imagesc(nx,nx,opd), axis xy image, colormap(jet), niceaxes
    hold on, plot(xc, yc, 'b'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('rb-drift: rms = %0.1f, pv = %0.1fnm', rms0)), 
    %xlabel(sprintf('rms = %0.3g, pv = %0.3g', stat2d(amp))), 
%print -dpng ../Figs/fig_opd

%return
end % -----------------------------------------------
%%
figure(1), clf, 

clear opdrb_nm

rms0 = 1; % nm

for ii = 1:10

    drb = randn(nact,1);
    wx  = dwdx * drb;
    opd = v2m(wx,indx);
    rmsi = stat2d(opd);
    opd = rms0 * opd / rmsi(1);
    rmsi = stat2d(opd);

    opdrb_nm(:,:,ii) = opd;

    subplot(3,4,ii), imagesc(opd), axis xy image, colormap(jet), niceaxes
    title(sprintf('%i: rms = %0.1fnm', ii, rmsi(1))), 
end
print -dpng ../Figs/fig_opd1

save dat_rb_wfe/opdrb_1000pm opdrb_nm

return
%end % -----------------------------------------------
%%

load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  

rms0 = stat2d(opdnm);

figure(1), clf, imagesc(nx,nx,opdnm, rms0(1)*3*[-1 1]), axis xy image, colormap(jet), niceaxes
    hold on, plot(xc, yc, 'b'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('rb-drift: rms = %0.1f, pv = %0.1fnm', rms0)), 
    %xlabel(sprintf('rms = %0.3g, pv = %0.3g', stat2d(amp))), 
print -dpng ../Figs/fig_opd


return
%end % -----------------------------------------------
%%

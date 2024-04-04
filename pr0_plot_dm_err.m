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

if 0

    load dat_dm_err/dm_err_1000pm dm1x_3d dm2x_3d

    aa = dm1x_3d(:,:,1); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = dm2x_3d(:,:,1); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2


return
end % -----------------------------------------------
%%
if 0

    load dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_10pm_run10 opdx opdi

    aa = opdi; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('pre-wfc: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = opdx; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-wfc: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2


return
end % -----------------------------------------------
%%

if 0

    Colx = ['ro-'; 'bs-'; 'gd-'; 'c^-'; 'm>-'; 'ko-'];

    runx = [1:21]'-1;
    drift = [10 30 100 300 1000 3000];
    ns = length(drift);

    fn = 'dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_%ipm_run1';

figure(1), clf
clear legx

for ii = 1:ns
    %load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp','dmrms');
    fni = sprintf(fn, drift(ii))
    load(fni, 'rms_opd');
    legx{ii} = sprintf('drift = %ipm', drift(ii));
    semilogy(runx, rms_opd(:,1)*1e3, Colx(ii,:)), hold on
end
hold off, grid
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iteration number')
    ylabel('rms of residual wfe [pm]')
    title('wfc with \lambda = 1500nm')
    legend(legx,  'Location', 'Northeast')
    %axis([0.1 3000 5e-11 1e-6])
print -dpng ../Figs/fig_rms

return
end % -----------------------------------------------
%%
if 0

runx = [1:10];
fn = 'dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_100pm_run%i';

rmsi = [];
rmso = [];

for ii = 1:length(runx)
    %load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rms_opd', 'rms_amp','dmrms');
    load(sprintf(fn, ii'), 'rms_opd', 'opdx');
    rmsi = [rmsi rms_opd(1,1)];
    rmso = [rmso rms_opd(end,1)];
end

%s1 = sprintf('segment rb drift rms = %ipm', aa1);

figure(1), clf, semilogy(runx,rmsi, 'ro-', runx, rmso, 'ro--'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('random error case number')
    ylabel('rms of differential wfe [pm]')
    title('wfc with \lambda = 1500nm')
    legend('pre-wfc', 'dm1+dm2 wfc',  'Location', 'Northwest')
    %axis([0.1 3000 5e-11 1e-6])
%print -dpng ../Figs/fig_cbx2b

figure(2), clf, aa = opdx(:,:,1)*1e3; imagesc(aa), axis image, colormap(jet), colorbar, title(sprintf('rms = %0.1f, pv = %0.1fpm', stat2d(aa)))
%figure(3), clf, aa = opdx(:,:,end); imagesc(aa), axis image, colorbar, title(sprintf('rms = %0.1f, pv = %0.1f', stat2d(aa)))

return
end % -----------------------------------------------
%%
load ~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_3000pm_run1_m5 rms_opd
    y1 = rms_opd;
load ~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_3000pm_run1 rms_opd
    runx = [1:21]'-1;

figure(1), clf    
    semilogy(runx, rms_opd(:,1)*1e3, 'ro-' , runx, y1(:,1)*1e3, 'bs-'), grid
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iteration number')
    ylabel('rms of residual wfe [pm]')
    title('wfc: \lambda = 1500nm, case of drift = 3nm')
    legend('act-reg factor \beta = -3', 'act-reg factor \beta = -5',  'Location', 'Northeast')
    %axis([0.1 3000 5e-11 1e-6])
print -dpng ../Figs/fig_rms
    
return
%end % -----------------------------------------------
%%

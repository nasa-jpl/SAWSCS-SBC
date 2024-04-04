% ------------------------------------------------------------------------------
% pr_psd_module.m
% Down-samples a surface map using its Zernike coefficients and
% a PSD model obtained from the surface map residual
% PSD Reference: E. Sidick, SPIE vol. 7390, paper 73900L (2009)
%
% Erkin Sidick. 3/22/2010
% Erkin Sidick. Modified on 5/26/2023
% ------------------------------------------------------------------------------

addpath('/home/esidick/pub/psd_module/util/');
addpath(genpath('/proj/exospin/users/joegreen/'));
    
% -----------------------------------------------------------------------
% Load opdo(rig), and divide it into 2 parts: opdo = opdf(it) + opdr(es)
% -----------------------------------------------------------------------
if 1
    % Load input OPD map:
    %load maps_psd/target_opd_amp_bb_00nm_bw20_500nm_fac0p5_gray24_sharp opd_target %amp_target ce_target
    %load maps_psd/target_opd_amp_bb_00nm_bw20_500nm_fac0p5_gray24_sharp opd_target %amp_target ce_target
    %load maps_psd/target_opd_amp_bb_00nm_bw20_1500nm_gray08_all_owa20 opd_target %amp_target ce_target
    %opdo = opd_target;

%load maps_psd/target_opd_amp_bb_00nm_bw20_1500nm_gray08_opd_owa12 opd_target %amp_target ce_target
%load maps_psd/target_opd_amp_bb_00nm_bw20_1500nm_gray08_dm_owa12 opd_target %amp_target ce_target
%load maps_psd/target_opd_amp_bb_00nm_bw20_1500nm_gray08_all_owa12 opd_target %amp_target ce_target
%load maps_psd/pos_efc_nom_opd_amp_gray_gap08_wfe30_1lam_1500nm opd_target %amp_target ce_target
load maps_psd/pos_efc_opdhi_opd_amp_gray_gap08_wfe30_1lam_1500nm opd_target %amp_target ce_target
%load maps_psd/pos_efc_flat_opd_amp_gray_gap08_wfe30_1lam_1500nm opd_target %amp_target ce_target
opd = opd_target;
%load ~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_flat_30nm_500nm_gray_gap08 opdx %dm1x opdx
%opd = opdx(:,:,1);

bandwidth = 32;
rolloff = 0;

idd = 1

% [opd_res, filter] = opdFilter(opd, bandwidth, rolloff)
%[opdlo, filter] = opdFilter(opd, bandwidth, rolloff);
[opdlo, filter] = opdFilter_Rect(opd, bandwidth, rolloff);
opdhi = opd - opdlo;

opdo = opd;
    
    rms0 = stat2d(opdo); cx = rms0(1) * 3 * [-1 1];
    figure(1), imagesc(opdo, cx), axis image, colormap(jet), %colorbar,
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-efc-1500nm: filtered, total'))
        xlabel(sprintf('rms = %0.1f PV = %0.1f nm',rms0))
    print -dpng ../Figs/fig_opd2

    rms0 = stat2d(opdlo); cx = rms0(1) * 3 * [-1 1];
    figure(10), imagesc(opdlo, cx), axis image, colormap(jet), %colorbar,
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-efc-1500nm: filtered, opd-lo'))
        xlabel(sprintf('rms = %0.1f PV = %0.1f nm',rms0))
    print -dpng ../Figs/fig_opd2a

    rms0 = stat2d(opdhi); cx = rms0(1) * 3 * [-1 1];
    figure(11), imagesc(opdhi, cx), axis image, colormap(jet), %colorbar,
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('post-efc-1500nm: filtered, opd-hi'))
        xlabel(sprintf('rms = %0.1f PV = %0.1f nm',rms0))
    print -dpng ../Figs/fig_opd2b

        %return

    % Remove extra zeros:
    [xmin, xmax, ymin, ymax] = find_boundary(opdo~=0, 0.5);
    opdo = opdo(ymin:ymax, xmin:xmax);
    [nr, nc] = size(opdo);
    Mgrid = min([nr nc]);

    % Make Mgrid odd if not already
    if mod(Mgrid, 2) == 0, Mgrid = Mgrid + 1; end

    opdo = pad(opdo, Mgrid);

if 0
	% Remove low-order Zernikes if needed
    nzern = 200;
    masko = opdo ~= 0;
    [opdf, zcoef] = zernike_fit(opdo, masko, [1:nzern]);
    opdr  = opdo; % - opdf;
end

% -----------------------------------------------------------------------
% Calculate the 2-PSD of the surface map residual
% -----------------------------------------------------------------------
if 1

    
 npad = Mgrid;
 opdr = opdo;

 ma = 64;
 dp = 400e-6; % act pitch in m
 
    % Big Circular Mask
    npix  = [-(Mgrid-1)/2:(Mgrid-1)/2];
    [x y] = meshgrid(npix);
    Mask  = sqrt(x.^2 + y.^2) < Mgrid/2;
    
    % PSD parameters
    Mcen = floor(Mgrid/2) + 1;      % Center of the gridded mask
    D  = ma * dp;               % Surface map diameter, mm
	dx = D / Mgrid;         % mm per pixel of mask;
    f0 = 1 / D;             % Spatial-frequency normalization factor, cyc/mm 
    du = 1 / D;             % Spatial frequency increment, cyc/mm
    A  = D * D;             % Square area of the surface map, mm^2

    % 2D-PSD:
    psd_2d = abs(fftshift(fft2((opdr)))) * dx * dx;
    psd_2d = psd_2d.^2 / A;       % nm^2 mm^2
    
%if 1    
    
% -----------------------------------------------------------------------
% Calculate the corresponding radially-averaged 1D-PSD
% -----------------------------------------------------------------------
    % Parameters for radially-averaged (1D) PSD
    [M N] = size(Mask);
    [v u] = meshgrid (1:N,1:M); 
    u = (u - (floor(N/2) + 1))/N;
    v = (v - (floor(M/2) + 1))/M;
    rf = sqrt (u.*u + v.*v) / dx;   % Radial frequency, cyc/mm

    % Sorted radial-frequency and 2D-PSD
    rn  = rf / f0;         % Normalized frequency
    [rs hs] = sort(rn(:));
    psd0   = psd_2d(:);
    psds  = psd0(hs);
    
    % For radially-averaged quantities:
    Nr  = round(Mgrid/16);  % Total # of 1D PSD frequency
    ri = logspace(log10(rs(2)), log10(rs(end)), Nr-1);
    
    % Obtain 1D-PSD by radial-average
    rf_1d  = [];
    psd_1d = [];
    for ii = 1:length(ri)-1
        kk = find(rs > ri(ii) & rs <= ri(ii+1));
        if ~isempty(kk)
            rf_1d  = [rf_1d mean(rs(kk))];
            psd_1d = [psd_1d mean(psds(kk))];
        end
    end

    xs = rs / Mgrid;
    xv = rf_1d / Mgrid;

    xx = [-128:127];

    xa = [-1 1 1 -1 -1] * 32;
    ya = [-1 -1 1 1 -1] * 32;

    %save maps_psd/psd_gap08_owa20_1500nm rf_1d psd_1d  

figure(2), clf
    loglog(rs, psds, 'r.', rf_1d, psd_1d, 'b'), grid
    parms = fun_newaxes(16, 2, 0);
    xlabel('spatial-freq (cyc/D)')
    ylabel('PSD (nm^2 m^2)')
    title('2D-PSD & Its Radial-Average') 
    axis([1 rs(end) 1e-10 1e-2])
    legend('2D-PSD_t_r_u_e', '1D-PSD_t_r_u_e', 'Location', 'Southwest')
%print -dpng ../Figs/fig_psd

    figure(3),imagesc(xx, xx, log10([psd_2d]+eps), [-10 -2]), axis xy image, colormap(jet), %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        hold on, plot(xa, ya, 'k--', 'Linewidth', 2), hold off
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title('psd: opdhi, rect, bw = 32 cyc/D, rolloff = 0')
        title('post-efc psd: filtered, total, bw32, roll0')
        %title('2d-psd')
        xlabel('spatial-freq [cyc/D]')
        ylabel('spatial-freq [cyc/D]')
print -dpng ../Figs/fig_psd2

return
end
% -----------------------------------------------------------------------
% Find a PSD Model from the 1D-PSD curve. HPF = Half-power frequency
% Note: You must adjust the values of LB, UB & xvar_init as needed
% before calling fminsearch.m
% -----------------------------------------------------------------------
global rf_fit psd_fit

%addpath('/home/esidick/AFTA/mcb1/matprog/')

% Determine the fitting region PSD-freq pair
LB = 15;     % Fitting region lower-bound frequency, user-selected
UB = 128;   % Fitting region upper-bound frequency, user-selected
kk = find(rf_1d>LB & rf_1d < UB);
rf_fit  = rf_1d(kk);
psd_fit = psd_1d(kk);

% Carry out a Model fit to 1D-PSD in the selected region
[psdmax hmax] = max(psd_1d);
xvar_init = [psdmax rf_1d(hmax) 2.3];  % Max, HPF & Exponent

psd_tem = xvar_init(1) ./ (1 + power(rf_fit/xvar_init(2), xvar_init(3)) );

figure(4), clf
    loglog(rf_1d, psd_1d, 'ro-', rf_fit, psd_tem, 'bs-'), grid
    parms = fun_newaxes(16, 2, 0);
%return

% Find the Max, HPF & Exponent of the PSD model
    xvar_final = fminsearch('fun_findpsd_parms', xvar_init);  
    psdmax = xvar_final(1)
    hpf    = xvar_final(2)
    powerx = xvar_final(3)
    psdx = psdmax ./ (1 + power(rf_1d/hpf, powerx) );
    psd_model = psdx * 0;
    psd_model(kk) = psdx(kk);

mgrid = npad;

    legx{1} = sprintf('PSDmax = %0.3g', xvar_final(1));
    legx{2} = sprintf('HPF = %0.3f', xvar_final(2));
    legx{3} = sprintf('Power-Law = %0.3f', xvar_final(3));
    
    k1 = find(rf_1d >= LB & rf_1d <= UB);

%    save maps_psd/psd_gap08_owa12 rs rf_1d k1 psd_1d psd_model xvar_final 

figure(5), clf
    %loglog(rf_1d, psd_1d, 'ro-', rf_1d, psd_model+min(psd_1d), 'bs-'), grid
    loglog(rf_1d, psd_1d, 'ro-', rf_1d(k1), psd_model(k1), 'bs-'), grid
    parms = fun_newaxes(14, 2, 5);
    xlabel('spatial-freq (cyc/D)')
    ylabel('1D-PSD (nm^2 m^2)')
    title(sprintf('gap24, 1D-PSD: LB = %i, UB = %i Cyc/D',LB,UB)) 
    legend('1D-PSD_t_r_u_e', 'Fit', 'Location', 'Northeast')
    axis([1 150 1e-7 1e-1])
    text(1.2, 1e-4, legx, 'FontSize', 14)
%print -dpng ../Figs/fig_psd4a

return
%end % ----------------------------------------------------
%%

load maps_psd/psd_gap08_owa12_1500nm rf_1d psd_1d 
    m0 = 1; x1 = rf_1d; y1 = psd_1d; legx{m0} = 'owa = 12';
load maps_psd/psd_gap08_owa16_1500nm rf_1d psd_1d 
    m0 = m0 + 1; x2 = rf_1d; y2 = psd_1d; legx{m0} = 'owa = 16';
load maps_psd/psd_gap08_owa20_1500nm rf_1d psd_1d 
    m0 = m0 + 1; x3 = rf_1d; y3 = psd_1d; legx{m0} = 'owa = 20';

    m0 = m0 + 1; legx{m0} = 'x = 32';

figure(1), clf
    loglog(x1,y1,'r',x2,y2,'b',x3,y3,'g'), grid
    parms = fun_newaxes(14, 2, 0);
    xlabel('spatial-freq (cyc/aper)')
    ylabel('1D-PSD (nm^2 m^2)')
    title(sprintf('radially-averaged psd')) 
    %axis([1 150 1e-7 1e-1])
    xline(32, 'k--', 'x = 32', 'Linewidth', 2, 'FontSize', 14)
    legend(legx, 'Location', 'Southwest')
print -dpng ../Figs/fig_psd

return
%end % ----------------------------------------------------
%%

load maps_psd/psd_gap00 rs rf_1d k1 psd_1d psd_model xvar_final 
    m0 = 1; x1 = rf_1d; y1 = psd_1d; legx{m0} = 'gap = 00mm';
load maps_psd/psd_gap08 rs rf_1d k1 psd_1d psd_model xvar_final 
    m0 = m0 + 1; x2 = rf_1d; y2 = psd_1d; legx{m0} = 'gap = 08mm';
load maps_psd/psd_gap16 rs rf_1d k1 psd_1d psd_model xvar_final 
    m0 = m0 + 1; x3 = rf_1d; y3 = psd_1d; legx{m0} = 'gap = 16mm';
load maps_psd/psd_gap24 rs rf_1d k1 psd_1d psd_model xvar_final 
    m0 = m0 + 1; x4 = rf_1d; y4 = psd_1d; legx{m0} = 'gap = 24mm';
load maps_psd/psd_gap24_sharp rs rf_1d k1 psd_1d psd_model xvar_final 
    m0 = m0 + 1; x5 = rf_1d; y5 = psd_1d; legx{m0} = 'gap = 24mm, sharp';

    m0 = m0 + 1; legx{m0} = 'x = 32';

figure(1), clf
    loglog(x1,y1,'r',x2,y2,'b',x3,y3,'g',x4,y4,'m',x5,y5,'k'), grid
    parms = fun_newaxes(14, 2, 0);
    xlabel('spatial-freq (cyc/aper)')
    ylabel('1D-PSD (nm^2 m^2)')
    title(sprintf('radially-averaged psd')) 
    axis([1 150 1e-7 1e-1])
    xline(32, 'k--', 'x = 32', 'Linewidth', 2, 'FontSize', 14)
    legend(legx, 'Location', 'Southwest')
print -dpng ../Figs/fig_psd

return
end % ----------------------------------------------------
%%
load maps_psd/target_opd_amp_bb_00nm_bw20_1500nm_gray08_all_owa12 opd_target %amp_target ce_target

opd = opd_target;
bandwidth = 32;
rolloff = 0;

idd = 1

% [opd_res, filter] = opdFilter(opd, bandwidth, rolloff)
%[opdlo, filter] = opdFilter(opd, bandwidth, rolloff);
[opdlo, filter] = opdFilter_Rect(opd, bandwidth, rolloff);
opdhi = opd - opdlo;

aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
figure(1), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('opd-raw', bandwidth)), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
%print -dpng ../Figs/fig_opd1

aa = opdlo; rms0 = stat2d(aa); %cx = rms0(1)*3*[-1 1];
figure(2), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('opd-low: bandwidth = %i cyc/aper', bandwidth)), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
print -dpng ../Figs/fig_opd2a

aa = opdhi; rms0 = stat2d(aa); %cx = rms0(1)*3*[-1 1];
figure(3), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('opd-high: bandwidth = %i cyc/aper', bandwidth)), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
print -dpng ../Figs/fig_opd2b

return
%end % ----------------------------------------------------
%%


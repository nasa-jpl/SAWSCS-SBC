% ------------------------------------------------------------------------------
% pr_psd_module.m
% Down-samples a surface map using its Zernike coefficients and
% a PSD model obtained from the surface map residual
% PSD Reference: E. Sidick, SPIE vol. 7390, paper 73900L (2009)
%
% Erkin Sidick. 3/22/2010
% ------------------------------------------------------------------------------

addpath('/home/esidick/AMD/psd_module/util/');
addpath('/home/joegreen/matlab/');
    

% -----------------------------------------------------------------------
% Calculate the 2-PSD of the surface map residual
% -----------------------------------------------------------------------
if 0
    
 npad = Mgrid;
 opdr = opdo;
 
    % Big Circular Mask
    npix  = [-(Mgrid-1)/2:(Mgrid-1)/2];
    [x y] = meshgrid(npix);
    Mask  = sqrt(x.^2 + y.^2) < Mgrid/2;
    
    % PSD parameters
    Mcen = floor(Mgrid/2) + 1;      % Center of the gridded mask
    D  = 26.7326;               % Surface map diameter, mm
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

figure(5), clf
    loglog(rs, psds, 'r.', rf_1d, psd_1d, 'b'), grid
    parms = fun_newaxes(16, 2, 0);
    xlabel('Radial Frequency (Cycles/D)')
    ylabel('PSD (nm^2 m^2)')
    title('gpi-2: 2D-PSD & Its Radially-Avgd 1D-PSD') 
    %axis([3 rs(end) 1e-6 1e5])
    legend('2D-PSD_t_r_u_e', '1D-PSD_t_r_u_e', 1)
print -dtiff ../Figs/fig_psd
return
end
% -----------------------------------------------------------------------
% Find a PSD Model from the 1D-PSD curve. HPF = Half-power frequency
% Note: You must adjust the values of LB, UB & xvar_init as needed
% before calling fminsearch.m
% -----------------------------------------------------------------------
global rf_fit psd_fit

addpath('/home/esidick/AFTA/mcb1/matprog/')

% Determine the fitting region PSD-freq pair
LB = 4;     % Fitting region lower-bound frequency, user-selected
UB = 25;   % Fitting region upper-bound frequency, user-selected
kk = find(rf_1d>LB & rf_1d < UB);
rf_fit  = rf_1d(kk);
psd_fit = psd_1d(kk);

% Carry out a Model fit to 1D-PSD in the selected region
[psdmax hmax] = max(psd_1d);
xvar_init = [psdmax rf_1d(hmax) 2.3];  % Max, HPF & Exponent

psd_tem = xvar_init(1) ./ (1 + power(rf_fit/xvar_init(2), xvar_init(3)) );

figure(1), clf
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

%end
%save -v6 gpi_1_psd_parms psdmax hpf powerx LB UB D mgrid npad   

%end

    s1 = sprintf('PSDmax = %0.1f, HPF = %0.2f, Power-Law = %0.3f', xvar_final);
    
    k1 = find(rf_1d >= LB & rf_1d <= UB);

figure(2), clf
    %loglog(rf_1d, psd_1d, 'ro-', rf_1d, psd_model+min(psd_1d), 'bs-'), grid
    loglog(rf_1d, psd_1d, 'ro-', rf_1d(k1), psd_model(k1), 'bs-'), grid
    parms = fun_newaxes(14, 2, 5);
    xlabel('Radial Frequency (Cycles/D)')
    ylabel('1D-PSD (nm^2 m^2)')
    title(sprintf('gpi-1: 1D-PSD: LB = %i, UB = %i Cyc/D',LB,UB)) 
    legend('1D-PSD_t_r_u_e', 'Fit', 4)
    text(1.2, 2e-3, s1, 'FontSize', 14)
print -dtiff ../Figs/fig_psd1
return
%end
% -----------------------------------------------------------------------
% Calculate the 99x99 grid Zernikes and PSD
% -----------------------------------------------------------------------
%save -v6 rdn/dat_psd_parms psdmax hpf powerx

load rdn/dat_psd_parms psdmax hpf powerx LB UB D mgrid npad   
load rdn/dat_dopd_nm dopd opdi zk masktb_307
zcoef = zk;
zcoef(1:3) = 0;
opdr = pad(dopd,npad);

    % Small Circular Mask
    mgrid = npad;
    nx  = [-(mgrid-1)/2:(mgrid-1)/2];
    [x y] = meshgrid(nx);
    mask  = sqrt(x.^2 + y.^2) < mgrid/2;
    
figure(1), clf
clear opdx

dum = randn(6000,1);

for ii = 1:12
    
    faci = randn(300,1);
    kk = find(abs(faci) <= 1.5);
    facii = faci(kk);
    zki = zcoef .* (1 + facii(1:60));
    
    % Zernike components
    [opdz] = zernike_compose(mask, zki);
    opdz = pad(opdz, npad);
    aa = opdz .* masktb_307;
    rmsz = stat2d(aa);
   
    % PSD component
    rms0 = stat2d(opdr);

    % PSD Model Parameters
	psd_parms  = [D/mgrid   % mm per pixel of mask;
           	      rms0(1)   % rms WFE level (nm)
            	  hpf       % Half-power frequency (cycles/mm)
        	      powerx    % Exponent in power law
                  LB        % Controllable-band Lower-Bound (cycles/mm)
                  UB];      % Controllable-band Upper-Bound (cycles/mm)

    %[opd_abr, psdf, dr] = fun_PSDAberrate_bandpass(mask, psd_parms);
    %opd_abr = pad(opd_abr, npad);
    
    psd_mask = (rn >=LB) .* (rn <= UB);
    psdx = psd_mask .* psdmax ./ (1 + power(rn/hpf, powerx) );
    A = 1;
    H = fftshift(sqrt(psdx * A));
    wn      = rand(N, M);
    ph      = angle((fft2(wn)));
    
    fopd = ifft2(exp(i*ph) .* H);
    opd_abr = real(fopd);
   
    bb = opd_abr .* masktb_307;
    rms1 = stat2d(bb);
    fac = rms0(1) / rms1(1);  % Re-scale PSD component
    opd1 = bb * fac;  % Re-scale PSD component
   
    % Total 99x99 grid surface map
    opdttl = opd1 + aa;
    rmst = stat2d(opdttl);
    
    bb1 = opd_abr * fac;
    %rmst = stat2d(bb1);
    opdx(:,:,ii) = opdz + bb1;
    %bb1 = opdx(:,:,ii);
    bb1 = opdx(:,:,ii);
    %rmst = stat2d(bb1);

    if ii == 1, cx = rmst(1) * 3 * [-1 1]; end
   
    figure(1), subplot(3,4,ii), imagesc(opdttl,cx), axis xy image, colormap('jet'), niceaxes
        title(sprintf('Run%i: %0.1f, %0.1fnm',ii,rmst)), drawnow
end
    print -dpng ../Figs/fig_amp

    figure(2), clf, imagesc(opdttl, cx), axis xy image, colormap('jet'), niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
    print -dpng ../Figs/fig_amp2

%end    
    xv = u(:,1) * D/ dx;
    
    figure(3), clf, imagesc(xv,xv,psdx), axis xy image, colormap('jet'), %niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
    %print -dpng ../Figs/fig_amp2a
    
    figure(4), clf, loglog(rn(154,:), psdx(154,:)+eps, 'ro-', rf_1d, psd_1d, 'b'), grid
    axis([1 200 1e-4 1e2])
    
save -v6 rdn/dat_opdx_nm5 opdx
    
    
return
% ----------------------------------------------------------------------

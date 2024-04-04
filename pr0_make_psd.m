
% -------------------------------------------------------------    
%dirmap = '../hlc20161228/mcblam550_10percentlyot0.44_0_0_0.8_dmsflatwf1_10mils_13_8_pi_bandw0.1_state322_nro_';

% If you haven’t get too far,  do 20180815 (multi nickel) and 20180919 (extended pmgi ), instead of my previous suggestion.

%It just came to me that your previous evaluation of new HLC design’s impact on LOWFS is based on FALCO design, 
% traditional style mask (not of a multi nickel or extended pmgi). If I’m correct on this, could you do a quick evaluation of 
% Dwight’s new design? He has uploaded his designs in directory: 
% s383/ home/moody/jtt/xxx,  do 20180815 (multi nickel) and 20180823 (extended pmgi ) and send me your result.

%pup = fitsread('/home/esidick/Afalco/falco20200916/erkin_iris/dat/pupil_10mm_512pix.fits');
%pup = fitsread('dat_amt/pupil_cir_512pix.fits');


%load dat_redding

if 0
    
    figure(1), clf

    for ii = 1:1
        opd = trial(ii).opd_final * 1e3;
        rms0 = stat2d(opd);
        opd = pad(opd, 256) * 20 / rms0(1);
        rms0 = stat2d(opd);
        cx = rms0(1) * 3 * [-1 1];

        %subplot(2,2,ii), 
        imagesc(opd, cx);  axis image; colormap(jet), niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('hi-freq error: trial(%i).opd-final', ii))
        xlabel(sprintf('rms = %0.1f, pv = %0.1f nm', rms0))
    print -dpng ../Figs/fig_opd
    end

return
end % -------------------------------------------------------------  
%%
np = 256;
nx = [-np/2:np/2-1];
[xm,ym] = meshgrid(nx);
rm = hypot(xm, ym);
mask = rm <= np/2;

D      = 6000; % mm
mgrid  = 512;
rms0   = 10; % nm
hpf    = 10/D;  % cyc/mm
powerx = 2.5;
LB     = 0;
UB     = 500; 

	dx = D / mgrid;         % mm per pixel of mask;
    f0 = 1 / D;             % Spatial-frequency normalization factor, cyc/mm 
    du = 1 / D;             % Spatial frequency increment, cyc/mm
    A  = D * D;             % Square area of the surface map, mm^2

    % PSD Model Parameters
	psd_parms  = [D/mgrid   % mm per pixel of mask;
           	      rms0      % rms WFE level (nm)
            	  hpf       % Half-power frequency (cycles/mm)
        	      powerx    % Exponent in power law
                  LB        % Controllable-band Lower-Bound (cycles/mm)
                  UB];      % Controllable-band Upper-Bound (cycles/mm)

    % Parameters for radially-averaged (1D) PSD
    [M N] = size(mask);
    [v u] = meshgrid (1:N,1:M); 
    u = (u - (floor(N/2) + 1))/N;
    v = (v - (floor(M/2) + 1))/M;
    rf = sqrt (u.*u + v.*v) / dx;   % Radial frequency, cyc/mm

   xv = u(1,:) * D / dx; % cyc/D
    
    rn  = rf / f0;         % Normalized frequency
    [rs hs] = sort(rn(:));
        
    psd_mask = (rn >=LB) .* (rn <= UB);
    psdx = mask *1 ./ (1 + power(rf/hpf, powerx) );

    A = 1;
    H = fftshift(sqrt(psdx * A));
    wn      = rand(N, M);
    ph      = angle((fft2(wn)));
    
    fopd = ifft2(exp(i*ph) .* H);
    opdh = mask .* real(fopd);

    rmsi = stat2d(opdh);
    opdh = opdh * 10 / rmsi(1);
    rmso = stat2d(opdh);

    figure(2), clf, imagesc(opdh), axis image, colormap(jet), %niceaxes 
        parms = fun_newaxes(14, 0, 0);         
        colorbar, parms = fun_newaxes(14, 0, 0); 
        title('psd error' )
        xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rmso))
    print -dpng ../Figs/fig_psd

    opdnm = opd + opdh;

    rmsi = stat2d(opdnm);
    opdnm = opdnm * 30 / rmsi(1);
    rmso = stat2d(opdnm);

    cx = rmso(1) * 3 * [-1 1];

    figure(3), clf, imagesc(opdnm, cx), axis image, colormap(jet), %niceaxes 
        parms = fun_newaxes(14, 0, 0);         
        colorbar, parms = fun_newaxes(14, 0, 0); 
        title('total aberration: zernike + psd' )
        xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rmso))
    print -dpng ../Figs/fig_opdtot

save hex_maps/dat_opdnm opdnm opdh mask    

return
% end % -------------------------------------------------------------  
%%







% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% Erkin Sidick at JPL, Cal-tech. 9/8/2022
% -------------------------------------------------------------------------
%

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;

if 1

load dat_mp_30mm_cir98_no_psd_2to12lamD_chg4_5lam mp  

%                  dm1 = mp.dm1.V; % nm
%                  dm2 = mp.dm2.V; % nm
%                pupil = mp.P1.full.mask;
%                lyot  = mp.P4.full.mask;
%               charge = mp.F3.VortexCharge;
%       psf_resolution = mp.Fend.res; % pix per lam/D
% bandpass_wavelengths = mp.full.lambdas;  % m

lam = 500e-9;
D   = 6;
lamD = 0.1 * lam / D;
xdeg = lamD * 180 / pi;
xmas = xdeg * 60 * 60e3

%--Compute image
tic; 
fprintf('Generating the PSF for a finite stellar size…') 
mp.full.pol_conds = [0]; %--Which polarization states to use when creating an image.
mp.full.TTrms = 0; % [mas]
mp.full.Dstar = xmas; % [mas]
mp.full.Dtel  = D; % [meters]
mp.full.TipTiltNacross = 7;

%end

if 0

load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pup1_mono_512pix pupil
    mp.P1.full.mask    = pupil; %(load your file here)
 load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_mono_pm_cir98_2to12lamD_chg4_5lam_run1_cbx_itr65 dm1 dm2 
    mp.dm1.V = dm1; % nm
    mp.dm2.V = dm2; % nm
end

ImBoth = falco_get_summed_image_TipTiltPol(mp);
Im = ImBoth;
%Im = falco_get_summed_image(mp);

fprintf('done. Time = %.2f s\n',toc)

%end



    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res; % lam/D

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = Im .* maskc; 
    cb = mean(nonzeros(psfn(:)));

    fs = 14;

    figure(1), clf
        imagesc(xv,xv,log10(psfn), [-12 -8]); axis xy image, colormap(jet); 
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('starD = 0.1, [2 12], 5\\lambda-efc, 5\\lambda-scoring: Cb = %0.3g', cb)]);
        axis(16*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfa2

toc

return
end % -----------------------------------------------
%%
% wt = [wt; ilam itt mp.full.lambda_weights_all(ilam) mp.full.wsTT(itt) Npol];
load dat/dat_star_size wt Iall wt_lam wt_wsTT

    aa = wt_lam; %rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa); colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('psf positions')
        ylabel('wavelength number')
        title(sprintf('wavelength weights: %0.3f - %0.3f', min(aa(:)), max(aa(:))));
   %print -dpng ../Figs/fig_dm1

    aa = wt_wsTT; %rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(aa); colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('psf positions')
        ylabel('wavelength number')
        title(sprintf('star intesnity weights: %0.3f - %0.3f', min(aa(:)), max(aa(:))));
   %print -dpng ../Figs/fig_dm2

   bb = unique(wt_wsTT(1,:))

   figure(3), clf, plot(wt_wsTT(1,:), 'ro-'), grid
   hold on, plot(sort(wt_wsTT(1,:)), 'bs-'), hold off

   x = [-3:3];

   kk = 1:7;
   map = zeros(7,7);
   map1 = zeros(7,7);
   for jj = 1:7
       for ii = 1:7
           xx = ii + jj;
           y = abs(x(jj)) + abs(x(ii));
           map1(jj,ii) = abs(x(jj)) + abs(x(ii));
           if y == 4
               map(jj,ii) = bb(1);
           elseif y == 3
               map(jj,ii) = bb(2);
           elseif y == 2
               map(jj,ii) = bb(3);
           elseif y < 2
               map(jj,ii) = bb(4);
           end
       end
   end


    figure(4), clf
        imagesc(x,x,map); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('x-offset positions')
        ylabel('y-offset positions')
        title(sprintf('star intesnity weights: %0.3f - %0.3f', min(aa(:)), max(aa(:))));
   %print -dpng ../Figs/fig_dm2

figure(5), clf, imagesc(x,x,map1); axis xy image, colormap(jet); colorbar

return

y = wt(:,3);
x = [1:length(y)]';

    figure(12), clf, semilogy(x, y, 'ro-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('data number', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('wavelength weights', 'Fontsize',fs,'Interpreter', 'Latex');
        title('for finite sitar size simulations','Fontsize',fs,'Interpreter', 'Latex');
    print -dpng ../Figs/fig_lam

return
%end % -----------------------------------------------
%%

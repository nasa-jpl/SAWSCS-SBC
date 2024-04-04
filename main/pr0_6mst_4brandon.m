% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%

addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;

%mp.Fend.res = 4;
% ----------------- load file -----------------------

load ../dat_brandon/dat_mp mp Im ImSimOffaxis xv thput

% ------------------ calculate contrast maps ---------------

%% Take initial broadband image 
Im = falco_get_summed_image(mp);

    %--Compute the current contrast level
    cb = mean(Im(mp.Fend.corr.maskBool))
    
    [mp, thput, ImSimOffaxis] = falco_compute_thput(mp); %<***** thput ***
    
% ------------------ plot results ---------------

    dm1 = mp.dm1.V;
    dm2 = mp.dm2.V;

    pupil = mp.P1.full.mask; 
    lyot  = mp.P4.full.mask; 

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2.5) & (rm <= 12);

    psfn = Im .* maskc; 
    cb = mean(nonzeros(psfn(:)))
    
    fs = 14;

    figure(1), clf
        imagesc(xv,xv,log10(Im.*maskc)); axis xy image, colormap(jet); colorbar
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('10mm, 3-12\\lambda/D: Cb = %0.3g', cb)]);
    
    figure(2), clf
    imagesc(xv,xv,log10(ImSimOffaxis)); axis xy image, colormap(jet), colorbar
    title(sprintf('thput = %0.3f', thput)), drawnow

    figure(3), clf
    imagesc(pupil); axis xy image, colormap(jet), colorbar
    title('pupil'), drawnow

    figure(4), clf
    imagesc(lyot); axis xy image, colormap(jet), colorbar
    title('lyot mask'), drawnow


return
%end % -----------------------------------------------
%%




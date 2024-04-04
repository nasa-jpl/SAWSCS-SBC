
function [rx, yy, ynom] = fun_rms_contrast(psf0, psfx, rm, lim1, lim2)

    psfi = (psfx - psf0).^2;
    %figure(10), clf, imagesc(xv, xv, log10(psfi)), axis xy image, colormap(jet), colorbar, drawnow
    psfv = psfi(:);
    cs = psf0(:);
    rx = [];
    yy = [];
    ynom = [];
    
    dx = 0.1;

    for ri = lim1:0.1:lim2
        rx = [rx ri];
        maski = (rm >= ri-dx) .* (rm <= ri+dx);
        kk = find(maski(:));
        yy = [yy sqrt(mean(psfv(kk)))];
        
        ynom = [ynom mean(cs(kk))];
    end

return
%end % ---------------------------------------------------------------
%


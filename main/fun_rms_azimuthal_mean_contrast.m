
function [rx, yy, ynom] = fun_rms_azimuthal_mean_contrast(psf0, dlamD, lim1, lim2, flag_rms)

    if nargin < 5
        flag_rms = 0;
    end

    [nr, nc] = size(psf0);
    xv = [-nr/2:nr/2-1] * dlamD;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    %[rv, hs] = sort(rm(:));

    psfv = psf0(:);
    %psfs = psfv(hs);

    rx = [];
    yy = [];
    
    dx = 0.1;

    for ri = lim1:0.1:lim2
        
        rx = [rx ri];
        maski = (rm >= ri-dx) .* (rm <= ri+dx);
        kk = find(maski(:));
        
        if flag_rms
            yy = [yy sqrt(mean(psfv(kk).^2))];
        else
            yy = [yy mean(psfv(kk))];
        end
        
    end

return
%end % ---------------------------------------------------------------
%


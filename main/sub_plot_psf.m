% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.
    
[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, 24);
    [xc2, yc2] = fun_circle(xv, 24);

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)))

    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn(:)));
    
    fs = 14;

    tc = fom(end)*10;

    figure(10), clf
        imagesc(xv,xv,log10(psfn), [-12 -8]); axis xy image, colormap(jet); 
        %imagesc(xv,xv,log10(aa), [-12 -8]); axis xy image, colormap(jet); 
        %hold on, plot(xc, yc, 'w', xc1, yc1, 'w', xc2, yc2, 'w', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('12nm: 5\\lambda, Csb = [%0.3g %0.3g], tc = %0.1f', cb1, cb, tc) '%']);
        axis(13*[-1 1 -1 1])
   figure(10), print -dpng ../Figs/fig_psf1

return
%end % -----------------------------------------------
%%


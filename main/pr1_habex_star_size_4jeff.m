% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
%
% Erkin Sidick at JPL, Cal-tech. 9/19/2022
% -------------------------------------------------------------------------
%

%%
% wt = [wt; ilam itt mp.full.lambda_weights_all(ilam) mp.full.wsTT(itt) Npol];
load dat/dat_star_size wt Iall wt_lam wt_wsTT
load dat/dat_star_size_xoffset_yoffset xoffset yoffset

    figure(1), clf
        imagesc(wt_lam); colormap(jet); %niceaxes
        colorbar, 
        xlabel('psf positions')
        ylabel('wavelength number')
        title(sprintf('wavelength weights: %0.3f - %0.3f', min(wt_lam(:)), max(wt_lam(:))));

   x = [-3:3];

   bb = unique(wt_wsTT(1,:));

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

    figure(2), clf
        imagesc(x,x,map); axis xy image, colormap(jet); 
        colorbar,   
        xlabel('x-offset positions')
        ylabel('y-offset positions')
        title(sprintf('star intesnity weights: %0.3f - %0.3f', min(map(:)), max(map(:))));

     figure(3), clf, plot(xoffset, yoffset, 'ro'), axis square, grid
        parms = fun_newaxes(16, 0, 10);  
        xlabel('X-Axis Field-Angle [\lambda/D');
        ylabel('Y-Axis Field-Angle [\lambda/D]');
        title(sprintf('Tip-Tilt Offset Points Used: %i', length(xoffset)) );

return
%end % -----------------------------------------------
%%

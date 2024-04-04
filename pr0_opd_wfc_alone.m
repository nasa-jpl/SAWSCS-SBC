    load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_bb_30nm_bw20_355nm_flat opd_target amp_target ce_target
    load ~esidick/Afalco/falco20200916/macos/IFopd/distorted_opd_amp_bb_30nm_bw20_355nm_flat opd_nm_distorted amp_distorted ce_distorted

if 0    
    aa = opd_nm_distorted; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    opdi = aa;
    figure(1), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('pre-wfc opd: rms = %0.1fnm', rms0(1)));
   print -dpng ../Figs/fig_opdi

    aa = amp_distorted; aa = aa / max(aa(:)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    ampi = aa;
    figure(2), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('pre-wfc amp: rms = %0.3g', rms0(1)));
   print -dpng ../Figs/fig_ampi

   disp('paused'), pause

   end

   flag = 1

    ce = ce_distorted ./ (ce_target + eps);
    lambda = 355e-9;

    np   = 512/2;
    maskopd = abs(pad(mp.P1.full.mask,np)) > 0;

    %opd = maskopd .* (1e9 * atan2(imag(ce),real(ce)) * lambda / 2 / pi);
    opd = maskopd .* (1e9 * angle(ce) * lambda / 2 / pi);
    

    aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    %opdt = aa; rms0 = stat2d(opdi-opdt);
    figure(1), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        %imagesc(opd); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('\\Deltaopd: rms = %0.1fnm', rms0(1)));
   print -dpng ../Figs/fig_opdt

    aa = amp_target; aa = aa / max(aa(:)); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    ampt = aa; rms0 = stat2d(ampi-ampt);
    figure(2), clf
        %imagesc(aa); axis xy image, colormap(jet); niceaxes
        imagesc(ampi-ampt); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('\\Deltaamp: rms = %0.3g', rms0(1)));
   print -dpng ../Figs/fig_ampt
   
   return

if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np   = 512/2;
indx = mp.wfc.indx;
nwfc = 10;

%flag = 5
%disp('dm1/dm2-wfc paused'), pause

ma = mp.dm1.Nact;
ma2 = ma * ma;
du0 = zeros(ma2,1);

maskopd = abs(pad(mp.P1.full.mask,np));

ndm1 = length(mp.wfc.km1);
ndm2 = length(mp.wfc.km2);
ndm  = ndm1 + ndm2;

t1 = clock;
tsum = 0;

clear dm1x dm2x opdx ampx
rmsx = [];

dm1 = mp.dm1.V;
dm2 = mp.dm2.V;

% ---------------------------------------

for jj = 1:nwfc+1

    ce = pad(EP3,np);
    pha_nm = 1e9 * atan2(imag(ce),real(ce)) * lambda / 2 / pi;
    %pha_nm = real(pha_nm);
    opd = (pha_nm  - mp.wfc.opd_target) .* maskopd; %mp.wfc.maskopd;
    amp = (abs(ce) - mp.wfc.amp_target) .* maskopd;

    if 0
        m = 256; 
        bb = opd(:);
        k1 = find(bb>0); ap = bb(k1); meanp = mean(ap)*100;
        k2 = find(bb<0); am = bb(k2); meanm = mean(am)*100;
        k1 = find(bb>meanp); bb(k1) = 0; 
        k2 = find(bb<meanm); bb(k2) = 0; 
        opdf = reshape(bb,m,m);
        maski = abs(opdf) > 0;

        opd = opd .* maski;
        amp = amp .* maski;
    end

    %opd1= (pha_nm  - mp.wfc.opd_target*0) .* maskopd; %mp.wfc.maskopd;
    
    %disp('paused'), pause(30)

    dw = m2v(opd, indx);
    da = m2v(amp, indx);
    dw = [dw(:); da(:)];
    ddu =  -mp.wfc.G * dw;

    du1 = du0;
    du1(mp.wfc.km1) = ddu(1:ndm1);

    du2 = du0;
    du2(mp.wfc.km2) = ddu(ndm1+1:ndm);

    del_dm1 = reshape(du1, ma, ma);
    del_dm1 = del_dm1';

    del_dm2 = reshape(du2, ma, ma);
    del_dm2 = del_dm2';

    dm1 = dm1 + del_dm1;
    mp.dm1.V = dm1;

    dm2 = dm2 + del_dm2;
    mp.dm2.V = dm2;

    opdx(:,:,jj) = opd;
    ampx(:,:,jj) = amp;
    dm1x(:,:,jj) = dm1;
    dm2x(:,:,jj) = dm2;

    rmsi = stat2d(opd)

    rmsx = [rmsx; rmsi];

    %figure(10), clf, plot(rmsx(:,1), 'ro-'), grid, drawnow

    delt = etime(clock, t1);
    t1   = clock;
    tsum = tsum + delt;
    count = [jj nwfc+1 round(delt) round(tsum)]

end % nwfc - loop

    %save('~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dmx_opdx_08nm dmx dm1x dm2x opdx opd_in
    %save(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx')

end %--END OF FUNCTION
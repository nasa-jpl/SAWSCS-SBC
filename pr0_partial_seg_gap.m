% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;
addpath(genpath('/proj/exospin/users/joegreen/'));

pp = 1


pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');

if 0

ps = size(pupil);

n = 256;
N = size(pupil,1);
nh = (n-1)/2;
nx = [-nh:nh];
nx = [-N/2:N/2-1];
[xm,ym] = meshgrid(nx+1, nx);
rm = hypot(xm,ym);
maskb = rm <= 252/2;
%maskb = pad(mask,N);
%maskb = maskb';

figure(1), clf, imagesc(pupil+maskb);  axis xy image; colormap(jet), %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    %title('opd at apod: falco, no dm'), 
    %title('\Deltaopd at apod: falco, dm2-act(32,32) = 5nm'), 
    %xlabel(sprintf('rms = %0.2f, pv = %0.1f nm', rmsa)), 
%print -dpng ../Figs/fig_opd1

return
%end % -----------------------------------------------
%%
load /home/sab/Projects/6MST/macos/data/MonoOPD.mat OPDLyot OPDnom

fac = sqrt((16/24))

masklyo = double(OPDLyot ~= 0);
maskpm  = double(OPDnom ~= 0);
opdnom_nm = OPDnom * 1e6;

aa = maskpm .* pupil;
bb = ((1 - aa) * fac) .* maskpm;
mask0 = aa + bb;

%mask0 = maskpm .* pupil; % for gap = 24mm only

%fitswrite(masklyo, 'maps_20230530/lyo.fits')
%fitswrite(mask0,  'maps_20230530/pupil_24mm.fits')
%save maps_20230530/opdnom_nm_00mm opdnom_nm

kc = 128 + [-5:5];
kr = 204 + [-5:5];

figure(1), clf, imagesc(mask0, [0 1]);  axis xy image; colormap(jet), %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    %title('opd at apod: falco, no dm'), 
    %title('pupil + lyot: gap = 00mm'), 
    %title('white = pupil, gray = lyot'), 
    title('pupil map: gap = 8mm'), 
    %xlabel(sprintf('rms = %0.2f, pv = %0.1f nm', rmsa)), 
%print -dpng ../Figs/fig_pup2

figure(2), clf, imagesc(mask0(kr,kc), [0 1]);  axis xy image; colormap(jet), %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    %title('opd at apod: falco, no dm'), 
    %title('pupil + lyot: gap = 00mm'), 
    title(sprintf('gap = 8mm: gap-int = %0.3f', fac^2)), 
    %xlabel(sprintf('rms = %0.2f, pv = %0.1f nm', rmsa)), 
%print -dpng ../Figs/fig_pup2a

return
    aa = OPDnom*1e6; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('Nom-OPD: rms = %0.2f, pv = %0.1fnm', rms0));
   %print -dpng ../Figs/fig_dm1
return
%end % -----------------------------------------------
%%
%figure(1), clf, imagesc(pupil); axis xy image, colormap(gray); %niceaxes
%figure(2), clf, imagesc(maskpm); axis xy image, colormap(gray); %niceaxes
%return
nbin = 32;
nbin1 = nbin + 2;
np = 256;
np1 = 260;
pup = pupil .* maskpm;
maskbig = pad(maskpm, np1);
pup1 = pad(pup, np1);
maskx = ones(np1,np1) - maskbig + pup1;
%figure(1), clf, imagesc(maskx); axis xy image, colormap(gray); %return %niceaxes
%figure(2), clf, imagesc(maskpm+pupil); axis xy image, colormap(jet); return %niceaxes

pupbig = binup(maskx, nbin);  % <<*****************

ns = size(pupbig,1)

figure(1), clf, imagesc(pupbig); axis xy image, colormap(jet); %niceaxes

%end

nstep = nbin * 2;
k0 = [-nstep/2:nstep/2+1];

masko = maskx;

for jj = nbin+1:ns-nbin-1
    kr = jj + k0;
    for ii = nbin+1:ns-nbin-1
        kc = ii + k0;
        mapi = pupbig(kr,kc);
        masko(jj,ii) = mean(mapi(:));
    end
    if mod(jj, 100) == 1, count = [jj ns-nbin], end
end
%end

masks = pad(bindown(pad(masko,np*nbin), nbin),np) .* maskpm;

figure(1), clf, imagesc(masko); axis xy image, colormap(jet); colorbar, %niceaxes
figure(2), clf, imagesc(masks); axis xy image, colormap(jet); colorbar, %niceaxes
        

%return
%end % -----------------------------------------------
%%
nx = [-np/2+0.5:np/2-0.5];
[xm,ym] = meshgrid(nx);
rm = hypot(xm,ym);
masknom = rm <= 252/2;
idd = 2
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm psdnm opdnm_only  
%load dat_dm_err1/opd30_erkin dopd_c_nm
maske = opdnm_only ~= 0;
%maskg = dopd_c_nm ~= 0;
[m1 m2 n1 n2] = find_boundary(maske,0.5); wid_mono = [m2-m1 n2-n1] + 1
[m1 m2 n1 n2] = find_boundary(pupil,0.5); wid_seg = [m2-m1 n2-n1] + 1

%end

figure(3), clf, imagesc(pupil+masknom); axis xy image, colormap(jet); colorbar, %niceaxes
figure(4), clf, imagesc(pupil.*masknom); axis xy image, colormap(jet); colorbar, %niceaxes

%end

%masknom = masknom / max(masknom(:));
pup0 = pupil .* masknom;
Itot = sum(masknom(:).^2);
Ipup = sum(pup0(:).^2) / Itot;
Igap = 1 - Ipup;

masks = pad(bindown(pad(masko,np*nbin), nbin),np) .* masknom;
masks = masks / max(masks(:));
Ix = 1 - sum(masks(:).^2) / Itot;
x0 = (Igap/Ix) * 24

%fitswrite(masks, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_15mm.fits');

%figure(2), clf, imagesc(masknom); axis xy image, colormap(jet); %colorbar, %niceaxes
figure(3), clf, imagesc(masks); axis xy image, colormap(jet); %colorbar, %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('no gap: total energy = 1'));
        %title(sprintf('with gap = 24mm: energy = %0.3f', Ipup));
        title(sprintf('Gap-2: energy = %0.3f, gap = %0.3fmm', Ix, x0));
   print -dpng ../Figs/fig_gap1

return
%end % -----------------------------------------------
%%
masknom = double(masknom / max(masknom(:)));
masks = pad(bindown(pad(masko,np*nbin), nbin),np) .* masknom;
masks = double(masks / max(masks(:)));

%fitswrite(masks, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');

maskgap = double((masks > 0) .* (masks < 1));
maskdel = masks .* maskgap;

%end

idd = 1

masks_ones = pup0 == 1;
masks_gray = (masknom - masks_ones) .* masknom;

Inom = sum(masks_ones(:).^2);
Ig24 = sum(masks_gray(:).^2) / Inom / 24;

end
masks1 = pup0;

masks_ones = masks == 1;
masks_gray = (masks - masks_ones); % * sqrt(8/81.7243);
Ig = sum(masks_gray(:).^2) / Inom;

maskg = masks_gray > 0;
[y, indx] = m2v(maskg);
y = m2v(masks_gray, indx);
[ys, hs] = sort(y);
x = [1:length(y)];

%figure(1), clf, plot(x, ys, 'ro'), grid, return

x0 = Ig / Ig24

pup_new = masks_ones + masks_gray;

%fitswrite(masks, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');

figure(1), clf, imagesc(pup0); axis xy image, colormap(jet); colorbar, title('pup0')%niceaxes
figure(2), clf, imagesc(masks); axis xy image, colormap(jet); colorbar, title('masks')%niceaxes
figure(4), clf, imagesc(masks_gray); axis xy image, colormap(jet); colorbar, title('masks_gray')%niceaxes

figure(3), clf, imagesc(masks+pup0); axis xy image, colormap(jet); %colorbar, %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('no gap: total energy = 1'));
        %title(sprintf('with gap = 24mm: energy = %0.3f', Ipup));
        title(sprintf('pupil: gap = %imm', round(x0)));
   %print -dpng ../Figs/fig_gap1


return
%end % -----------------------------------------------
%%


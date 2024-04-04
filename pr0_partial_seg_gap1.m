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
np = 256;
nx = [-np/2+0.5:np/2-0.5];
[xm,ym] = meshgrid(nx);
rm = hypot(xm,ym);
masknom = rm <= 252/2;
pup0 = pupil .* masknom;

%figure(1), clf, imagesc(pupil); axis xy image, colormap(gray); %niceaxes
%figure(2), clf, imagesc(maskpm); axis xy image, colormap(gray); %niceaxes
%return

nbin = 12;

nbin1 = nbin + 2;
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
masks = masks / max(masks(:));

%end

masks1 = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');
%ps = size(masks1), return

%end
x = [-np/2:np/2-1];

figure(1), clf, imagesc(x,x,masks1); axis xy image, colormap(jet); colorbar, %niceaxes
figure(2), clf, imagesc(x,x,masks1); axis xy image, colormap(jet); colorbar, %niceaxes
hold on, contour(x,x, pup0, 1, 'w'), hold off
        

return
%end % -----------------------------------------------
%%
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

%end

masks = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');


%fitswrite(masks, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');

maskgap = double((masks > 0) .* (masks < 1));
maskdel = masks .* maskgap;

%end

idd = 1

masks_ones = pup0 == 1;
masks_gray = (masknom - masks_ones) .* masknom;

Inom = sum(masks_ones(:).^2);
Ig24 = sum(masks_gray(:).^2) / Inom / 24;

%end
masks1 = pup0;

masks_ones = masks == 1;
masks_gray = (masks - masks_ones); % * sqrt(8/81.7243);
Ig = sum(masks_gray(:).^2) / Inom;

maskg = masks_gray > 0;
[y, indx] = m2v(maskg);
y = m2v(masks_gray, indx);
[ys, hs] = sort(y);
x = [1:length(y)];

figure(1), clf, plot(x, ys, 'ro'), grid, %return
    parms = fun_newaxes(14, 0, 5); 
    xlabel('vector index')
    ylabel('sorted gray pix values')
%print -dpng ../Figs/fig_gap

x0 = Ig / Ig24

pup_new = masks_ones + masks_gray;

%fitswrite(masks, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/masknom_nbin32.fits');

%end

x = [-np/2:np/2-1];

figure(1), clf, imagesc(pup0); axis xy image, colormap(jet); colorbar, title('pup0')%niceaxes
figure(2), clf, imagesc(masks); axis xy image, colormap(jet); colorbar, title('masks')%niceaxes
figure(4), clf, imagesc(masks_gray); axis xy image, colormap(jet); colorbar, title('masks_gray')%niceaxes

figure(3), clf, imagesc(x,x,masks); axis xy image, colormap(jet); %colorbar, %niceaxes
        %hold on, contour(x,x,pup0, 1, 'b'), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('no gap: total energy = 1'));
        title(sprintf('equivalent gap = %0.3fmm', x0));
        %title(sprintf('gray-scale gap: blue = orig gap'));
   print -dpng ../Figs/fig_gap1

   tx = 'print -dpng ../Figs/fig_gap1'

return
end % -----------------------------------------------
%%
Colx = 'rbgcmkrbgcmk';
x = [1:length(y)];
a = [2000:500:6000];
xfac = [15.95 18:2:28 29.5];
%xfac = 15.95;
%xfac = 20;
xfac = 29.3
ns = length(xfac);

a0 = 10000;

clear legx ymat
legx{1} = 'data';

figure(1), clf, plot(x, ys, 'ro'), grid, %return

for ii = 1:ns
    x2 = x*xfac(ii);
    %y = x2 ./ (x2 + a(ii));
    y = (10/xfac(ii))*x2 ./ (x2 + a0);
    ymat(ii,:) = y;
    hold on, plot(x, y, Colx(ii)), hold off
    %legx{ii+1} = sprintf('a = %i', a(ii));
    legx{ii+1} = sprintf('xfac = %0.2f', xfac(ii));
end
    parms = fun_newaxes(14, 2, 5); 
    xlabel('vector index')
    ylabel('sorted gray pix values')
    title(sprintf('data & multipliers: a = %i',a0))
    legend(legx, 'Location', 'Southeast')
%print -dpng ../Figs/fig_gap

%return
%end % -----------------------------------------------
%%
y = m2v(masks_gray, indx);
[ys, hs] = sort(y);
ns = length(y);
x = [1:ns];

gx = [];

idd = 1

%for ii = 1:length(a)
for ii = 1:length(xfac)

    y2 = y*0;
    fac = ymat(ii,:);
    y2(hs) = ys .* fac(:);
    masks_grayi = y2; % * sqrt(8/81.7243);
    Ig = sum(masks_grayi.^2) / Inom;
    x0 = Ig / Ig24;
    gx = [gx x0];
end
    
%y1 = y(hs);
%figure(2), clf, plot(x, y, 'ro', x, y2, 'bs'), grid
figure(2), clf, plot(xfac, gx, 'ro-'), grid
    parms = fun_newaxes(14, 2, 5); 
    xlabel('parameter xfac')
    ylabel('equivalent gap width [mm]')
    title(sprintf('equivalent gap: a = %i',a0))
%print -dpng ../Figs/fig_gap1

gapi = v2m(masks_grayi, indx);
pup  = masks_ones + gapi;

%pup = double(masknom);
%x0 = 0;

fitswrite(pup, '/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');

xs = [-np/2:np/2-1];

figure(3), clf, imagesc(xs,xs,pup); axis xy image, colormap(jet); %colorbar, %niceaxes
        %hold on, contour(x,x,pup0, 1, 'b'), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('no gap: total energy = 1'));
        title(sprintf('equivalent gap = %0.1fmm', x0));
        %title(sprintf('gray-scale gap: blue = orig gap'));
   print -dpng ../Figs/fig_gap

   tx = 'print -dpng ../Figs/fig_gap'

return
%end % -----------------------------------------------
%%

%% EXAMPLE MATLAB RANDOM GENERATION

   addpath ~esidick/matlab_erkin/
   addpath /proj/exospin/users/joegreen/matlab/

   mp.Fend.res = 4;
   fs = 16;

if 1

    load /proj/jwst2/6MST/data/6MST_wfs_seg512.mat
    indx_pup = indx;
    opdnom = v2m(wnom, indx_pup);
    pup = opdnom ~= 0;

    %load /proj/jwst2/6MST/data/6MST_lyot512.mat
    load /proj/jwst2/6MST/data/6MST_lyot512.mat
    indx_lyo = indx;
end

if 1

%whos, return

    ns = size(pup,1);
    nx = [-ns/2+1:ns/2];

%pup = lyot;
[m1 m2 n1 n2] = find_boundary(pup>0, 0.5);
wid1 = m2-m1+1;

lyo = v2m(lyot, indx_lyo); %mp.P4.full.croppedMask;
lyo = (lyo ~= 0);
[m1 m2 n1 n2] = find_boundary(lyo>0, 0.5);
wid2 = m2-m1+1;

figure(1), clf, imagesc(pup*2-lyo), axis xy image, colormap(gray), niceaxes
    parmx = fun_newaxes(fs, 0, 0);  
    colorbar, parmx = fun_newaxes(fs, 0, 0);  
    %title(sprintf('pupil: dia = %i', wid1))
    title(sprintf('pupil + lyot'))
%print -dtiff ../Figs/fig_pup1

figure(2), clf, imagesc(nx,nx, lyo), axis xy image, colormap(gray), %niceaxes
    parmx = fun_newaxes(fs, 0, 0);  
    colorbar, parmx = fun_newaxes(fs, 0, 0);  
    title(sprintf('lyo: dia = %ipix = %0.3fD', wid2, wid2/wid1))
%print -dtiff ../Figs/fig_lyo

aa = opdnom*1e6; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
figure(3), clf, imagesc(aa), axis xy image, colormap(jet), %niceaxes
    parmx = fun_newaxes(fs, 0, 0);  
    colorbar, parmx = fun_newaxes(fs, 0, 0);  
    title(sprintf('nominal opd: rms = %0.1f, pv = %0.1fnm', rms0))
%print -dtiff ../Figs/fig_opd

%fitswrite(single(pup), 'maps/pup512_ogap3_igap1')
fitswrite(single(lyo), 'maps/lyo512_ogap3_igap1_v2.fits')
%fitswrite(single(aa),  'maps/opd512_ogap3_igap1_nom_nm')

return
end % --------------------------------------------------------
%

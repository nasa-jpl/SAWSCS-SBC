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

% It includes a pupil mask. 6DOF RB sensitivities for 19 PM segments and m2
%load /proj/jwst2/6MST/data/dwdx_6mst_19HexSeg_m2_urad_um_12222023 OPDMask_g dwdx indx wnom
%whos, return

pup = pad(OPDMask_g, 512);
fitswrite(single(pup), 'hex_maps/pupil.fits')

figure(1), clf, imagesc(pup);  axis image; colormap(jet), %niceaxes
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title('19-segment pm'), 
    %xlabel(sprintf('rms = %0.3g, pv = %0.3g', stat2d(amp))), 
print -dpng ../Figs/fig_pup

return
%end % -----------------------------------------------
%%



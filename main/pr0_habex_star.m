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

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%tval = tx(end)
%cb0 = cbx(end)
%figure(1), clf, semilogy(cbx,'ro-'), grid, %return

mp.Fend.res = 4;

tic

if 0


%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

% lamx = mp.sbp_centers in [m]


%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

% mp.Fend.res = 4; %3; %--Sampling [ pixels per lambda0/D]


%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design_erkin3_habex

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;


if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
mp.compact.flagGenPupil = false;

%pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_cirle512.fits');
%pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_atsa512.fits');
%apod  = fitsread('/home/esidick/Atsa/Jeff/map/apod_atsa512.fits');
%pupil = fitsread('/home/esidick/Atsa/Jeff1/Data/map_pupil.fits');
%apod  = fitsread('/home/esidick/Atsa/Jeff1/Data/map_apod.fits');

%pupil = fitsread('/home/esidick/2022habex/dat/pupil_15mm.fits');
%apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_sharp_v2.fits');
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_soft.fits');
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_00mm_512pix_soft.fits');
load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pupil_512pix_gap30mm pupil
apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
lyo   = 'inside setup/falco_flesh_out_workspace';

amp_save = pupil;

% save dat/dat_wfc_residual_psd del_psd opdz1 psd_nm indx_opd dm1

%load ~esidick/2022habex/psd/dat_9psd_opd_m opdm opdx  
%load ~esidick/2022habex/dat/dat_wfc_residual_psd psd_nm dm1  
load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
%load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo1000 del_psd psd_nm indx_opd dm1

%load ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo1000 del_psd psd_nm indx_opd dm1

dm1_save = dm1 * 1e9;
%figure(1), clf, imagesc(del_psd), axis xy image, colormap(jet), colorbar
%figure(2), clf, imagesc(dm1_save), axis xy image, colormap(jet), colorbar
%return
%rms0 = stat2d(del_psd), disp('paused'), pause

%pupil = pupil .* exp(1i*2*pi*del_psd*1e-9/mp.lambda0);

mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

mp.full.flagGenApod = false;
mp.compact.flagGenApod = false;
mp.P3.compact.mask     = apod * 0 + 1;
mp.P3.full.mask        = apod * 0 + 1;

% for apodizer, 
mp.flagApod = false;    %--Whether to use an apodizer or not
%mp.flagApod = true;    %--Whether to use an apodizer or not

end

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%% Step 5: Perform the Wavefront Sensing and Control
% falco000/setup

    [mp, out] = falco_flesh_out_workspace(mp);

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization

mp.Nitr = 10; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use

%end
%% -------------------------------------------------------------------
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg4_5lam_run1_cbx_itr140 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/z1z2z3_30mm_cir98_no_psd_2to12lamD_chg4_4lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/z1z2z3_30mm_cir98_no_psd_2to12lamD_chg6_4lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value


fs = 14;

mx = [1:length(cbx)] - 1;
s1 = sprintf('last Cb = %0.3g', cbx(end));
    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

if 0
    load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
    %load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo1000 del_psd psd_nm indx_opd dm1
    %load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pup1_mono_512pix pupil
    load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pupil_512pix_gap30mm pupil

    %pupil = fitsread('/home/esidick/2022habex/dat/pupil_00mm_512pix_soft.fits');
    %pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_soft.fits');
    %pupil = pupil .* exp(1i*2*pi*del_psd*1e-9/mp.lambda0);

    n = size(pupil,1);
    nx  = [-n/2:n/2-1] / (n/2);
    [xxm, yym] = meshgrid(nx);
    xoff = 0;
    yoff = 0;
    pha = (xoff .* xxm + yoff .* yym);    
    %pupil = pupil .* exp(1i*pi*pha);
    
    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)
    mp.F3.VortexCharge = 6; %4; %--Charge of the vortex mask

end

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    lim = 12;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

xx = [0.01:0.01:0.1];
ns = length(xx);
mx = [0];
cbxx = [0.248]*1e-10; % chrg4

xx = [0.01 0.05 0.1];
ns = length(xx);
mx = [];
cbxx = [];

%cbxx = [0.1208 3.6332]*1e-10;  %chrg6

clear psfm

for ii = 1:ns
    xi = xx(ii);

if 1

lam = 500e-9;
D   = 6;
lamD = xi * lam / D;
xdeg = lamD * 180 / pi;
xmas = xdeg * 60 * 60e3

%--Compute image
tic; 
fprintf('Generating the PSF for a finite stellar size…') 
mp.full.pol_conds = [0]; %--Which polarization states to use when creating an image.
mp.full.TTrms = 0; % [mas]
mp.full.Dstar = xmas; % [mas]
mp.full.Dtel  = D; % [meters]
mp.full.TipTiltNacross = 7;

Im = falco_get_summed_image_TipTiltPol(mp);

    psfm(:,:,ii) = Im;

end
%end

%mp.thput_eval_x = 0.01; %6; % x location [lambda_c/D] in dark hole at which to evaluate throughput
%mp.thput_eval_y = 0; % y location [lambda_c/D] in dark hole at which to evaluate throughput
%Im = falco_get_summed_image(mp);

%[mp, thput, ImSimOffaxis] = falco_compute_thput(mp); %<***** thput ***

%end

%Imnom = Im;
%Im = ImSimOffaxis;
%Im = Im * max(Imnom(:)) / max(Im(:));

    psfn = Im .* maskc; 
    cb = mean(nonzeros(psfn(:)))

    mx = [mx xi];
    cbxx = [cbxx cb];

    figure(1), clf, semilogy(mx, cbxx, 'ro-'), grid, drawnow

    count = [ii ns]

end
title('charge 4')
print -dpng ../Figs/fig_star_tiptilt_chrg6

    %end

%load ../dat_paper/dat_fig27_10lamD_chrg6_Im_v2 Im cb maskc xv
%save ../dat_paper/dat_fig27_05lamD_chrg6_Im Im cb maskc xv
%save ../dat_paper/dat_fig27_10lamD_chrg6_Im Im cb maskc xv
%save ../dat_paper/chrg4_cbxx_vs_lamD mx cbxx
load ../dat_paper/z1z2z3_chrg4_cbxx_vs_lamD mx cbxx
save ../dat_paper/z1z2z3_chrg4_cbxx_vs_lamD mx cbxx psfm
%save ../dat_paper/z1z2z3_chrg6_cbxx_vs_lamD mx cbxx psfm
%cb
end

Im = psfm(:,:,2);
    
    fs = 14;

    figure(2), clf
        %imagesc(xv,xv,log10(Im)); axis xy image, colormap(jet); %niceaxes
        imagesc(xv,xv,log10(Im.*maskc), [-10 -6]); axis xy image, colormap(jet); niceaxes
        %imagesc(xv,xv,log10(aa), [-12 -8]); axis xy image, colormap(jet); 
        %hold on, plot( xc1, yc1, 'b',  'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Charge 4', 'FontSize', 20);
        %title([sprintf('10mm, psd, chrg4, bit = 16, 2-%i\\lambda/D, v2-r1-itr225: Cb = %0.3g', lim, cb)]);
        %title([sprintf('5\\lambda, \\lambda/100, chrg4, 2-%i, itr70: Cb = %0.3g', lim, cb)]);
        %title([sprintf('5\\lambda-control, 7\\lambda-scoring: Cb = %0.3g', cb)]);
        %title([sprintf('Star-D = 0.01\\lambda/D: Cb = %0.3g', cb1)]);
        title([sprintf('0.05\\lambda/D')], 'FontSize', 20);
        %title([sprintf('xoff = 3\\lambda/D: Cb = %0.3g', cb1)]);
        %xline(-3, 'w', 'Linewidth',2);
        %yline(0, 'w','Linewidth',2);
        %axis(13*[-1 1 -1 1])
print -dpng ../Figs_paper/fig31_chrg4_05lamD
%print -dpng ../Figs/fig_psf3

return

    figure(1), clf, semilogy(mx, cbx, 'ro-'), grid, 
    %figure(12), clf, semilogy(mx1, cbx(mx1+1), 'ro-', mx2, cbx(mx2+1), 'bs-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('30mm, 0.98D, $\Delta$psd-$\lambda$/100, 5$\lambda$-efc, chrg4, [2 12] $\lambda$/D, ','Fontsize',fs,'Interpreter', 'Latex');
        legend(s1, 'Location', 'Northeast')
    %print -dpng ../Figs/fig_efc1


return
%end % -----------------------------------------------
%%


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

%EXAMPLE_defaults_LUVOIRB_VC_design_erkin1
EXAMPLE_defaults_LUVOIRB_VC_design_erkin3_habex

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%%--[OPTIONAL] Start from a previous FALCO trial's DM settings
% fn_prev = 'ws_Series0002_Trial0001_HLC_WFIRST20180103_2DM48_z1_IWA2.7_OWA10_6lams575nm_BW12.5_EFC_30its.mat';
% temp = load(fn_prev,'out');
% mp.dm1.V = temp.out.DM1V;
% mp.dm2.V = temp.out.DM2V;
% clear temp

% %--DEBUGGING
% mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
% mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
% mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
%%

if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
mp.compact.flagGenPupil = false;

%pupil = fitsread('/home/esidick/2022habex/dat/pupil_15mm.fits');
%apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_sharp_v2.fits');
pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_soft.fits');
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_00mm_512pix_soft.fits');
apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
lyo   = 'inside falco_flesh_out_workspace';

load ~esidick/2022habex/psd/dat_9psd_opd_m opdm opdx  
%figure(1), clf, imagesc(opdx), axis xy image, colormap(jet), colorbar
%return

pupil = pupil .* exp(1i*4*pi*opdx/mp.lambda0);

mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

mp.full.flagGenApod = false;
mp.compact.flagGenApod = false;
mp.P3.compact.mask = apod * 0 + 1;
mp.P3.full.mask = apod * 0 + 1;

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

    mp.controller = 'plannedefc';

       SetA = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
            %repmat([1, 1j, 12, 1, 1], [5, 1]);

       mp.ctrl.sched_mat = [...
           repmat(SetA,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetC2,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

a = 'figure(1), print -dpng ../Figs/fig_psfx'

%end

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_1to12lamD_chg4_7lam_run5_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg4_7lam_run5_cbx_itr170 dm1 dm2 tx cbx fom Im beta_value
    
ms = length(cbx);

id = 0;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;

    figure(11), clf, semilogy(cbx, 'ro-'), grid, %pause
    if id == 0, cbx = []; end %cbx([ms-1 ms]);

%[mp, out, cbx] = falco_wfsc_loop_erkin1(mp, out, cbx, dm1a, dm2a);

%end % ------------------------------------------
%%

Im = falco_get_summed_image(mp);
[mp, thput, ImSimOffaxis] = falco_compute_thput(mp); %<***** thput ***


[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    
    %psfn = Im .* maskc;
    %cb = mean(nonzeros(psfn(:)))

    %--Compute the current contrast level
    cb = mean(Im(mp.Fend.corr.maskBool))

end

    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1] = fun_circle(xv, 12);
    
    figure(1), clf
    imagesc(xv,xv,log10(Im.*maskc+eps), [-6 -2]); axis xy image, colormap(jet), 
    %hold on, plot(xc, yc, 'b', xc1, yc1, 'b', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('pre-efc: Cb = %0.2g', cb)), drawnow
    print -dpng ../Figs/fig_psf

    %close all

    figure(2), clf
    imagesc(xv,xv,log10(ImSimOffaxis), [-5 -1]); axis xy image, colormap(jet), colorbar
    %hold on, plot(xc, yc, 'y', xc1, yc1, 'r', 'Linewidth', 2), hold off
    title(sprintf('falco: with psd, pre-efc, thput = %0.3f', thput)), drawnow
    print -dpng ../Figs/fig_thput


return
%end % -----------------------------------------------
%

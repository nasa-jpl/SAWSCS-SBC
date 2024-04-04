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

%{
load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg6_7lam_run2_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value

    bit = 10;
    delz = 2000/(2^bit-1), % nm

    dm1a = dm1;
    dm2a = dm2;

    ni = round(dm1a / delz);
    dm1 = ni * delz;

    ni = round(dm2a / delz);
    dm2 = ni * delz;
figure(1), clf, imagesc([dm1a dm1 dm1a-dm1]), axis image, colormap(jet), colorbar, 
figure(2), clf, imagesc([dm1a-dm1]), axis image, colormap(jet), colorbar, return
%}

mp.Fend.res = 4;

tic

if 1


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
EXAMPLE_defaults_LUVOIRB_VC_design_erkin2_habex

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

%pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_cirle512.fits');
%pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_atsa512.fits');
%apod  = fitsread('/home/esidick/Atsa/Jeff/map/apod_atsa512.fits');
%pupil = fitsread('/home/esidick/Atsa/Jeff1/Data/map_pupil.fits');
%apod  = fitsread('/home/esidick/Atsa/Jeff1/Data/map_apod.fits');

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

pupil = pupil .* exp(1i*2*2*pi*opdx/mp.lambda0);

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

%end % -------- end of main ---------------------------------

%save ../dat_02jun2022/dat_main_mp mp

if 0

load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_25lamD_7lam_run5_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    mp.dm1.V = dm1*1;
    mp.dm2.V = dm2*1;
    
    cb0 = cbx(end)
    
    ms = length(cbx);

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;
    
    figure(1), clf, semilogy(cbx, 'ro-'), grid, 
    %cbx = []; %cbx([ms-1 ms]);

    sub_plot_psf

    return

end

%disp('paused-main'), pause
%return
% end % main if ---------------------------------------

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

        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -2, 12, 0, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];

       mp.ctrl.sched_mat = [...
           repmat(SetC,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetC2,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

a = 'figure(1), print -dpng ../Figs/fig_psfx'

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_06lamD_v4_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_16lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_20lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_22lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v2_7lam_run4_cbx_itr60 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v3_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v2_7lam_run5_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_25lamD_7lam_run5_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_2p5to12lamD_chg4_7lam_run2_cbx_itr78 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_1to12lamD_chg4_7lam_v2_run5_cbx_itr140 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg4_7lam_run3_cbx_itr85 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg6_7lam_run2_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value

if 1

    dm1a = dm1;
    dm2a = dm2;

    bit = 16;
    delz = 2000/(2^bit-1); % nm

    ni = round(dm1a / delz);
    dm1 = ni * delz;

    ni = round(dm2a / delz);
    dm2 = ni * delz;

end

%cbx(1) = 4.0281e-12;
    
ms = length(cbx);

id = 1;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;
    
    figure(11), clf, semilogy(cbx, 'ro-'), grid, 
    if id == 0, cbx = []; end %cbx([ms-1 ms]);

%sub_plot_psf

%disp('paused'), pause

%end

[mp, out, cbx, Im] = falco_wfsc_loop_erkin1_copy(mp, out, cbx, dm1a, dm2a);

figure(1), print -dpng ../Figs/fig_psfx1
figure(103), print -dpng ../Figs/fig_time1

end % ------------------------------------------
%%

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_1to12lamD_chg4_7lam_v2_run5_cbx_itr165 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg4_7lam_run2_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg6_7lam_run2_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value

fs = 14;

mx = [1:length(cbx)] - 1;

s1 = sprintf('last Cb = %0.3g', cbx(end));

    figure(12), clf, semilogy(mx, cbx, 'ro-'), grid, 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('10mm, 0.98D, psd, [2 12] $\lambda$/D','Fontsize',fs,'Interpreter', 'Latex');
        legend(s1, 'Location', 'Northeast')
    print -dpng ../Figs/fig_efc1

%Im = mp.Im;

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;
%Im = falco_get_summed_image(mp);

%sub_plot_psf


[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

lim = 12;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    maskc2 = (rm >= 2) & (rm <= lim);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, lim);
    [xc2, yc2] = fun_circle(xv, 24);

    psfn = Im .* maskc2; aa = psfn;
    cb = mean(nonzeros(psfn(:)))
    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn(:)));
    
    fs = 14;

    figure(10), clf
        imagesc(xv,xv,log10(Im.*maskc), [-11 -8]); axis xy image, colormap(jet); 
        %imagesc(xv,xv,log10(aa), [-12 -8]); axis xy image, colormap(jet); 
        %hold on, plot( xc1, yc1, 'b',  'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('10mm, 2-%i\\lambda/D, 7-\\lambda, psd, r2-itr260: Cb = %0.3g', lim, cb)]);
        axis(16*[-1 1 -1 1])
   print -dpng ../Figs/fig_psf1

   return

    aa = dm1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = dm2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2
   

return
%end % -----------------------------------------------
%%


% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

%clear all;

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_1ring_10mm_cir98_24lamD_7lam_512_run3_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
%tval = tx(end)
%figure(1), clf, semilogy(cbx,'ro-'), grid, return


tic

if 0


%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;

% lamx = mp.sbp_centers in [m]


%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false

% mp.Fend.res = 4; %3; %--Sampling [ pixels per lambda0/D]


%% Step 2: Load default model parameters

%EXAMPLE_defaults_LUVOIRB_VC_design_erkin1
%EXAMPLE_defaults_LUVOIRB_VC_design_erkin2_habex
sub_setup_habex1

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
%pupil = fitsread('/home/esidick/2022habex/dat/pupil_30mm_512pix_sharp_v2.fits');
pupilf = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_soft.fits');
%pupilc = fitsread('/home/esidick/2022habex/dat/pupil_20mm_512pix_sharp.fits');
%pupilf = fitsread('/home/esidick/2022habex/dat/pupil_20mm_2048pix_sharp.fits');
%pupilf = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_1ring.fits');
pupilc = pupilf;
%apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
%apod = pupil * 0 + 1;
lyo   = 'inside falco_flesh_out_workspace';

mp.P1.compact.mask = pupilc; %(load your file here)
mp.P1.full.mask    = pupilf; %(load your file here)

mp.full.flagGenApod = false;
mp.compact.flagGenApod = false;
mp.P3.compact.mask = pupilc * 0 + 1;
mp.P3.full.mask = pupilf * 0 + 1;

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

    [mp, out] = falco_flesh_out_workspace(mp);

%end % -------- end of main ---------------------------------

%save ../dat_02jun2022/dat_main_mp mp

if 0

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_05per_3lam_run1_cbx_itr36 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run1_cbx_itr11 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im 
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im 
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_30mm_7lam_run1_cbx_itr15 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir96_7lam_run1_cbx_itr15 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir95_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_3lamD_7lam_run5_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_24lamD_7lam_run2_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_24lamD_7lam_run2_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_0512_run1_cbx_itr10 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_1024_run1_cbx_itr5 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_1024_run1_cbx_itr5 dm1 dm2 tx cbx fom Im beta_value

id = 1;

    mp.dm1.V = dm1*id;
    mp.dm2.V = dm2*id;
    
    cb0 = cbx(end)
    
ms = length(cbx);

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;
    
    figure(1), clf, semilogy(cbx, 'ro-'), grid, 
    %cbx = []; %cbx([ms-1 ms]);

Im = falco_get_summed_image(mp);

%end
end



%sub_plot_psf

%return

%end

%disp('paused-main'), pause
%end
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
            3, -2, 12, 0, 0;...
            5, -3, 12, 0, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            2, -1, 12, 0, 0;...
            ];
                
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            4, -3, 12, 0, 0;...
            4, -4, 12, 0, 0;...
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            1, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            4, -1, 12, 0, 0;...
            1, -3, 12, 0, 0;...
            4, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];

       mp.ctrl.sched_mat = [...
           repmat(SetD,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetC2,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

a = 'figure(1), print -dpng ../Figs/fig_psfx'

%load /home/esidick/Afalco/falco20200916/dat_erkin/run1_cbx_itr10 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run2_cbx_itr25 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run3_cbx_itr45 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run4_cbx_itr50 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/jeff_run1_cbx_itr6 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/jeff_run2_cbx_itr23 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/dsn2_650_run1_cbx_itr43 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/dsn3_bbtt_run1_cbx_itr73 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/dsn4_bbtt_run1_cbx_itr9 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/dsn4_bbtt_run1_cbx_itr49 dm1 dm2 tx cbx 

%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_05mm_run1_cbx_itr2 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_15mm_7lam_run1_cbx_itr55 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_7lam_run1_cbx_itr30 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_30mm_7lam_run1_cbx_itr55 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_50mm_7lam_run4_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_7lam_run6_cbx_itr80 dm1 dm2 tx cbx fom Im 

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run2_cbx_itr26 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_30mm_v2_7lam_run2_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
%ps = [length(beta_value); length(cbx)], beta_value, return

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_30mm_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir95_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir96_7lam_run1_cbx_itr15 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir97_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir98_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir99_7lam_run4_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir100_7lam_run3_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%cbx1 = cbx;
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_cir96_7lam_run2_cbx_itr15 dm1 dm2 tx cbx fom Im beta_value
%cbx2 = cbx;

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_3lamD_7lam_run3_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_7lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_6lamD_v3_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_12lamD_7lam_run2_cbx_itr20 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_24lamD_v2_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_6lamD_7lam_run5_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_16lamD_7lam_run2_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_22lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_03lamD_v2_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_12lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_03lamD_v3_7lam_run3_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_16lamD_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_512_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_1024_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_1ring_10mm_cir98_24lamD_7lam_512_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_1ring_v3_10mm_cir98_24lamD_7lam_512_run2_cbx_itr150 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_run3_cbx_itr60 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_v2_run2_cbx_itr73 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_3to12lamD_chg4_7lam_run1_cbx_itr25 dm1 dm2 tx cbx fom Im beta_value

ms = length(cbx);

id = 1;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;

    cb0 = cbx(end)
    tt = tx(end)
    
    figure(11), clf, semilogy(cbx, 'ro-'), grid, 
    if id == 0, cbx = []; end %cbx([ms-1 ms]);

%sub_plot_psf

%disp('paused'), pause

%end

[mp, out, cbx] = falco_wfsc_loop_erkin1_copy(mp, out, cbx, dm1a, dm2a);

figure(1), print -dpng ../Figs/fig_psfx1
figure(103), print -dpng ../Figs/fig_time1

end % ------------------------------------------------------------

%Im = mp.Im;
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_run3_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value

fs = 14;

mx = [1:length(cbx)] - 1;

s1 = sprintf('last Cb = %0.3g', cbx(end));

    figure(12), clf, semilogy(mx, cbx, 'ro-'), grid, 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('dark-hole radii = [3 12] $\lambda$/D','Fontsize',fs,'Interpreter', 'Latex');
        legend(s1, 'Location', 'Northeast')
    print -dpng ../Figs/fig_efc1


[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

lim = 12;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 3) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    maskc2 = (rm >= 3) & (rm <= lim);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, lim);
    [xc2, yc2] = fun_circle(xv, 24);

    psfn = Im .* maskc2; aa = psfn;
    cb = mean(nonzeros(psfn(:)))
    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn(:)));
    
    fs = 14;

    figure(10), clf
        imagesc(xv,xv,log10(Im.*maskc2), [-13 -8]); axis xy image, colormap(jet); 
        %imagesc(xv,xv,log10(aa), [-12 -8]); axis xy image, colormap(jet); 
        %hold on, plot(xc, yc, 'w', xc1, yc1, 'w',  'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('10mm, 3-12\\lambda/D, cir-0.98, r1-itr45: Cb = %0.3g', cb)]);
        axis(15*[-1 1 -1 1])
   print -dpng ../Figs/fig_psf1

return
%end % -----------------------------------------------
%%


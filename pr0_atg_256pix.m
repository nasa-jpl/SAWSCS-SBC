% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

%addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;
addpath(genpath('/proj/exospin/users/joegreen/'));

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%tval = tx(end)
%cb0 = cbx(end)
%figure(1), clf, semilogy(cbx,'ro-'), grid, %return

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
%EXAMPLE_defaults_LUVOIRB_VC_design_erkin2_habex
%EXAMPLE_defaults_LUVOIRB_VC_design_erkin3_habex

% -------------------------------
    sub_setup_atg
% -------------------------------

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = true; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;


mp.wfc.model_id = 1;

if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
mp.compact.flagGenPupil = false;

%pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/pup512_ogap3_igap1.fits');
%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/dat_atg/pup256_cir.fits');

lyo   = 'inside setup/falco_flesh_out_workspace';

if 1

%pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

mp.full.flagGenApod = false;
mp.compact.flagGenApod = false;
%mp.P3.compact.mask     = apod * 0 + 1;
%mp.P3.full.mask        = apod * 0 + 1;

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

load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/d200_cir98_2to12lamD_chg6_1lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value 

id = 0;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    
end

%return    
%end % ------------------------------------------
%%


fs = 14;

%return
end % -----------------------------------------------
%%

mp.wfc.model_id = 7; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = atg

load /home/esidick/Afalco/falco20200916/macos/dat_atg/dat_V1_nm V1_nm
load /home/esidick/Afalco/falco20200916/macos/dat_atg/dat_V2_nm V2_nm

%dmx = dm2 * 0;
%dmx(32,32) = 5; % [nm]
    mp.dm1.V = V1_nm * 0;
    mp.dm2.V = V2_nm * 1;

Im = falco_get_summed_image(mp);

    figure(2), clf
        imagesc(log10(Im)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title([sprintf('post-dm1/dm2-wfc: \\sigma = 1nm, Cb = %0.3g ',  cb0)]);
        %axis(12*[-1 1 -1 1])
   %figure(1), print -dpng ../Figs/fig_psfo

return
%end % -----------------------------------------------
%%







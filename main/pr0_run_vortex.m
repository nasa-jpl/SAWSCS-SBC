% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

%clear all;

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;

%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

%EXAMPLE_defaults_try_running_FALCO
%EXAMPLE_defaults_VC_simple_erkin
EXAMPLE_defaults_VC_simple_orig

%% Step 3: Overwrite default values as desired

%%--Special Computational Settings
mp.flagParfor = false; %--whether to use parfor for Jacobian calculation
mp.flagPlot = true;

%--Record Keeping
mp.SeriesNum = 1;
mp.TrialNum = 1;

%--Use just 1 wavelength for initial debugging/testing of code
mp.fracBW = 0.01;       %--fractional bandwidth of the whole bandpass (Delta lambda / lambda0)
mp.Nsbp = 1;            %--Number of sub-bandpasses to divide the whole bandpass into for estimation and control
mp.Nwpsbp = 1;          %--Number of wavelengths to used to approximate an image in each sub-bandpass

mp.Nitr = 3; %--Number of wavefront control iterations

%% Step 4: Generate the label associated with this trial

mp.runLabel = ['Series',num2str(mp.SeriesNum,'%04d'),'_Trial',num2str(mp.TrialNum,'%04d_'),...
    mp.coro,'_',mp.whichPupil,'_',num2str(numel(mp.dm_ind)),'DM',num2str(mp.dm1.Nact),'_z',num2str(mp.d_dm1_dm2),...
    '_IWA',num2str(mp.Fend.corr.Rin),'_OWA',num2str(mp.Fend.corr.Rout),...
    '_',num2str(mp.Nsbp),'lams',num2str(round(1e9*mp.lambda0)),'nm_BW',num2str(mp.fracBW*100),...
    '_',mp.controller];

%return
%% Step 5: Perform the Wavefront Sensing and Control

if 1
    [mp, out] = falco_flesh_out_workspace(mp);
    %save ~esidick/Afalco/falco20200916/dat_erkin/dat_mp mp out
else
    load ~esidick/Afalco/falco20200916/dat_erkin/dat_mp mp out
end

return

[mp, out] = falco_wfsc_loop_erkin(mp, out);

return
% --------------------------------------------------------------
%
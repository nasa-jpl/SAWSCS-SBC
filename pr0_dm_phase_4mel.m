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

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
%tval = tx(end)
%cb0 = cbx(end)
%figure(1), clf, semilogy(cbx,'ro-'), grid, %return

mp.Fend.res = 4;

%tic, sub_wfc, toc, return

idd = 1
%disp('paused'), pause

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
    sub_setup_macos
% -------------------------------

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

%if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
mp.compact.flagGenPupil = false;

%pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/pup512_ogap3_igap1.fits');
%pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
%pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_00mm.fits');
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps_20230530/pupil_gray_gap08.fits');

%pupil = double(pupil > 0);

lyo   = 'inside setup/falco_flesh_out_workspace';

%opdnom_nm = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/opd512_ogap3_igap1_nom_nm.fits');
%load maps/opd_dm1_nm_20221219.mat

%load ~esidick/2023_6mst/macos/dat/dat_amp_opdc amp opdc
%amp_macos = amp; opd_macos_nm = opdc;

%load ~esidick/2023_6mst/6mst/dat/macos_opdnom_20230118 opdnom
%load ~esidick/2023_6mst/6mst/dat/macos_opdnm_20230125 opdnm

%load ~esidick/2022habex/dat_macos/opdnom_256pix_macos opdnm
%dopd_c_nm = opdnm;

%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10 dopd_c_nm 
%load dat_dm_err1/opd30_erkin dopd_c_nm
%load maps_20230530/opdnom_nm_00mm opdnom_nm
%opdnm = opdnm * 1 / 3;
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
%dopd_c_nm = opdnom_nm;
%dopd_c_nm = opdnm;

bandwidth = 32;
rolloff = 0;
% [opd_res, filter] = opdFilter(opd, bandwidth, rolloff)
%[opdlo, filter] = opdFilter(opd, bandwidth, rolloff);
%[opdlo, filter] = opdFilter_Rect(opdnm, bandwidth, rolloff);
%opdhi = opdnm - opdlo;
%opdnm = opdhi;

%opdnm = dopd_c_nm;

mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc

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

    mp.controller = 'plannedefc';

       %SetA = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
       SetA = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
            %repmat([1, 1j, 12, 1, 1], [5, 1]);

        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            3, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            2, -1, 12, 0, 0;...
            ];

        SetB1 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            1, -6, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            4, -1, 12, 0, 0;...
            ];

        SetB2 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            3, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            2, -1, 12, 0, 0;...
            ];

        SetB3 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            8, -1, 12, 0, 0;...
            ];

        SetB4 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            1, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            4, -1, 12, 0, 0;...
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
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            6, -1, 12, 0, 0;...
            ];
        
        SetE = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            ];
if 0
       mp.ctrl.sched_mat = [...
           repmat(SetB,[1,1]);...  % [2, 1] repeat 2 times
           repmat(SetB1,[1,1]);...
           repmat(SetB2,[1,1]);...
           repmat(SetB3,[1,1]);...
           %repmat(SetB4,[1,1]);...
           ];
else
       mp.ctrl.sched_mat = [...
           repmat(SetB1,[1,1]);...  % [2, 1] repeat 2 times
           repmat(SetB,[1,1]);...
           %repmat(SetB3,[1,1]);...
           %repmat(SetB1,[1,1]);...
           %repmat(SetB4,[1,1]);...
           ];
end
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dm1_psd_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr205 dm1 dm2 tx cbx fom Im beta_value

end
    
id = 0;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

flag = 1;

    mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, 4 = dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

Im = falco_get_summed_image(mp);

load ~esidick/Afalco/falco20200916/macos/dat_mel/dat_ep3_ce_at_apod_nom amp opd ce
    
    aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('opd: dm-errors, rms = %0.1fnm', rms0(1)));
   %print -dpng ../Figs/fig_opd_dm_rb1

    aa = amp; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('amp: dm-errors, rms = %0.3g', rms0(1)));
   %print -dpng ../Figs/fig_amp_dm_rb1


return    
%end % ------------------------------------------
%%

if 0

%close all

tc = fom(end)*10;

fs = 14;

mx = [1:length(cbx)] - 1;
s1 = [sprintf('itr%i: last Cb = %0.3g, Tc = %0.2f', mx(end), cbx(end), tc) '%'];

ns = length(cbx);

    figure(11), clf, semilogy(mx, cbx, 'ro-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('gap08, bw20, rms30, flattened, 4\lambda')
        legend(s1, 'Location', 'Northeast')
    %print -dpng ../Figs/fig_efc

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, 3);
    [xc2, yc2] = fun_circle(xv, 12);

end

flag = 1
mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
    %disp('paused'), pause
% ------------------------------------------

if 1
    id_dm = 1
    mp.dm1.V = dm1 * id_dm;
    mp.dm2.V = dm2 * id_dm;

    Im = falco_get_summed_image(mp);

end
% ------------------------------------------

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    fs = 14;
    cx = [-11 -8];
    cx1 = [-7 -4];

%end
%close all

    legx{1} = ['noflat, post 30->10nm-20%-12\lambda: ' sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%'];
    %legx{2} = 'With WF Flattening';
    t1 = ['gap08, bw20, rms30, 1\lambda'];
    t2 = [sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%'];

    figure(12), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title(['post-30nm-20%-12\lambda, ' sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%']);
        title(t1);
        xlabel(t2)
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfo

return

if 0

%Im1 = falco_get_summed_image(mp);

    psfn1 = Im1 .* maskc; 
    cb1 = mean(nonzeros(psfn1(:)));

    figure(13), clf
        imagesc(xv,xv,log10(psfn1)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title(['pre-30nm-20%-12\lambda-flat, ' sprintf('Cb = %0.3g',  cb1)]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfi
end

   %return
   %end

%load dat_dm_err1/dm_err_1nm_uniform dm1x dm2x

    %close all

    aa = dm1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = dm2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2

   ns = length(bvalo);
   ni = length(cbx) - ns + 1;
   nx = [ni:length(cbx)];

   return

    figure(10), clf, semilogy(nx, bvalo, 'ro-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('\beta-value');
        title('rms = 00nm, bw = 20%, 4\lambda')
    print -dpng ../Figs/fig_bval

return
%end % -----------------------------------------------
%%




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
disp('paused'), pause

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

lyo   = 'inside setup/falco_flesh_out_workspace';

%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10 dopd_c_nm 
%load dat_dm_err1/opd30_erkin dopd_c_nm
%load maps_20230530/opdnom_nm_00mm opdnom_nm
load /home/esidick/Afalco/falco20200916/macos/maps_20230530/opdnm_30mm opdnm %psdnm opdnm_only  
%dopd_c_nm = opdnom_nm;
%dopd_c_nm = opdnm;

%load maps_psd/post_efc_opd_amp_30nm_bw20_1500nm_gray08_owa12 opd_target amp_target ce_target
%opdnm = opd_target;
%pupil = amp_target;

if 0
load /home/esidick/Afalco/falco20200916/dat_bb_256pix/bb_bw20_30nm_2to12lamD_chg6_4lam_itr300 dm1 dm2  
load dat_dm_err1/wfc_dm1_flat_30nm_500nm dm1x opdx
dm1_save = dm1x(:,:,end);

    aa = opdnm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    %aa = opdx(:,:,end); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    %aa = dm1-dm1_save; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        %imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        imagesc(2*pupil-double((aa~=0))); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('post30: rms = %0.1f, pv = %0.1fnm, rat = %0.2f', rms0, rms0(2)/rms0(1)));
    %print -dpng ../Figs/fig_opd30a
   return
end

%opdnm = dopd_c_nm;

mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc

pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
%pupil = ce_target;

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
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            3, -1, 12, 0, 0;...
            ];

        SetB1 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            1, -6, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            4, -1, 12, 0, 0;...
            ];

        SetB2 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            2, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            3, -1, 12, 0, 0;...
            ];

        SetB3 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            1, -6, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            4, -1, 12, 0, 0;...
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

       mp.ctrl.sched_mat = [...
           repmat(SetB1,[1,1]);...  % [2, 1] repeat 2 times
           repmat(SetB,[1,1]);...
           repmat(SetB3,[1,1]);...
           repmat(SetB2,[1,1]);...
           %repmat(SetB4,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%end

%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_00mm_bw20_30nm_2to12lamD_chg6_4lam_itr400 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_16mm_bw20_30nm_2to12lamD_chg6_4lam_itr420 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_24mm_bw20_30nm_2to12lamD_chg6_4lam_itr440 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/sharp_gap_24mm_bw20_30nm_2to12lamD_chg6_4lam_itr360 dm1 dm2 tx cbx fom Im bvalo 

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to16lamD_chg6_4lam_itr370 dm1 dm2 tx cbx fom Im bvalo 
%load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to20lamD_chg6_4lam_itr420 dm1 dm2 tx cbx fom Im bvalo 

%bvali = beta_value;
bvali = bvalo;

id = 1

    %mp.dm1.V = dm1 * id + dm1_save;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    %mp.dm1.V = dm1_save;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;

    figure(11), clf, semilogy(cbx, 'ro-'), grid, %disp('paused'), pause

    if id == 0, cbx = []; bvali = []; end %cbx([ms-1 ms]);

%sub_plot_psf

%disp('paused'), pause

mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

%[mp, out, cbx, Im, bvalo] = falco_wfsc_loop_erkin1(mp, out, cbx, dm1a, dm2a, bvali);
%fom = out.thput * 10;

%out.thput = mp.thput_vec;
%out.Nitr = mp.Nitr;
%out.InormHist = InormHist;


%figure(1), print -dpng ../Figs/fig_psfx
%figure(103), print -dpng ../Figs/fig_time


a = 'figure(1), print -dpng ../Figs/fig_psfx'

%tic, sub_target_amp_opd, toc, return
%tic, sub_wfc, toc, return
%tic, sub_wfc_flatten, toc, return
%tic, sub_psf_eval, toc, return

%return    
end % ------------------------------------------
if 1
%%

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 

tc = fom(end)*10;

fs = 14;

mx = [1:length(cbx)] - 1;
s1 = [sprintf('itr%i: last Cb = %0.3g, Tc = %0.2f', mx(end), cbx(end), tc) '%'];

ns = length(cbx);

    figure(11), clf, semilogy(mx, cbx, 'ro-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('gap08, bw20, rms30, owa16, 4\lambda')
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

flag = 7
mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
    %disp('paused'), pause
% ------------------------------------------

if 1
    id = 0;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    Im = falco_get_summed_image(mp);

    id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

%else

    dm1 = mp.dm1.V;
    dm2 = mp.dm2.V;

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
    t1 = ['from ce-500: gap08, bw20, rms30, 4\lambda'];
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
   print -dpng ../Figs/fig_psfo

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
load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case2_20nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case3_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 

load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw05_10nm_2to12lamD_chg6_3lam_itr255 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw10_10nm_2to12lamD_chg6_4lam_itr350 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw15_10nm_2to12lamD_chg6_4lam_itr250 dm1 dm2 tx cbx fom Im beta_value 
%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/bb_bw20_10nm_2to12lamD_chg6_4lam_itr150 dm1 dm2 tx cbx fom Im beta_value 


dm1_save = dm1;
dm2_save = dm2;

pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
%load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10 dopd_c_nm 
load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
load dat_dm_err1/opd30_erkin dopd_c_nm

opdnm = dopd_c_nm + opd_rb_drift_100pm * 1e-3;

pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

%%
%end

if 0

    mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

    load dat_dm_err1/wfc_dm1_dm2_case6_30nm dm1x dm2x opdx ampx
        dm1 = dm1x(:,:,end);
        dm2 = dm2x(:,:,end);

    load dat_dm_err1/dm_err_1nm dm1x dm2x
        dm1 = dm1_save + dm1x;
        dm2 = dm2_save + dm2x;

id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

%Im = falco_get_summed_image(mp);

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    fs = 14;
    cx = [-12 -8];
    %cx = [-10 -6];

    %close all

    figure(15), clf
        imagesc(xv,xv,log10(psfn),cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title([sprintf('post-wfc-30nm: Cb = %0.3g ',  cb0)]);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfi3

return
end % -----------------------------------------------
%
if 0
    %load IFopd/G_dm1_lim5_m1 G km1
    %load IFopd/G_dm2_lim5_m1 G km2
    load IFopd/G_dm1_dm2_lim5_m2_1064nm G km1 km2 indx
    %load IFopd/dwddm_dm1_nm_over_nm_full indx maskopd
    mp.wfc.indx = indx;
    mp.wfc.km1 = km1;
    mp.wfc.km2 = km2;
    mp.wfc.G = G;
end

%return
%end % -----------------------------------------------
%%

%load IFopd/du_dm1_dm2_20230126 dm1 dm2 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd opd_target amp_target
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm opd_target amp_target 
%load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm_355nm opd_target amp_target 
load ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_case6_30nm_1064nm opd_target amp_target 

mp.wfc.opd_target = opd_target;
mp.wfc.amp_target = amp_target;

    aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-opd: rms = %0.2f, pv = %0.1fnm', rms0));

    aa = amp_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa); axis xy image, colormap(jet);  colorbar
        title(sprintf('target-amp: rms = %0.2f, pv = %0.1f', rms0));
    %return


%mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_opdx_10nm';
mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1064nm';

%end


%load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/d200_cir98_2to12lamD_chg6_1lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value 
%load dat_dm_err1/dm_err_10nm dm1x dm2x

if 0
    ma = 64;
    del = 1; % nm
    dm1x = del * randn(ma,ma);
    dm2x = del * randn(ma,ma);
    dm1x(find(dm1x >  del)) = del;
    dm1x(find(dm1x < -del)) =-del;
    dm2x(find(dm2x >  del)) = del;
    dm2x(find(dm2x < -del)) =-del;
    %save dat_dm_err1/dm_err_1nm dm1x dm2x
else
    load dat_dm_err1/dm_err_1nm dm1x dm2x
end

dm1 = dm1_save + dm1x;
dm2 = dm2_save + dm2x;

%end

%load dat_dm_err1/wfc_dm1_dm2_case2_20nm dm1x dm2x opdx ampx
%dm1 = dm1x(:,:,end);
%dm2 = dm2x(:,:,end);

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;

mp.wfc.model_id = 6; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

    Im = falco_get_summed_image(mp);

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    %end

    cx = [-12 -8];

    figure(2), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title(sprintf('post-wfc-10nm: \\sigma = 1nm, Cb = %0.3g ',  cb0));
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfo1

   return

    aa = dm1x; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('\\Deltadm1: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1

    aa = dm2x; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('\\Deltadm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2

return
%end % -----------------------------------------------
%%

close all

mp.wfc.model_id = 1; % =1 stantard, =2 dm1_opd, =3 dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp

load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value 
dm1_save = dm1;
dm2_save = dm2;

load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/segDrift_6mst_0323 opd_rb_drift_100pm
load dat_dm_err1/opd30_erkin dopd_c_nm

%opdnm = dopd_c_nm;
opdnm = dopd_c_nm + opd_rb_drift_100pm * 1e-3;

load dat_dm_err1/dm_err_1nm dm1x dm2x
del_dm1 = dm1x;
del_dm2 = dm2x;

%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm dm1x dm2x opdx ampx
load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_355nm dm1x dm2x opdx ampx
%load /home/esidick/Afalco/falco20200916/macos/dat_dm_err1/wfc_dm1_dm2_case7c_30nm_1064nm dm1x dm2x opdx ampx


cc = 1;

if 0 %cc == 1
    %load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case4_10nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
    %load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd10_erkin dopd_c_nm 
    load dat_dm_err1/opd10_erkin dopd_c_nm
    load dat_dm_err1/wfc_dm1_dm2_case4_10nm dm1x dm2x opdx ampx
%elseif cc == 2
    load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case5_20nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
    %load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd20_erkin dopd_c_nm 
    load dat_dm_err1/opd20_erkin dopd_c_nm
    load dat_dm_err1/wfc_dm1_dm2_case5_20nm dm1x dm2x opdx ampx
%else
    load /home/esidick/Afalco/falco20200916/dat_macos_20230209_256pix/case6_30nm_cir98_2to12lamD_chg6_1lam_run1_cbx_itr100 dm1 dm2  
    %load /proj/jwst2/jzlou/git-6MST/old-6MST-jzlouStuff/macos/opd_wfc_4Erkin/opd30_erkin dopd_c_nm 
    load dat_dm_err1/opd30_erkin dopd_c_nm
    load dat_dm_err1/wfc_dm1_dm2_case6_30nm dm1x dm2x opdx ampx
end

    pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
    %opdnm = dopd_c_nm;
    pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

    ns = size(opdx,3);

    %mp.dm1.V = dm1_save + del_dm1;
    %mp.dm2.V = dm2_save + del_dm2;
    mp.dm1.V = dm1x(:,:,end);
    mp.dm2.V = dm2x(:,:,end);

tic
    Im = falco_get_summed_image(mp);
toc

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    %end

    fs = 14;

    cx = [-9 -5];
    cx = [-12 -8];

    figure(2), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('post-355: dm/rb-errors, \\sigma = 1nm, Cb = %0.3g ',  cb0)]);
        axis(12*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfi2a
   return

   mx  = [0];
   cbx = [cb0];
   rms_opd = [];
   rms_amp = [];

   clear psfx

   for ii = 1:ns

        mp.dm1.V = dm1x(:,:,ii);
        mp.dm2.V = dm2x(:,:,ii);

        Im = falco_get_summed_image(mp);
        psfn = Im .* maskc; 
        cbi = mean(nonzeros(psfn(:)));

        psfx(:,:,ii) = psfn;

        mx = [mx ii];
        cbx = [cbx cbi];

        aa = opdx(:,:,ii); rms_opd = [rms_opd; stat2d(aa)];
        aa = ampx(:,:,ii); rms_amp = [rms_amp; stat2d(aa)];

        figure(1), clf, semilogy(mx, cbx, 'ro-'), grid, drawnow

   end

    save dat_dm_err1/post_wfc_case6_30nm mx cbx rms_opd rms_amp psfx

return
%end % -----------------------------------------------
%%








load /home/esidick/Afalco/falco20200916/dat_macos_20230125/sixMST_cir97_2to12lamD_chg6_1lam_run1_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value

load IFopd/dwddm_dm1_nm_over_nm_full indx maskopd
mp.wfc.indx = indx;
mp.wfc.maskopd = maskopd;

load IFopd/target_opd target_opd
%load ~esidick/2023_6mst/6mst/dat/macos_opdnm_20230125 opdnm
load dat/dat_ep3_opdnm_target_20230125 opdnm_target

    mp.dm1.V = dm1 * 0;
    mp.dm2.V = dm2 * 0;

    mp.wfc.opd_target = target_opd * 0;

    load ~esidick/Afalco/falco20200916/macos/dat_dm_err/test_opd opd
    aa = real(opd); rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
    
%return
%end % -----------------------------------------------
%%

load ~esidick/2023_6mst/6mst/dat/macos_opdnm_20230125 opdnm
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/pup512_ogap3_igap1.fits');
pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

mp.wfc.model_id = 3; % =1 stantadr, =2 dm1_opd, =3 dm1_wfc
%load dat_dm_err/wfc_dmx_opdx_01nm dmx dm2x opdx
load dat_dm_err/flatten_init_wfc_dmx_opdx dmx opdx
mp.dm1.V = dmx(:,:,end)*0;
mp.dm2.V = dm2*0;

%load dat_dm_err/input_dmx_01nm dm1x dm2x
    %save dat_dm_err/input_dmx_06nm dm1x dm2x del

fac = 1;
del = 0.4;

if 0
    dm1x = dm1 + del * rand(size(dm1i));
    dm2x = dm2 + del * rand(size(dm2i));

    mp.dm1.V = dm1x * fac;
    mp.dm2.V = dm2x * fac;

    %save dat_dm_err/input_dmx_01nm dm1x dm2x
    %save dat_dm_err/input_dmx_04nm dm1x dm2x del
end

    %mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dmx_opdx_04nm'; % to save wfc dm1 commands
    mp.wfc.fn = '~esidick/Afalco/falco20200916/macos/dat_dm_err/flatten_init_wfc_dmx_opdx'; % to save wfc dm1 commands

    Imo = falco_get_summed_image(mp);

    %load ~esidick/Afalco/falco20200916/macos/IFopd/dm1_opdtarget dm1 opd_target
    
    psfn = Imo .* maskc; 
    cb = mean(nonzeros(psfn(:)))

    psfn1 = Imo .* maskc1;
    cb1 = mean(nonzeros(psfn1(:)));

%end
cx = [-12 -8];

    figure(2), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('post-wf-flatten: \\sigma = %0.1fnm, 1\\lambda, Cb = %0.3g', del, cb)]);
        title([sprintf('post-wf-flatten: 1\\lambda, Cb = %0.3g', cb)]);
        axis(12*[-1 1 -1 1])
   figure(2), print -dpng ../Figs/fig_psf

return

   dm1 = dm1b;
   dm2 = dm2b;

   return
   %end
    aa = (dm1b+dm2b)*2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title('(dm1 + dm2) x 2')
        xlabel(sprintf('rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm1
return

%end
%load IFopd/target_opd target_opd
%opdnm = target_opd;
load ~esidick/2023_6mst/6mst/dat/macos_opdnm_20230125 opdnm
pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/pup512_ogap3_igap1.fits');
pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);
mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

    figure(2), clf, imagesc([target_opd opdnm]), axis xy image, colormap(jet), colorbar, niceaxes, 
%return

mp.wfc.model_id = 1; % =1 stantadr, =2 dm1_opd, =3 dm1_wfc


load /home/esidick/Afalco/falco20200916/dat_macos_20230125/sixMST_cir97_2to12lamD_chg6_1lam_run1_cbx_itr50 dm1 dm2 cbx
cb0 = cbx(end); 
load dat_dm_err/input_dmx_08nm dm1x dm2x %del
    mp.dm1.V = dm1x;
    mp.dm2.V = dm2x;

    del = 0.1;

    Imo = falco_get_summed_image(mp);
    psfn = Imo .* maskc; 
    cb = mean(nonzeros(psfn(:)));

load dat_dm_err/wfc_dmx_opdx_08nm dmx %dm1x dm2x opdx opd_in

%load dat_dm_err/wfc_dmx_opdx_02nm dmx opdx
%load IFopd/wfc_dmx_opdx dmx opdx
ns = size(opdx,3);


jj = 0;

figure(2), clf
    figure(2), subplot(4,6,jj+1), imagesc(log10(psfn+eps), [-12 -8]), axis xy image, colormap(jet), colorbar, niceaxes, 
    %figure(2), clf, imagesc(log10(psfn+eps), [-12 -8]), axis xy image, colormap(jet), colorbar, niceaxes, 
        title(sprintf('%i: Cb = %0.3g', jj, cb)), drawnow

mx = [0];
cbx= [cb];

for jj = 1: ns

    mp.dm1.V = dmx(:,:,jj);

    Imo = falco_get_summed_image(mp);

    psfn = Imo .* maskc; 
    cb = mean(nonzeros(psfn(:)));

    cbx = [cbx cb];

    figure(2), subplot(4,6,jj+1), imagesc(log10(psfn+eps), [-12 -8]), axis xy image, colormap(jet), colorbar, niceaxes, 
    %figure(2), clf, imagesc(log10(psfn+eps), [-12 -8]), axis xy image, colormap(jet), colorbar, niceaxes, 
        title(sprintf('%i: Cb = %0.3g', jj, cb)), drawnow
end
   %print -dpng ../Figs/fig_psf

%save dat_dm_err/cbx_08nm cbx

%end

mx = [1:length(cbx)] - 1;
%mx = [0 mx];
%cbx = [8.6e-10 cbx];

figure(10), clf, semilogy(mx, cbx, 'ro-', mx, mx*0+cb0,'b'), grid, 
    parms = fun_newaxes(14, 2, 5);  
    xlabel('wfc iterations')
    title('\sigma = 0.6nm: Mean Contrast ')
    legend('data','nominal','Location','Northwest')
    axis([0 20 7e-11 1e-9])
%print -dpng ../Figs/fig_cb06

return
%end % -----------------------------------------------
%%







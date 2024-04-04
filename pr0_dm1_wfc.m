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

if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
mp.compact.flagGenPupil = false;

pupil = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/pup512_ogap3_igap1.fits');

lyo   = 'inside setup/falco_flesh_out_workspace';

%opdnom_nm = fitsread('/home/esidick/Afalco/falco20200916/macos/maps/opd512_ogap3_igap1_nom_nm.fits');
%load maps/opd_dm1_nm_20221219.mat

%load ~esidick/2023_6mst/macos/dat/dat_amp_opdc amp opdc
%amp_macos = amp; opd_macos_nm = opdc;

%load ~esidick/2023_6mst/6mst/dat/macos_opdnom_20230118 opdnom
load ~esidick/2023_6mst/6mst/dat/macos_opdnm_20230125 opdnm

load maps/du_dm1_dm2_nm_20221219.mat dm1 dm2 
dm1i = dm1 * 1;
dm2i = dm2 * 1;
%opdnom_nm = opd_dm1_nm;
opdnom_nm = opdnm;

load ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_opdnm_target_20230125 opdnm_target

opdnm_macos  = opdnm;
opdnm_target = real(opdnm_target);

    aa = opdnm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('macos-opdnm: no dm, no diff'));
        xlabel(sprintf('rms = %0.2f, pv = %0.1fnm', rms0));
   %print -dpng ../Figs/fig_opd1

    aa = opdnm_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(aa); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('target: dm1-dm2 distance = 0.2m'));
        xlabel(sprintf('rms = %0.2f, pv = %0.1fnm', rms0));
   %print -dpng ../Figs/fig_opd1
%return

mp.wfc.model_id = 1; % =1 stantadr, =2 dm1_opd, =3 dm1_wfc



if 1



if 0
np = 512;
load ~esidick/Afalco/falco20200916/macos/IFopd/dat_EP3_target EP3
lambda = 500e-9;

    ce = pad(EP3,np);
    pha_nm = 1e9 * atan2(imag(ce),real(ce)) * lambda / 2 / pi;
    opd_target = pha_nm .* pupil;

    %save ~esidick/Afalco/falco20200916/macos/IFopd/dm1_opdtarget dm1 opd_target

    del = 0.1; % nm
    dm1x = dm1i + del * rand(size(dm1i));
    dm2x = dm2i + del * rand(size(dm2i));

end

id = 0;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    mp.model_id = 1;

%save maps/du_dm1_dm2_nm_opd_20221219.mat dm1 dm2 opd_dm1_nm
%return



%load maps/opd_dm_post_maint_4falco 
%opdnom_nm = opd_dm1_new_nm;

%rms0 = stat2d(opdx) *1e9
%dm1_save = dm1 * 1e9;
%figure(1), clf, imagesc(opdx), axis xy image, colormap(jet), colorbar
%figure(2), clf, imagesc(dm1_save), axis xy image, colormap(jet), colorbar
%return

%pupil = pupil .* exp(1i*2*2*pi*opdx/mp.lambda0);
pupil = pupil .* exp(1i*2*pi*opdnm*1e-9/mp.lambda0);

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

    mp.controller = 'plannedefc';

       SetA = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
            %repmat([1, 1j, 12, 1, 1], [5, 1]);

        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -2, 12, 0, 0;...
            5, -3, 12, 0, 0;...
            5, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -2, 12, 0, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetD = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            ];

       mp.ctrl.sched_mat = [...
           repmat(SetC,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetC2,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_20230125/sixMST_cir97_2to12lamD_chg6_1lam_run1_cbx_itr10 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_macos_20230125/d200_cir97_2to12lamD_chg6_1lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value

id = 1;

    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    %mp.dm1.V = dm1_save;

    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;
    
    figure(11), clf, semilogy(cbx, 'ro-'), grid, %pause
    if id == 0, cbx = []; end %cbx([ms-1 ms]);

%sub_plot_psf

%disp('paused'), pause

%end

%[mp, out, cbx, Im] = falco_wfsc_loop_erkin1(mp, out, cbx, dm1a, dm2a);

figure(1), print -dpng ../Figs/fig_psfx
figure(103), print -dpng ../Figs/fig_time


a = 'figure(1), print -dpng ../Figs/fig_psfx'

%end

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/opdz_30mm_cir98_rms_30nm_2to12lamD_chg6_4lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/nom_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr10 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
    
    figure(11), clf, semilogy(cbx, 'ro-'), grid, %pause

%return    
end % ------------------------------------------
%%

%load /home/esidick/Afalco/falco20200916/dat_macos_20230125/sixMST_cir97_2to12lamD_chg6_1lam_run1_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value

%save IFopd/du_dm1_dm2_20230126 dm1 dm2 opdnm_macos opdnm_target

fs = 14;

tc = fom(end)*10;

mx = [1:length(cbx)] - 1;
s1 = [sprintf('last Cb = %0.3g, Tc = %0.1f', cbx(end), tc) '%'];

ns = length(cbx);

    figure(1), clf, semilogy(mx, cbx, 'ro-'), grid, 
    %figure(12), clf, semilogy(mx1, cbx(mx1+1), 'ro-', mx2, cbx(mx2+1), 'bs-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        %title('10mm, 0.98D, $\Delta$psd-$\lambda$/100, chrg4, [2 12] $\lambda$/D, ','Fontsize',fs,'Interpreter', 'Latex');
        %title('opdz-rms=30nm, 30mm, 0.98D','Fontsize',fs,'Interpreter', 'Latex');
        %title('gap = 35/12mm, 0.98D, 1\lambda','Fontsize',fs,'Interpreter', 'Latex');
        %title('Nominal','Fontsize',fs,'Interpreter', 'Latex');
        %title('10mm, 0.98D, no psd, 4$\lambda$-efc, chrg6, [2 12] $\lambda$/D ','Fontsize',fs,'Interpreter', 'Latex');
        title('Cb vs Itr: dm1/dm2-d = 0.2m')
        legend(s1, 'Location', 'Northeast')
    %print -dpng ../Figs/fig_efc

if 0
    %load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pupil_512pix_gap30mm pupil
    load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pup1_mono_512pix pupil
    %load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
    load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo1000 del_psd psd_nm indx_opd dm1

    %pupil = pupil .* exp(1i*2*pi*del_psd*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

end

%return
end % -----------------------------------------------
%%

load IFopd/du_dm1_dm2_20230126 dm1 dm2 

mp.wfc.model_id = 1; % =1 stantadr, =2 dm1_opd, =3 dm1_wfc

dmx = dm2 * 0;
dmx(32,32) = 5; % [nm]

    mp.dm1.V = dm1 * 0;
    mp.dm2.V = dmx * 1;

Im = falco_get_summed_image(mp);

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, 3);
    [xc2, yc2] = fun_circle(xv, 12);

    psfn = Im .* maskc; 
    cb0 = mean(nonzeros(psfn(:)));

    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn1(:)));

    %tc = fom(end)*10;

    %save ../dat_paper/fig26_Im_5c_9s_lamo100_5lam_30mm Im xv cb maskc 
    %save ../dat_paper/fig26_Im_5c_9s_lamo100_5lam_30mm_v2 Im xv cb maskc 

%end
    
    fs = 14;
    cx = [-12 -8];

    figure(2), clf
        imagesc(xv,xv,log10(psfn)); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        title([sprintf('post-efc: 200mm, 1\\lambda, Cb = %0.3g ',  cb0)]);
        axis(12*[-1 1 -1 1])
   %figure(1), print -dpng ../Figs/fig_psfo

   return

    aa = dm1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(13), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %title(sprintf('dm1+dm2: rms = %0.2f, pv = %0.1fnm', rms0));
        title(sprintf('dm1: rms = %0.2f, pv = %0.1fnm', rms0));
   %print -dpng ../Figs/fig_dm1

    aa = dm2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   %print -dpng ../Figs/fig_dm2

return
%end % -----------------------------------------------
%%


%load ~esidick/Afalco/falco20200916/macos/IFopd/G_matrix G
%load IFopd/G_matrix G km1
load IFopd/G_matrix_m1 km1
load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full dwddm indx maskopd
M = dwddm(:,km1);  % nm/nm
clear dwddm1
        
    A = M'*M;

    sv = diag(A);
    jmax = max(sv);
    sn = sv / jmax;
    
    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;

%end

    gi  = 0; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

%end

%load IFopd1/G_matrix G km1
load IFopd/G_matrix_lim5_m1_v2 G km1

mp.wfc.G = G;
mp.wfc.km1 = km1;
clear G

%return
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

%end






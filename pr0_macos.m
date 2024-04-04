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
load maps/opd_dm1_nm_20221219.mat

%load ~esidick/2023_6mst/macos/dat/dat_amp_opdc amp opdc
%amp_macos = amp; opd_macos_nm = opdc;

%load ~esidick/2023_6mst/6mst/dat/macos_opdnom_20230118 opdnom

load maps/du_dm1_dm2_nm_20221219.mat dm1 dm2 
dm1b = dm1 * 1;
dm2b = dm2 * 1;
opdnom_nm = opd_dm1_nm;
%opdnom_nm = opdnom;

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
pupil = pupil .* exp(1i*2*pi*opdnom_nm*1e-9/mp.lambda0);

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
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -1, 12, 0, 0;...
            ];
        
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -2, 12, 0, 0;...
            5, -3, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            5, -4, 12, 0, 0;...
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

a = 'figure(1), print -dpng ../Figs/fig_psfx'

%end

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/opdz_30mm_cir98_rms_30nm_2to12lamD_chg6_4lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/nom_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr10 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec19_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr40 dm1 dm2 tx cbx fom Im beta_value

dm1 = dm1b;
dm2 = dm2b;
%dm1 = dm1_new_nm_4falco;
%dm2 = dm2_new_nm_4falco;

if 0

    dm1a = dm1;
    dm2a = dm2;

    bit = 16;
    delz = 2000/(2^bit-1); % nm

    ni = round(dm1a / delz);
    dm1 = ni * delz;

    ni = round(dm2a / delz);
    dm2 = ni * delz;

end

    
ms = length(cbx);

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

%[mp, out, cbx, Im, dm1, dm2, fom] = falco_wfsc_loop_erkin1(mp, out, cbx, dm1a, dm2a);

%figure(1), print -dpng ../Figs/fig_psfx
%figure(103), print -dpng ../Figs/fig_time

%save dat/dat_mp_20221220 mp
%return

%end % ------------------------------------------
%%
%load /home/esidick/Afalco/falco20200916/dat_macos_part1/dec28_30mm_cir97_2to12lamD_chg6_1lam_run1_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value

fs = 14;

tc = fom(end)*10;

mx = [1:length(cbx)] - 1;
s1 = [sprintf('last Cb = %0.3g, Tc = %0.1f', cbx(end), tc) '%'];

%mx1 = mx(1:261);
%mx2 = mx(262:end);

%s1 = sprintf('no digit-errors: last Cb = %0.3g', cbx(261));
%s2 = sprintf(' with digit-errors:last Cb = %0.3g', cbx(end));
ns = length(cbx);
%kk = 152:ns;
%mx = [1:length(kk)];

    figure(12), clf, semilogy(mx, cbx, 'ro-'), grid, 
    %figure(12), clf, semilogy(mx1, cbx(mx1+1), 'ro-', mx2, cbx(mx2+1), 'bs-'), grid, 
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        %title('10mm, 0.98D, $\Delta$psd-$\lambda$/100, chrg4, [2 12] $\lambda$/D, ','Fontsize',fs,'Interpreter', 'Latex');
        %title('opdz-rms=30nm, 30mm, 0.98D','Fontsize',fs,'Interpreter', 'Latex');
        %title('gap = 35/12mm, 0.98D, 1\lambda','Fontsize',fs,'Interpreter', 'Latex');
        %title('Nominal','Fontsize',fs,'Interpreter', 'Latex');
        %title('10mm, 0.98D, no psd, 4$\lambda$-efc, chrg6, [2 12] $\lambda$/D ','Fontsize',fs,'Interpreter', 'Latex');
        title('EFC of 2022-12-28: With Drifts')
        legend(s1, 'Location', 'Northeast')
    %print -dpng ../Figs/fig_efc

%Im = mp.Im;


%lyot  = fitsread('/home/esidick/2022habex/dat/lyot_cir98_512pix_soft.fits');
%load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/lyo_512pix_098D lyo
 %   lyot = lyo;
%mp.P4.compact.mask = lyot; %(load your file here)
%mp.P4.full.mask    = lyot; %(load your file here)

%save ../dat_bijan/dat_mp_30mm_cir98_no_psd_2to12lamD_chg4_5lam mp maskc  

if 0
    %load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pupil_512pix_gap30mm pupil
    load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/pup1_mono_512pix pupil
    %load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
    load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo1000 del_psd psd_nm indx_opd dm1

    %pupil = pupil .* exp(1i*2*pi*del_psd*1e-9/mp.lambda0);

    mp.P1.compact.mask = pupil; %(load your file here)
    mp.P1.full.mask    = pupil; %(load your file here)

end

end

    mp.dm1.V = dm1b * 0;
    mp.dm2.V = dm2b * 0;

Im = falco_get_summed_image(mp);

%sub_plot_psf

%end

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

%end

lim = 12;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    %maskc2 = (rm >= 2) & (rm <= lim);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, 3);
    [xc2, yc2] = fun_circle(xv, 12);

    %end

    psfn = Im .* maskc; 
    cb = mean(nonzeros(psfn(:)));
    %cba = cbx(end);

    %pcb = [cba cb]

    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn1(:)));

    tc = fom(end)*10;

    %save ../dat_paper/fig26_Im_5c_9s_lamo100_5lam_30mm Im xv cb maskc 
    %save ../dat_paper/fig26_Im_5c_9s_lamo100_5lam_30mm_v2 Im xv cb maskc 

%end
    
    fs = 14;
    cx = [-12 -8];

    figure(10), clf
        imagesc(xv,xv,log10(Im.*maskc), cx); axis xy image, colormap(jet); 
        %imagesc(xv,xv,log10(Im)); axis xy image, colormap(jet); 
        %hold on, plot(xc,yc,'b', xc1,yc1,'b', xc2,yc2,'b')
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('post-efc: Cs(2.5-3.5) = %0.3g, Cb(2-12) = %0.3g', cb1, cb) ]);
        %title([sprintf('1\\lambda, Csb = [%0.2f %0.2f]*1e-10, tc = %0.1f', cb1*1e10, cb*1e10, tc) '%']);
        title([sprintf('post-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        axis(12*[-1 1 -1 1])
   figure(10), print -dpng ../Figs/fig_psf

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

    aa = dm2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(14), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('dm2: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_dm2

    aa = opdnom_nm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(15), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('input-opd: rms = %0.2f, pv = %0.1fnm', rms0));
   print -dpng ../Figs/fig_opd
   
   %save maps/du_dm1_dm2_nm_20221219 dm1 dm2
   %save maps/du_dm1_dm2_nm_20221228 dm1 dm2 opd_dm1_new_nm

return
%end % -----------------------------------------------
%%


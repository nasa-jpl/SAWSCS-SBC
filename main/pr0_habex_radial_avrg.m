% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

%clear all;

tic

if 1

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
pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_sharp_v2.fits');
apod  = fitsread('/home/esidick/2022habex/dat/apod.fits');
lyo   = 'inside falco_flesh_out_workspace';

mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)

mp.full.flagGenApod = false;
mp.compact.flagGenApod = false;
mp.P3.compact.mask = apod * 0 + 1;
mp.P3.full.mask = apod * 0 + 1;

% for apodizer, 
%mp.flagApod = false;    %--Whether to use an apodizer or not
mp.flagApod = true;    %--Whether to use an apodizer or not

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

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im 
%    mp.dm1.V = dm1*1;
%    mp.dm2.V = dm2*1;

%tic

    Im = falco_get_summed_image(mp);

fs = 16;

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2.4) & (rm <= 24);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    
    [xc, yc]  = fun_circle(xv, 2.4);
    [xc1, yc1]  = fun_circle(xv, 3.5);
    [xc2, yc2] = fun_circle(xv, 24);
    
    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));
    psfn = Im .* maskc1;
    cb1 = mean(nonzeros(psfn(:)));
    
    psfm = Im; % .* maskc1;

end

%[rx, yy, ynom] = fun_rms_contrast(psf0, psfx, rm, lim1, lim2)
[rx, yy, ynom] = fun_rms_contrast(Im, Im+1, rm, 2, 25);

    figure(1), clf
        imagesc(xv,xv,log10(aa),[-11 -8]); axis xy image, colormap(jet); 
        %hold on, plot(xc, yc, 'w', xc1, yc1, 'w', xc2, yc2, 'w', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %caxis([-11 -6]);                
        %xlabel('Field-Angle [$\lambda_{central}$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title(['7\lambda-15%: ' sprintf('Cs = %0.3g, Cb = %0.3g', cb1, cb)]);
   %print -dpng ../Figs/fig_psf1b

cb3 = ynom;

save ../dat_bb_wfc2/mean_cb_vs_angle_10mm_20per rx ynom

figure(2), clf, semilogy(rx, ynom, 'r-'), grid
toc

return

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

if 1         

    mp.controller = 'plannedefc';

       SetA = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
            %repmat([1, 1j, 12, 1, 1], [5, 1]);

        SetB = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            2, -2, 12, 0, 0;...
            3, -3, 12, 0, 0;...
            5, -4, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            ];
        
        SetC = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -3, 12, 0, 0;...
            3, -4, 12, 0, 0;...
            2, -5, 12, 0, 0;...
            8, -2, 12, 0, 0;...
            ];
        
       mp.ctrl.sched_mat = [...
           repmat(SetB,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetC2,[1,1]);...
           ];
%end
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%end

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
load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im 
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_15per_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
%ps = [length(beta_value); length(cbx)], beta_value, return

ms = length(cbx);

    mp.dm1.V = dm1 * 1;
    mp.dm2.V = dm2 * 1;
    dm1a = mp.dm1.V;
    dm2a = mp.dm2.V;

    %cbx = []; %cbx([ms-1 ms]);
    
%end

[mp, out] = falco_wfsc_loop_erkin1(mp, out, cbx, dm1a, dm2a);

figure(1), print -dpng ../Figs/fig_psfx
%figure(103), print -dpng ../Figs/fig_time

end

return
%end % -----------------------------------------------
%

%figure(103), print -dpng ../Figs/fig_time

%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_05mm_run1_cbx_itr20 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_10mm_run1_cbx_itr21 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_15mm_run1_cbx_itr21 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_02jun2022/dat_15mm_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc1/dat_05mm_7lam_run1_cbx_itr15 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc1/dat_15mm_7lam_run1_cbx_itr15 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc1/dat_05mm_1lam_run1_cbx_itr20 dm1 dm2 tx cbx fom Im

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_7lam_run2_cbx_itr40 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_30mm_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_50mm_7lam_run5_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im 
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc2/dat_10mm_15per_7lam_run4_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value

    mp.dm1.V = dm1*0;
    mp.dm2.V = dm2*0;

    aa = mean(Im(mp.Fend.corr.maskBool))
    
    tic
    %InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring region of dark hole.
    %Im = falco_get_summed_image(mp);
    toc

%end

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)))

    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn1(:)))
    
    psfm = Im; % .* maskc1;

    cx = [-11 -8];
    cx1 = [-8 -3];

    %end

    tx = '7\lambda-20%: ';

    figure(2), clf
        imagesc(xv,xv,log10(aa),cx); axis xy image, colormap(jet); 
        %hold on, plot(xc, yc, 'w', xc1, yc1, 'w', xc2, yc2, 'w', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %caxis([-11 -6]);                
        %xlabel('Field-Angle [$\lambda_{central}$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title(sprintf('15mm, \\Delta\\lambda=100nm, 7\\lambdas: C(2.4 - 24) = %0.2g', cb),'Fontsize',fs,'Interpreter', 'Latex');
        %title(['post-10mm:' tx sprintf('Cs = %0.3g, Cb = %0.3g', cb1, cb)]);
        title([tx sprintf('Cs = %0.3g, Cb = %0.3g', cb1, cb)]);
   print -dpng ../Figs/fig_psf

   return

%{   
    Im = falco_get_summed_image(mp);
    psfn = Im .* maskc;
    cb = mean(nonzeros(psfn(:)))
    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn1(:)))
    
    figure(3), clf
        imagesc(xv,xv,log10(Im),cx1); axis xy image, colormap(jet); 
        hold on, plot(xc, yc, 'w', xc1, yc1, 'w', xc2, yc2, 'w', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title(['pre' tx sprintf('Cs = %0.3g, Cb = %0.3g', cb1, cb)]);
   print -dpng ../Figs/fig_psf0

   return

%}
   %end
    x = [0:length(cbx)-1];
    figure(3), clf, semilogy(x, cbx, 'ro-'), grid
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iteration', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean of Broadband NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Contrast vs WFC Iteration','Fontsize',fs,'Interpreter', 'Latex');
        legend(sprintf('Last Cb = %0.3g', cbx(end)), 'Location', 'Northeast'), drawnow
print -dtiff ../Figs/fig_cbx
    
    aa = dm1; rmsa = stat2d(aa); cx = rmsa(1)*3*[-1 1]; %[min(aa(:)) max(aa(:))];
    figure(4),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM1 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    print -dpng ../Figs/fig_dm1

    aa = dm2; rmsa = stat2d(aa); cx = rmsa(1)*3*[-1 1]; %[min(aa(:)) max(aa(:))];
    figure(5),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM2 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    print -dpng ../Figs/fig_dm2

return
% --------------------------------------------------------------------












if 0
%--Compute image
    tic; fprintf('Generating the PSF for a finite stellar size…') 
    mp.full.pol_conds = [0]; %--Which polarization states to use when creating an image.
    mp.full.TTrms = 0; % [mas]
    mp.full.Dstar = 1; % [mas]
    mp.full.Dtel  = 6; % [meters]
    mp.full.TipTiltNacross = 7; 
    Im = falco_get_summed_image_TipTiltPol(mp);
    fprintf('done. Time = %.2f s\n',toc)

else

    tic
    %InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring region of dark hole.
    Im = falco_get_summed_image(mp);
    toc

end

%end
    
    psfn = Im .* maskc;
    cb = mean(nonzeros(psfn(:)));
    
    psfm = Im; % .* maskc1;

    figure(1), clf
        imagesc(xv,xv,log10(psfm),[-11 -6]); axis xy image, colormap(jet); 
        hold on, plot(xc, yc, 'y', xc1, yc1, 'y', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %caxis([-11 -6]);                
        %xlabel('Field-Angle [$\lambda_{central}$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title(sprintf('Apodized 512pix Vortex NI: Cb = %0.2g', cb),'Fontsize',fs,'Interpreter', 'Latex');
        %title(['$\lambda$0 = 650nm, 10per-BB:' sprintf('Cb = %0.2g', cb)],'Fontsize',fs,'Interpreter', 'Latex');
        title(['Part-3: Dstar = 1mas, 10\%-BB:' sprintf('Cb = %0.2g', cb)],'Fontsize',fs,'Interpreter', 'Latex');
        %title(['dsn4-run1-itr49, BBTT, $\lambda$0 = 650nm, 10\%-BB: ' sprintf('Cb = %0.2g', cb)],'Fontsize',fs,'Interpreter', 'Latex');
   %print -dpng ../Figs/fig_psf

   return

if 0
    
    psfv = radial_average(Im);
    rv   = radial_average(rm);
    
    figure(2), clf, semilogy(rv, psfv, 'ro-'), grid
        xlabel('Radial Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Radial Average of NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Radial Average of 10\% Broadband NI','Fontsize',fs,'Interpreter', 'Latex');
    %print -dtiff ../Figs/fig_psfv
    
    x = [0:length(cbx)-1];
    figure(3), clf, semilogy(x, cbx, 'ro-'), grid
        parms = fun_newaxes(14, 2, 5);  
        xlabel('EFC Iteration', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean of Broadband NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Contrast vs WFC Iteration','Fontsize',fs,'Interpreter', 'Latex');
    %print -dtiff ../Figs/fig_cbx
    
    aa = dm1; rmsa = stat2d(aa); cx = rmsa(1)*3*[-1 1]; %[min(aa(:)) max(aa(:))];
    figure(4),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM1 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    %print -dpng ../Figs/fig_dm1

    aa = dm2; rmsa = stat2d(aa); cx = rmsa(1)*3*[-1 1]; %[min(aa(:)) max(aa(:))];
    figure(5),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM2 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    %print -dpng ../Figs/fig_dm2
end
    
rv0   = rv;    
cb0   = cb;
psf0  = Im;
psfv0 = psfv;
dm1_save = dm1;
dm2_save = dm2;

return
%end % -----------------------------------------------
%

cc = 1





bb = mp.P3.compact.mask + mp.P1.compact.mask;
aa = mp.P3.full.mask + mp.P1.full.mask;
%aa = mp.P3.full.mask;
psize = size(aa)

s1 = sprintf('ATSA Pupil + Lyot(%0.2f-%0.2f)', mp.P3.IDnorm, mp.P3.ODnorm);

    figure(10),clf, imagesc(aa), axis xy image, colormap(gray), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        %colorbar, parms = fun_newaxes(14, 0, 0);  
        title(s1, 'Interpreter','LaTex')
        %xlabel('$\lambda_{central}$/D','Fontsize',20,'Interpreter','LaTex');
    print -dpng ~esidick/Afalco/falco20200916/Figs/fig_img

    figure(11),clf, imagesc(bb), axis xy image, colormap(gray), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        %colorbar, parms = fun_newaxes(14, 0, 0);  
        title('Compact: ATSA Pupil + Lyot', 'Interpreter','LaTex')
        %xlabel('$\lambda_{central}$/D','Fontsize',20,'Interpreter','LaTex');
    %print -dpng ~esidick/Afalco/falco20200916/Figs/fig_img
return
%end % -----------------------------------------------
%

figure(2), clf, semilogy(rv0, psfv0, 'r-', 'Linewidth', 2), hold on, drawnow

Colx = 'rbgcmk';

clear psfx psfvx legx

xx = [0.02:0.02:0.1]; %gain error
xx = [0.02:0.02:0.1]; % piston error
rng('shuffle');

cbb = [];
psfvx = [];

legx{1} = sprintf('\\delta = 0, Cb = %0.2g', cb0);
%legx{1} = sprintf('\\Deltah = 0pm, Cb = %0.2g', cb0);

for jj = 1:length(xx)
    
    del0 = xx(jj);
    psfs = 0;
    clear psfm
        
    tic
    for ii = 1:10
        
        dm1 = (1 + del0 * randn(size(dm1_save))) .* dm1_save;
        dm2 = (1 + del0 * randn(size(dm2_save))) .* dm2_save;
        
        %dm1 = dm1_save + del0 * randn(64,64);
        %dm2 = dm2_save + del0 * randn(64,64);

        mp.dm1.V = dm1;
        mp.dm2.V = dm2;
    
        Im = falco_get_summed_image(mp);
        
        psfm(:,:,ii) = Im;
        psfx(:,:,ii,jj) = Im;
    end
    toc
    
    psf  = mean(psfm, 3);
    psfn = psf .* maskc;
    cb   = mean(nonzeros(psfn(:)));
    cbb  = [cbb cb];
    
    legx{1+jj} = sprintf('\\delta = %0.3f, Cb = %0.2g', xx(jj), cb);
    %legx{1+jj} = sprintf('\\Deltah = %ipm, Cb = %0.2g', round(xx(jj)*1e3), cb);
    
    psfvi = radial_average(psf);
    psfvx = [psfvx; psfvi];
    
    figure(1), clf
        imagesc(xv,xv,log10(Im), [-11 -5]); axis xy image, colormap(jet); 
        hold on, plot(xc, yc, 'b', xc1, yc1, 'b', 'Linewidth', 2), hold off
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title(['NI: $\delta$ = ' sprintf('%0.2f, Cb = %0.2g', del0, cb)],'Fontsize',fs,'Interpreter', 'Latex');
        %title(['NI: $\Delta$h = ' sprintf('%ipm, Cb = %0.2g', round(xx(jj)*1e3), cb)],'Fontsize',fs,'Interpreter', 'Latex');
   print('-dpng', sprintf('../Figs/fig_psf%i', jj))
   
    figure(2), semilogy(rv, psfvi, Colx(jj+1), 'Linewidth', 2), hold on, drawnow
    
end


 %end % main if ---------------------------------------

 dd = 1

 fs = 14;


    figure(2), grid %hold off, grid
        xlabel('Radial Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Radial Average of NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Effects of DM Gain-Errors on 10\% Broadband NI','Fontsize',fs,'Interpreter', 'Latex');
        legend(legx, 'Location', 'Southwest')
        axis([0 12 1e-11 1e-6])
    print -dtiff ../Figs/fig_psfv
   
 %save ../dat_erkin/part3_psfx_gain_error psfvx legx xx cbb psf0 cb0 xv  
 %save ../dat_erkin/part3_psfx_piston_error psfvx legx xx cbb psf0 cb0 xv  
    
return
%end % -----------------------------------------------
%





xx = [-0.05:0.01:0.05]; % lam/D

cbx = [];
mx  = [];
%psfvx = [];

figure(1), clf
figure(2), clf

%legx{1} = sprintf('\\Deltah = 0pm, Cb = %0.2g', cb0);

for jj = 1:length(xx)
    
    source_x_offset = xx(jj);

    %----------------------------
        Im = falco_get_summed_image(mp);
    %----------------------------
            
    psfn = psf .* maskc;
    cb   = mean(nonzeros(psfn(:)));
    cbx  = [cbx cb];
    mx   = [mx xx(jj)];

    psfc = psf .* maskc1;

    %legx{1+jj} = sprintf('\\Deltah = %ipm, Cb = %0.2g', round(xx(jj)*1e3), cb);
    %psfi  = pad(psf, npad);
    %psfvi = radial_average(psfi);
    %psfvx = [psfvx; psfvi];
    
    figure(1), subplot(3,4,jj)
        imagesc(xv,xv,log10(psfc), [-11 -5]); axis xy image, colormap(jet); niceaxes
        title(sprintf('$\\Delta$x = %0.2f', xx(jj)), 'Interpreter', 'Latex')
        axis(12*[-1 1 -1 1])
        drawnow

    figure(2), semilogy(mx, cbx, 'ro-', 'Linewidth', 2), grid, drawnow

end

fs = 14;

    figure(1), print -dpng ../Figs/fig_psf

%end

%cba = cbx;

load ../dat/dat_cb_vs_dx_nb xx cbx


    figure(2), semilogy(mx, cbx, 'ro-', mx, cba, 'bs-', 'Linewidth', 2), grid, drawnow
        xlabel('Source X-Tilt Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Cb vs Source X-Tilt: NB $\lambda$ = 550nm','Fontsize',fs,'Interpreter', 'Latex');
        legend('Flat-DM', 'With DM', 'Location', 'Southwest')
        %axis([0 14 5e-11 1e-6])
    print -dtiff ../Figs/fig_cb

return
%end % -----------------------------------------------
%


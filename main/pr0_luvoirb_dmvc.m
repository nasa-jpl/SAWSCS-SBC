% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

%clear all;



if 1

%% Step 1: Define Necessary Paths on Your Computer System

%--Required packages are FALCO and PROPER. 
% Add FALCO to the MATLAB path with the command:  addpath(genpath(full_path_to_falco)); savepath;
% Add PROPER to the MATLAB path with the command:  addpath(full_path_to_proper); savepath;

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;


%%--Output Data Directories (Comment these lines out to use defaults within falco-matlab/data/ directory.)
% mp.path.config = ; %--Location of config files and minimal output files. Default is [mp.path.falco filesep 'data' filesep 'brief' filesep]
% mp.path.ws = ; % (Mostly) complete workspace from end of trial. Default is [mp.path.falco filesep 'data' filesep 'ws' filesep];
% mp.flagSaveWS = false;  %--Whether to save out entire (large) workspace at the end of trial. Default is false


%% Step 2: Load default model parameters

EXAMPLE_defaults_LUVOIRB_VC_design_erkin


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

if 1  % Not overall control !!!
    
mp.estimator = 'perfect';
mp.full.flagGenPupil = false;
%mp.false.flagGenPupil = false;
mp.compact.flagGenPupil = false;

%pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_cirle512.fits');
pupil = fitsread('/home/esidick/Atsa/Jeff/map/pupil_atsa512.fits');

mp.P1.compact.mask = pupil; %(load your file here)
mp.P1.full.mask    = pupil; %(load your file here)
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

%return

%end

%% Wavefront Control: Controller Specific
% Controller options: 
%  - 'gridsearchEFC' for EFC as an empirical grid search over tuning parameters
%  - 'plannedEFC' for EFC with an automated regularization schedule
%  - 'SM-CVX' for constrained EFC using CVX. --> DEVELOPMENT ONLY
mp.controller = 'gridsearchEFC';
%mp.controller = 'plannedefc';

% % % GRID SEARCH EFC DEFAULTS     
%--WFSC Iterations and Control Matrix Relinearization
mp.Nitr = 100; %--Number of estimation+control iterations to perform
mp.relinItrVec = 1:mp.Nitr;  %--Which correction iterations at which to re-compute the control Jacobian
mp.dm_ind = [1 2]; %--Which DMs to use

if 0         
       SetA2 = [1, 1j, 12, 1, 1];  %--DMs 1 & 2. Relinearize every iteration.
        
        SetB2 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            5, -3, 12, 0, 0;...
            3, -5, 12, 0, 0;...
            5, -2, 12, 0, 0;...
            2, -1, 12, 0, 0;...
            ];
        
        SetC2 = ... %--DMs 1, 2, & 9. At first iteration only, relinearize and compute the new optimal Beta.
            [0, 0, 0, 1, 0;...
            10, -2, 1289, 0, 0;...
            ];
        
       mp.ctrl.sched_mat = [...
           repmat(SetB2,[2,1]);...
           %repmat(SetC2,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

end

%end

%load /home/esidick/Afalco/falco20200916/dat_erkin/run1_cbx_itr10 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run2_cbx_itr25 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run3_cbx_itr45 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run6_cbx_itr125 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run7_cbx_itr150 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run8_cbx_itr180 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run9_cbx_itr220 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run10_cbx_itr225 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run11_cbx_itr265 dm1 dm2 tx cbx 
%load /home/esidick/Afalco/falco20200916/dat_erkin/run12_cbx_itr295 dm1 dm2 tx cbx 
load /home/esidick/Afalco/falco20200916/dat_erkin/run13_cbx_itr325 dm1 dm2 tx cbx 
    mp.dm1.V = dm1;
    mp.dm2.V = dm2;
    dm1a = dm1;
    dm2a = dm2;

%[mp, out] = falco_wfsc_loop_erkin(mp, out, cbx, dm1a, dm2a);

%return
%end % -----------------------------------------------
%

tic
InormHist = zeros(mp.Nitr,1); % Measured, mean raw contrast in scoring region of dark hole.
Im = falco_get_summed_image(mp);
toc

%end

close all

fs = 16;

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 10);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1] = fun_circle(xv, 10);
    
    psfn = Im .* maskc;
    cb = mean(nonzeros(psfn(:)));
    
    cx = [-11 -5];
    
    figure(1), clf
        imagesc(xv,xv,log10(Im), cx); axis xy image, colormap(jet); 
        hold on, plot(xc, yc, 'r', xc1, yc1, 'r', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %caxis([-11 -6]);                
        %xlabel('Field-Angle [$\lambda_{central}$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title(sprintf('Apodized 512pix Vortex NI: Cb = %0.2g', cb),'Fontsize',fs,'Interpreter', 'Latex');
        title(['559nm-NI: $\delta$ = 0,' sprintf('Cb = %0.2g', cb)],'Fontsize',fs,'Interpreter', 'Latex');
        axis(11.5*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psf1
    
    psfv = radial_average(Im);
    rv = radial_average(rm);
    
    figure(2), clf, semilogy(rv, psfv, 'ro-'), grid
        xlabel('Radial Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Radial Average of NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('550nm: Radial Average of 10\% Broadband NI','Fontsize',fs,'Interpreter', 'Latex');
    %print -dtiff ../Figs/fig_psfv1
    
return    
    
if 0    
    x = [0:length(cbx)-1];
    figure(3), clf, semilogy(x, cbx, 'ro-'), grid
        xlabel('EFC Iteration', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean of Broadband NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Contrast vs WFC Iteration','Fontsize',fs,'Interpreter', 'Latex');
    print -dtiff ../Figs/fig_cbx
    
    aa = dm1; rmsa = stat2d(aa); cx = 3*rmsa(1)*[-1 1]; %[min(aa(:)) max(aa(:))];
    figure(4),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM1 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    print -dpng ../Figs/fig_dm1

    aa = dm2; rmsa = stat2d(aa); %cx = [min(aa(:)) max(aa(:))];
    figure(5),clf, imagesc(aa, cx), axis xy image, colormap(jet), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM2 Piston Values: RMS = %0.1f, PV = %0.1fnm', rmsa), 'Interpreter','LaTex')
    print -dpng ../Figs/fig_dm2
end

rv0   = rv;    
cb0   = cb;
psf0  = Im;
psfv0 = psfv;
dm1_save = dm1;
dm2_save = dm2;
    
%return
%end % -----------------------------------------------
%

figure(2), clf, semilogy(rv, psfv, 'r-', 'Linewidth', 2), hold on, drawnow

Colx = 'rbgcmk';

clear psfx psfvx legx

xx = [0.002:0.002:0.01];
rng('shuffle');

cbb = [];
psfvx = [];

legx{1} = sprintf('\\delta = 0, Cb = %0.2g', cb0);

for jj = 1:length(xx)
    
    del0 = xx(jj);
    psfs = 0;
    clear psfm
        
    for ii = 1:10
        dm1 = (1 + del0 * randn(size(dm1_save))) .* dm1_save;
        dm2 = (1 + del0 * randn(size(dm2_save))) .* dm2_save;

        mp.dm1.V = dm1;
        mp.dm2.V = dm2;
    
        tic
        Im = falco_get_summed_image(mp);
        toc
        
        psfm(:,:,ii) = Im;
        psfx(:,:,ii,jj) = Im;
    end
    
    psf  = mean(psfm, 3);
    psfn = psf .* maskc;
    cb   = mean(nonzeros(psfn(:)));
    cbb  = [cbb cb];
    
    legx{1+jj} = sprintf('\\delta = %0.3f, Cb = %0.2g', xx(jj), cb);
    
    psfvi = radial_average(psf);
    psfvx = [psfvx; psfvi];
    
    figure(1), clf
        imagesc(xv,xv,log10(Im), [-11 -5]); axis xy image, colormap(jet); colorbar
        hold on, plot(xc, yc, 'r', xc1, yc1, 'r', 'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title(['NI: $\delta$ = ' sprintf('%0.3f, Cb = %0.2g', del0, cb)],'Fontsize',fs,'Interpreter', 'Latex');
   print('-dpng', sprintf('../Figs/fig_psf%i', jj))
   
    figure(2), semilogy(rv0, psfvi, Colx(jj+1), 'Linewidth', 2), hold on, drawnow
    
end

%end

    figure(2), hold off, grid
        xlabel('Radial Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Radial Average of NI', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Radial Average of 10\% Broadband NI','Fontsize',fs,'Interpreter', 'Latex');
        legend(legx, 'Location', 'Northeast')
        axis([0 23 1e-11 1e-6])
    print -dtiff ../Figs/fig_psfv
   
 %save ../dat_erkin/dat_no_apod512_psfx psfvx legx xx cbb cb0 psf0 xv  
    
return
end % -----------------------------------------------
%
dm0 = randn(64,64);
yy = dm0(:);
yy = yy';
xx = [-3.5:0.1:3.5];
y1 = hist(yy, xx);

    figure(2), plot(xx, y1, 'ro-'), grid
        xlabel('Values of randn(64,64)', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Histogram', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('Distribution of randn(64,64): STD = %0.2f',std(yy)), 'Fontsize',fs,'Interpreter', 'Latex');
        %legend(legx, 'Location', 'Northeast')
        %axis([0 23 1e-11 1e-6])
    print -dtiff ../Figs/fig_std



return
%end % -----------------------------------------------
%





bb = mp.P3.compact.mask + mp.P1.compact.mask;
aa = mp.P4.full.mask + mp.P1.full.mask;
%aa = mp.P3.full.mask;
psize = size(aa)

s1 = sprintf('ATSA Pupil + Apod(%0.2f-%0.2f)', mp.P3.IDnorm, mp.P3.ODnorm);
s1 = sprintf('ATSA Pupil + Lyot(%0.2f-%0.2f)', mp.P4.IDnorm, mp.P4.ODnorm);

    figure(1),clf, imagesc(aa), axis xy image, colormap(gray), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        %colorbar, parms = fun_newaxes(14, 0, 0);  
        title(s1, 'Interpreter','LaTex')
        %xlabel('$\lambda_{central}$/D','Fontsize',20,'Interpreter','LaTex');
    print -dpng ~esidick/Afalco/falco20200916/Figs/fig_img1

    figure(2),clf, imagesc(bb), axis xy image, colormap(gray), %niceaxes, %colorbar, title('Frame-1')
        parms = fun_newaxes(14, 0, 0);  
        %colorbar, parms = fun_newaxes(14, 0, 0);  
        title('Compact: ATSA Pupil + Lyot', 'Interpreter','LaTex')
        %xlabel('$\lambda_{central}$/D','Fontsize',20,'Interpreter','LaTex');
    %print -dpng ~esidick/Afalco/falco20200916/Figs/fig_img
return
%end % -----------------------------------------------
%



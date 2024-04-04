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

idd = 2
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

       mp.ctrl.sched_mat = [...
           repmat(SetB,[1,1]);...  % [2, 1] repeat 2 times
           %repmat(SetB,[1,1]);...
           %repmat(SetB3,[1,1]);...
           %repmat(SetB2,[1,1]);...
           %repmat(SetB4,[1,1]);...
           ];
        
[mp.Nitr, mp.relinItrVec, mp.gridSearchItrVec, mp.ctrl.log10regSchedIn, mp.dm_ind_sched] = falco_ctrl_EFC_schedule_generator(mp.ctrl.sched_mat);

%return    
%end % ------------------------------------------
%%

load /home/esidick/Afalco/falco20200916/dat_bb_256pix_part2/gray_gap_08mm_bw20_30nm_2to12lamD_chg6_4lam_itr410 dm1 dm2 tx cbx fom Im bvalo 

tc = fom(end)*10;

fs = 14;

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    maskc1 = (rm >= 2.5) & (rm <= 3.5);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, 3);
    [xc2, yc2] = fun_circle(xv, 12);


flag = 7
mp.wfc.model_id = flag; % =1 stantard, 2 = dm1_opd, 3 = dm2_opd/amp, =4 dm1_wfc, 5 = dm2_wfc, 6 = dm1_dm2_wfc, 7 = target_opd_amp
    %disp('paused'), pause
% ------------------------------------------

    id = 1;
    mp.dm1.V = dm1 * id;
    mp.dm2.V = dm2 * id;

    ri = 32;
    ci = 32;

    load maps_psd/paper_opd_amp_nom opd_target amp_target ce_target Im
    %Im = falco_get_summed_image(mp);
    %save maps_psd/paper_opd_amp_nom opd_target amp_target ce_target Im
    ampnom = amp_target;
    opdnom = opd_target;
    cenom  = ce_target;
    psf0 = Im .* maskc;
    cb0 = mean(nonzeros(psf0(:)))


    dval = 0.5; % [nm]

    mp.dm1.V(ri,ci) = mp.dm1.V(ri,ci) + dval; % [nm]
    Im = falco_get_summed_image(mp);
    load maps_psd/paper_opd_amp opd_target amp_target ce_target

    pupil = fitsread('/home/esidick/2022habex/dat_macos/pup256_macos.fits');
    maskopd = abs(pupil);

    ce = ce_target ./ (cenom + eps);
    opd_target = maskopd .* (1e9 * atan2(imag(ce),real(ce)) * 500e-9 / 2 / pi);

    amp = amp_target - ampnom;
    opd = opd_target;

aa = amp; rms0 = stat2d(aa); cx = rms0(1)*5*[-1 1];


figure(1), clf, imagesc(aa, cx);  axis xy image; colormap(jet), niceaxes
    hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('amp: poke act(%i,%i) = 0.5nm', ri, ci)), 
    xlabel(sprintf('rms = %0.3g, pv = %0.1f', rms0)), 
%print -dpng ../Figs/fig_amp1

aa = opd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

figure(2), clf, imagesc(aa);  axis xy image; colormap(jet), niceaxes
    hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
    parms = fun_newaxes(14, 0, 0);  
    colorbar, parms = fun_newaxes(14, 0, 0);  
    title(sprintf('opd: poke act(%i,%i) = 0.5nm', ri, ci)), 
    xlabel(sprintf('rms = %0.1f, pv = %0.1fnm', rms0)), 
%print -dpng ../Figs/fig_opd1

    psfn = Im .* maskc; 
    cb1 = mean(nonzeros(psfn(:)));

    fs = 14;
    cx = [-11 -8];

    t1 = ['gap08, bw20, rms30, opdhi, 4\lambda'];
    t2 = [sprintf('poke act(%i,%i) = 0.5nm: Cb = %0.3g',  ri,ci, cb1)];
    t1 = [sprintf('nominal: Cb = %0.3g',  cb0)];

    figure(3), clf
        imagesc(xv,xv,log10(psfn), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title(['post-30nm-20%-12\lambda, ' sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%']);
        title(t2);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfo

    figure(4), clf
        imagesc(xv,xv,log10(psf0), cx); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        %xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        %title([sprintf('pre-efc: 1\\lambda, [Cs Cb] = [%0.3g %0.3g]', cb1, cb)]);
        %title(['post-30nm-20%-12\lambda, ' sprintf('Cb = %0.3g, Tc = %0.2f',  cb0, tc) '%']);
        title(t1);
        axis(12*[-1 1 -1 1])
   %print -dpng ../Figs/fig_psfo

%return
end % -----------------------------------------------
%%
idd = 2
%disp('paused'), pause

    load maps_psd/paper_opd_amp_nom opd_target amp_target ce_target Im
    %Im = falco_get_summed_image(mp);
    %save maps_psd/paper_opd_amp_nom opd_target amp_target ce_target Im
    ampnom = amp_target;
    opdnom = opd_target;
    cenom  = ce_target;
    psf0 = Im .* maskc;
    cb0 = mean(nonzeros(psf0(:)))

dval = 1; % [nm]

s1 = [30:60:330];
s2 = [15:30:345];

rv = [32];
cv = [32];

for ii = 1:length(s1)
    ri = round(32 + 13 * cosd(s1(ii)));
    ci = round(32 + 13 * sind(s1(ii)));
    rv = [rv ri];
    cv = [cv ci];
end

for ii = 1:length(s2)
    ri = round(32 + 24 * cosd(s2(ii)));
    ci = round(32 + 24 * sind(s2(ii)));
    rv = [rv ri];
    cv = [cv ci];
end

ns = length(rv)

if 0

figure(2), clf, plot(rv, cv, 'ro'), axis square, grid
    parms = fun_newaxes(14, 0, 5);  
    xlabel('act number')
    ylabel('act number')
    xline(32, 'g--', 'Linewidth',2);
    yline(32, 'g--','Linewidth',2);
    axis([1 64 1 64])
   %print -dpng ../Figs/fig_act
%return
end

video = VideoWriter('dat_img/erkin_movie_dm2.avi'); %create the video object
%video = VideoWriter('dat_img/my_movie_dm', 'mpeg-4'); %create the video object
video.FrameRate = 1;
open(video); %open the file for writing

for ii = 1:ns

    ri = rv(ii);
    ci = cv(ii);

    mp.dm1.V = dm1;
    mp.dm2.V = dm2;

    %mp.dm1.V(ri,ci) = mp.dm1.V(ri,ci) + dval; % [nm]
    mp.dm2.V(ri,ci) = mp.dm2.V(ri,ci) + dval; % [nm]
    Im = falco_get_summed_image(mp);
    psfn = Im .* maskc; 
    cb1 = mean(nonzeros(psfn(:)));

    %load maps_psd/paper_opd_amp_nom opd_target amp_target ce_target 
    %save maps_psd/paper_opd_amp_nom opd_target amp_target ce_target Im
    
    load maps_psd/paper_opd_amp opd_target amp_target ce_target

    ce = ce_target ./ (cenom + eps);
    opd_target = maskopd .* (1e9 * atan2(imag(ce),real(ce)) * 500e-9 / 2 / pi);

    amp = amp_target - ampnom;
    opd = opd_target;

    sub_plot_spie

        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(video,frame);
        
end
      
close(video); %close the file

return
%end % -----------------------------------------------
%%


close all

fig_size = [10 10 1200 1200];
handles.master = figure('Color',[1 1 1],'Position',fig_size);
set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')

    amp = amp_target - ampnom;
    aa = amp; rms0 = stat2d(aa); cx = rms0(1)*5*[-1 1];

            fst = 14; %--Font size for titles in the subplots

            h_psf = subplot(2,2,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.5, [1 1]*0.45])
            imagesc(aa, cx); niceaxes
            hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            xlabel(sprintf('rms = %0.3g, pv = %0.3g', rms0),'FontSize',14);
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('amp: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

aa = opd_target; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];

            h_psf = subplot(2,2,3); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.02, [1 1]*0.45])
            imagesc(aa); niceaxes
            hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('phase: rms = %0.2f, pv = %0.2fnm', rms0),'FontSize',fst);
            %title(sprintf('phase: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

            h_psf = subplot(2,2,2); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.5, [1 1]*0.45])
            imagesc(xv,xv,log10(psf0), [-11 -8]);  
            axis xy equal tight; colorbar(h_psf); colormap(jet);
            axis(12*[-1 1 -1 1])
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('nominal NI: Cb = %0.3g', cb0))
            
            h_psf = subplot(2,2,4); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.02, [1 1]*0.45])
            imagesc(xv,xv,log10(psfn), [-11 -8]);  
            axis xy equal tight; colorbar(h_psf); colormap(jet);
            axis(12*[-1 1 -1 1])
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('distorted NI: Cb = %0.3g', cb1))
return
%end % -----------------------------------------------
%%







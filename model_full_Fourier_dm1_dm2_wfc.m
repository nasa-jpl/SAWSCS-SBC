% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function Eout = model_full_Fourier(mp, lambda, Ein, normFac)
% 
% Function to run the full-knowledge optical model and return the final
% E-field.
% - Not used by the estimator and controller.
% - Only used to create simulated intensity images.
%
% "Fourier" layout:  
% - For design or basic modeling only. 
% - Instead use a full model with Fresnel propagations for more realistic simulations.
%
% ---------------
% INPUTS:
% - mp = structure of model parameters
% - lambda = wavelength in meters
% - Ein = 2D input E-field at entrance
% - normFac = intensity normalization factor 
%
%
% OUTPUTS:
% - Eout = 2-D electric field at final plane of optical layout
% - varargout{1}==Efiber = E-field at final plane when a single mode fiber
% is used
%
% REVISION HISTORY:
% --------------
% Modified on 2019-04-18 by A.J. Riggs to use varargout for Efiber instead
% of having Efiber as a required output. 
% Modified on 2019-04-05 by A.J. Riggs to have the normalization be
%   computed by moving the source off-axis instead of removing the FPM.
% Modified on 2019-02-14 by G. Ruane to handle scalar vortex FPMs
% Modified on 2019-02-14 by A.J. Riggs to be the "Fourier" layout for all
% types of coronagraphs.
% Modified on 2017-10-17 by A.J. Riggs to have model_full.m be a wrapper. All the 
%  actual full models, including this one, have been moved to sub-routines for clarity.
% Modified by A.J. Riggs from hcil_simTestbed.m to model_full.m.
% Modified on 2015-02-18 by A.J. Riggs from hcil_model.m to hcil_simTestbed.m to include
%  extra errors in the model to simulate the actual testbed for fake images.

function [Eout, varargout] = model_full_Fourier_dm1_dm2_wfc(mp, lambda, Ein, normFac)

mirrorFac = 2; % Phase change is twice the DM surface height.
NdmPad = mp.full.NdmPad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Change model values if the full model has a different value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(mp,'full'))
    if(isfield(mp.full,'dm1'))
        if(isfield(mp.full.dm1,'xc'));  mp.dm1.xc = mp.full.dm1.xc;  end % x-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1,'yc'));  mp.dm1.yc = mp.full.dm1.yc;  end % y-center location of DM1 surface [actuator widths]
        if(isfield(mp.full.dm1,'V0'));  mp.dm1.V = mp.dm1.V + mp.full.dm1.V0;  end % Add some extra starting command to the voltages  [volts]
    end
    if(isfield(mp.full,'dm2'))
        if(isfield(mp.full.dm2,'xc'));  mp.dm2.xc = mp.full.dm2.xc;  end % x-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2,'yc'));  mp.dm2.yc = mp.full.dm2.yc;  end % y-center location of DM2 surface [actuator widths]
        if(isfield(mp.full.dm2,'V0'));  mp.dm2.V = mp.dm2.V + mp.full.dm2.V0;  end % Add some extra starting command to the voltages  [volts]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masks and DM surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np   = 512/2;
indx = mp.wfc.indx;
nwfc = 20

%flag = 5
%disp('dm1/dm2-wfc paused'), pause

ma = mp.dm1.Nact;
ma2 = ma * ma;
du0 = zeros(ma2,1);

maskopd = abs(pad(mp.P1.full.mask,np));

ndm1 = length(mp.wfc.km1);
ndm2 = length(mp.wfc.km2);
ndm  = ndm1 + ndm2;

t1 = clock;
tsum = 0;

clear dm1x dm2x opdx ampx
rms_opd = [];
rms_amp = [];
dmrms = [];
dmsum = 0;

dm1 = mp.dm1.V;
dm2 = mp.dm2.V;

ce_target = mp.wfc.ce_target;

% ---------------------------------------

for jj = 1:nwfc+1

if(any(mp.dm_ind==1));  DM1surf = falco_gen_dm_surf(mp.dm1,mp.dm1.dx,NdmPad); end
if(any(mp.dm_ind==2));  DM2surf = falco_gen_dm_surf(mp.dm2,mp.dm2.dx,NdmPad); else; DM2surf=zeros(NdmPad); end

pupil = padOrCropEven(mp.P1.full.mask,NdmPad);
Ein = padOrCropEven(Ein,NdmPad);

if(mp.useGPU)
    pupil = gpuArray(pupil);
    Ein = gpuArray(Ein);
    if(any(mp.dm_ind==1)); DM1surf = gpuArray(DM1surf); end
    if(any(mp.dm_ind==2)); DM2surf = gpuArray(DM2surf); end
end

if(mp.flagDM1stop); DM1stop = padOrCropEven(mp.dm1.full.mask, NdmPad); else; DM1stop = 1; end
if(mp.flagDM2stop); DM2stop = padOrCropEven(mp.dm2.full.mask, NdmPad); else; DM2stop = 1; end

if(mp.flagDMwfe)
    if(any(mp.dm_ind==1));  Edm1WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm1.wfe,NdmPad,'extrapval',0)); else; Edm1WFE = ones(NdmPad); end
    if(any(mp.dm_ind==2));  Edm2WFE = exp(2*pi*1i/lambda.*padOrCropEven(mp.dm2.wfe,NdmPad,'extrapval',0)); else; Edm2WFE = ones(NdmPad); end
else
    Edm1WFE = ones(NdmPad);
    Edm2WFE = ones(NdmPad);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagation: entrance pupil, 2 DMs, (optional) apodizer, vortex FPM, LS, and final focal plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Define pupil P1 and Propagate to pupil P2
EP1 = pupil.*Ein; %--E-field at pupil plane P1
EP2 = propcustom_relay(EP1,mp.Nrelay1to2,mp.centering); %--Forward propagate to the next pupil plane (P2) by rotating 180 degrees mp.Nrelay1to2 times.

%--Propagate from P2 to DM1, and apply DM1 surface and aperture stop
if( abs(mp.d_P2_dm1)~=0 ); Edm1 = propcustom_PTP(EP2,mp.P2.full.dx*NdmPad,lambda,mp.d_P2_dm1); else; Edm1 = EP2; end  %--E-field arriving at DM1
Edm1 = Edm1WFE.*DM1stop.*exp(mirrorFac*2*pi*1i*DM1surf/lambda).*Edm1; %--E-field leaving DM1

%--Propagate from DM1 to DM2, and apply DM2 surface and aperture stop
Edm2 = propcustom_PTP(Edm1,mp.P2.full.dx*NdmPad,lambda,mp.d_dm1_dm2); 
Edm2 = Edm2WFE.*DM2stop.*exp(mirrorFac*2*pi*1i*DM2surf/lambda).*Edm2;

%--Back-propagate to pupil P2
if( mp.d_P2_dm1 + mp.d_dm1_dm2 == 0 ); EP2eff = Edm2; else; EP2eff = propcustom_PTP(Edm2,mp.P2.full.dx*NdmPad,lambda,-1*(mp.d_dm1_dm2 + mp.d_P2_dm1)); end %--Back propagate to pupil P2

%--Re-image to pupil P3
EP3 = propcustom_relay(EP2eff,mp.Nrelay2to3,mp.centering);

    ce = pad(EP3,np);

if 0
    txt = ' target'
    ce_target = ce;
    pha = atan2(imag(ce),real(ce)) * lambda / 2 / pi;
    opd_target = pha .* maskopd; %mp.wfc.maskopd;
    amp_target = (abs(ce) ) .* maskopd;
    %save ~esidick/Afalco/falco20200916/macos/IFopd/target_opd_amp_gap08_bb_30nm_bw20_1500nm_post_efc_no_ota opd_target amp_target ce_target
end

    cer = ce ./ (ce_target + eps);
    pha_nm = 1e9 * atan2(imag(cer),real(cer)) * lambda / 2 / pi;
    %pha_nm = real(pha_nm);
    opd = pha_nm .* maskopd; %mp.wfc.maskopd;
    amp = (abs(ce) - mp.wfc.amp_target) .* maskopd;

    opdx(:,:,jj) = opd; % pre-wfc
    ampx(:,:,jj) = amp; % pre-wfc
    dm1x(:,:,jj) = mp.dm1.V; 
    dm2x(:,:,jj) = mp.dm2.V;

    rmsi = stat2d(opd)

    rms_opd = [rms_opd; rmsi];
    rms_amp = [rms_amp; stat2d(amp)];
    dmrms = [dmrms; stat2d(dmsum)];
    
    if 0
        m = 256; 
        bb = opd(:);
        k1 = find(bb>0); ap = bb(k1); meanp = mean(ap)*100;
        k2 = find(bb<0); am = bb(k2); meanm = mean(am)*100;
        k1 = find(bb>meanp); bb(k1) = 0; 
        k2 = find(bb<meanm); bb(k2) = 0; 
        opdf = reshape(bb,m,m);
        maski = abs(opdf) > 0;

        opd = opd .* maski;
        amp = amp .* maski;
    end

    %opd1= (pha_nm  - mp.wfc.opd_target*0) .* maskopd; %mp.wfc.maskopd;
    
    %disp('paused'), pause(30)

    dw = m2v(opd, indx);
    da = m2v(amp, indx);
    dw = [dw(:); da(:)];
    ddu =  -mp.wfc.G * dw;

    du1 = du0;
    du1(mp.wfc.km1) = ddu(1:ndm1);

    du2 = du0;
    du2(mp.wfc.km2) = ddu(ndm1+1:ndm);

    del_dm1 = reshape(du1, ma, ma);
    del_dm1 = del_dm1';

    del_dm2 = reshape(du2, ma, ma);
    del_dm2 = del_dm2';

    dmsum = dmsum + [del_dm1 del_dm2];

    dm1 = dm1 + del_dm1;
    mp.dm1.V = dm1;

    dm2 = dm2 + del_dm2;
    mp.dm2.V = dm2;

    %figure(10), clf, plot(rmsx(:,1), 'ro-'), grid, drawnow

    delt = etime(clock, t1);
    t1   = clock;
    tsum = tsum + delt;
    count = [jj nwfc+1 round(delt) round(tsum)]

end % nwfc - loop

    opdi = opdx(:,:,1);
    ampi = ampx(:,:,1);
    %opdx = opd; % post-wfc
    ampx = amp; % post-wfc
    dm1x = mp.dm1.V; 
    dm2x = mp.dm2.V;

    %save('~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dmx_opdx_08nm dmx dm1x dm2x opdx opd_in
%save(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'ampx', 'rmsx', 'dmrms') %pre-7/26/2023
    save(mp.wfc.fn, 'dm1x', 'dm2x', 'opdi', 'ampi', 'opdx', 'ampx', 'rms_opd', 'rms_amp','dmrms')

    %save ~esidick/Afalco/falco20200916/macos/dat_dm_err/test_opd opd
    %figure(1), clf, imagesc(pha_nm), axis xy image, colormap(jet), niceaxes,  colorbar, drawnow, print -dpng ../Figs/fig_opd1
    %figure(2), clf, imagesc(mp.wfc.opd_target), axis xy image, colormap(jet), niceaxes,  colorbar, drawnow, print -dpng ../Figs/fig_opd2
    
    
% ======================================================================

%save /home/esidick/2023_6mst/macos/dat/falco_efield_ep3 EP3
%save ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_field_no_dm EP3
%save ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_field_with_dm_v2 EP3
%save ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_field_no_error EP3
% ======================================================================

%--Apply the apodizer mask (if there is one)
if(mp.flagApod)
    EP3 = mp.P3.full.mask.*padOrCropEven(EP3, mp.P3.full.Narr); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Select propagation based on coronagraph type   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
switch upper(mp.coro)
    case{'VORTEX','VC','AVC'}
        if(mp.flagApod==false)
            EP3 = padOrCropEven(EP3,2^nextpow2(mp.P1.full.Narr)); %--Crop down if there isn't an apodizer mask
        end
        % Get FPM charge 
        if(numel(mp.F3.VortexCharge)==1)
            % single value indicates fully achromatic mask
            charge = mp.F3.VortexCharge;
        else
            % Passing an array for mp.F3.VortexCharge with
            % corresponding wavelengths mp.F3.VortexCharge_lambdas
            % represents a chromatic vortex FPM
            charge = interp1(mp.F3.VortexCharge_lambdas,mp.F3.VortexCharge,lambda,'linear','extrap');
        end
        EP4 = propcustom_mft_Pup2Vortex2Pup( EP3, charge, mp.P1.full.Nbeam/2, 0.3, 5, mp.useGPU );  %--MFTs
        EP4 = padOrCropEven(EP4,mp.P4.full.Narr);
        
        %{
        a = charge
        b = mp.P1.full.Nbeam
        c = mp.useGPU
        d = mp.P4.full.Narr
        ps = size(EP3)
        ps1 = size(EP4)
        %}
        %save ~esidick/Afalco/falco20200916/macos/dat/dat_ep3_field EP3
        %save ~esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_ep3_field_psd_lamo100_7lam EP3
        %disp('paused inside model_full_Fourier'), %pause(1)
        

    case{'SPLC','FLC'}
        %--MFT from SP to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering); %--E-field incident upon the FPM
        EF3 = mp.F3.full.mask.amp.*EF3inc; %--Apply FPM
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); %--E-field incident upon the Lyot stop

    case{'LC','APLC','RODDIER'}
        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        if(strcmpi(mp.coro,'roddier'))
            FPM = mp.F3.full.mask.amp.*exp(1i*2*pi/lambda*(mp.F3.n(lambda)-1)*mp.F3.t.*mp.F3.full.mask.phzSupport);
            EF3 = (1-FPM).*EF3inc; %--Apply (1-FPM) for Babinet's principle later
        else
            EF3 = (1-mp.F3.full.mask.amp).*EF3inc;
        end
        % Use Babinet's principle at the Lyot plane. This is the term without the FPM.
        EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4subtrahend = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
        %--Babinet's principle at P4
        EP4 = padOrCropEven(EP4noFPM,mp.P4.full.Narr) - EP4subtrahend;

    case{'HLC'}
        %--Complex transmission of the points outside the FPM (just fused silica with optional dielectric and no metal).
        t_Ti_base = 0;
        t_Ni_vec = 0;
        t_PMGI_vec = 1e-9*mp.t_diel_bias_nm; % [meters]
        pol = 2;
        [tCoef, ~] = falco_thin_film_material_def(lambda, mp.aoi, t_Ti_base, t_Ni_vec, t_PMGI_vec, lambda*mp.FPM.d0fac, pol);
        transOuterFPM = tCoef;

        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        EF3 = (transOuterFPM-mp.FPM.mask).*EF3inc; %- transOuterFPM instead of 1 because of the complex transmission of the glass as well as the arbitrary phase shift.
        % Use Babinet's principle at the Lyot plane.
        EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
        EP4noFPM = transOuterFPM*padOrCropEven(EP4noFPM,mp.P4.full.Narr); %--Apply the phase and amplitude change from the FPM's outer complex transmission.
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
        %--Babinet's principle at P4
        EP4 = EP4noFPM-EP4subtra;

    case{'EHLC'}
        
        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
        EF3 = mp.FPM.mask.*EF3inc; %--Apply FPM
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4 = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     

    case{'FOHLC'}
        %--FPM representation (idealized as amplitude and phase)
        DM8amp = falco_gen_HLC_FPM_amplitude_from_cube(mp.dm8,'full');
        DM8ampPad = padOrCropEven( DM8amp,mp.full.Nfpm,'extrapval',1);
        DM9surf = falco_gen_HLC_FPM_surf_from_cube(mp.dm9,'full');
        DM9surfPad = padOrCropEven( DM9surf,mp.full.Nfpm);
        FPM = DM8ampPad.*exp(2*pi*1i/lambda*DM9surfPad);

        %--MFT from apodizer plane to FPM (i.e., P3 to F3)
        EF3inc = propcustom_mft_PtoF(EP3, mp.fl,lambda,mp.P2.full.dx,mp.F3.full.dxi,mp.F3.full.Nxi,mp.F3.full.deta,mp.F3.full.Neta,mp.centering);
        % Apply (1-FPM) for Babinet's principle later
        EF3 = (1 - FPM).*EF3inc;
        % Use Babinet's principle at the Lyot plane.
        EP4noFPM = propcustom_relay(EP3,mp.Nrelay3to4,mp.centering); %--Propagate forward another pupil plane 
        EP4noFPM = padOrCropEven(EP4noFPM,mp.P4.full.Narr);
        %--MFT from FPM to Lyot Plane (i.e., F3 to P4)
        EP4subtra = propcustom_mft_FtoP(EF3,mp.fl,lambda,mp.F3.full.dxi,mp.F3.full.deta,mp.P4.full.dx,mp.P4.full.Narr,mp.centering); % Subtrahend term for Babinet's principle     
        %--Babinet's principle at P4
        EP4 = EP4noFPM - EP4subtra;

    otherwise
        error('model_full_Fourier.m: Modely type\t %s\t not recognized.\n',mp.coro);
end

%--Remove the FPM completely if normalization value is being found
if(normFac==0)
    switch upper(mp.coro)
        case{'VORTEX','VC','AVC'}
            EP4 = propcustom_relay(EP3,mp.Nrelay3to4, mp.centering);
            EP4 = padOrCropEven(EP4,mp.P4.full.Narr);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Back to common propagation any coronagraph type   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--Apply the Lyot stop
EP4 = mp.P4.full.croppedMask.*EP4; %padOrCropEven(EP4,mp.P4.full.Narr);

%--MFT from Lyot Stop to final focal plane (i.e., P4 to Fend)
EP4 = propcustom_relay(EP4,mp.NrelayFend,mp.centering); %--Rotate the final image 180 degrees if necessary
EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.full.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

%{
%function [Efoc] = propcustom_mft_PtoF(Epup, f,lambda,dx,dxi,Nxi,deta,Neta,varargin)
 dx = mp.P4.full.dx;
 dxi = mp.Fend.dxi;
 Nxi = mp.Fend.Nxi;
 deta = mp.Fend.deta;
 Neta = mp.Fend.Neta;
 save /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_4jeff dx dxi Nxi deta Neta
 disp('paused-model-full-Fourier'), pause
%}

%--Don't apply FPM if normalization value is being found
if(normFac==0)
    Eout = EFend;  %--Don't normalize if normalization value is being found
else
    Eout = EFend/sqrt(normFac); %--Apply normalization
end

if(mp.useGPU); Eout = gather(Eout); end

if(mp.flagFiber)
    if(mp.flagLenslet)
        Efiber = cell(mp.Fend.Nlens,1);
        sbpIndex = find(mp.sbp_centers == lambda);
        
        for nlens = 1:mp.Fend.Nlens
            EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering,'xfc',mp.Fend.x_lenslet_phys(nlens),'yfc',mp.Fend.y_lenslet_phys(nlens));
            Elenslet = EFend.*mp.Fend.lenslet.mask;
            EF5 = propcustom_mft_PtoF(Elenslet,mp.lensletFL,lambda,mp.Fend.dxi,mp.F5.dxi,mp.F5.Nxi,mp.F5.deta,mp.F5.Neta,mp.centering);
            Efiber{nlens} = mp.F5.fiberMode(:,:,sbpIndex).*sum(sum(mp.F5.fiberMode(:,:,sbpIndex).*conj(EF5)));
        end
        
        Efiber = permute(reshape(cell2mat(Efiber)', mp.F5.Nxi, mp.F5.Neta, mp.Fend.Nlens), [2,1,3]);
        varargout{1} = Efiber;
        
    else  %Fibers placed in the focal plane with no lenslets
        EFend = propcustom_mft_PtoF(EP4,mp.fl,lambda,mp.P4.compact.dx,mp.Fend.dxi,mp.Fend.Nxi,mp.Fend.deta,mp.Fend.Neta,mp.centering);

        sbpIndex = find(mp.sbp_centers == lambda);
        
        Efiber = zeros(mp.Fend.Nxi, mp.Fend.Neta);
        for i=1:mp.Fend.Nfiber
            Eonefiber = mp.Fend.fiberMode(:,:,sbpIndex,i).*sum(sum(mp.Fend.fiberMode(:,:,sbpIndex,i).*conj(EFend)));
            Efiber = Efiber + Eonefiber;
        end
        
        varargout{1} = Efiber;

        figure(901);
        imagesc(log10(abs(Efiber).^2)); axis equal tight; colorbar;
    end
end

end %--END OF FUNCTION
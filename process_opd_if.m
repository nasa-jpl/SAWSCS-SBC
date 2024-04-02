% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Script to perform a DM-apodized VC (DMVC) simple design run.

addpath ~esidick/matlab_erkin/
addpath(genpath('/home/esidick/Afalco/falco20200916/')); %savepath;
addpath('~esidick/Afalco/proper_v3.0.1_matlab_22aug17/'); %savepath;

flag = 1
%disp('paused'), pause

if 1

    %load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_1500nm_30nm_bw20_gap08 dwddm1 indx maskopd
    load ~sab/Projects/6MST/Afalco/falco20200916/macos/IFhex_opd/dwddm_dm1_nm_over_nm_1500nm_opd1_27nm_bw20_gap08 dwddm1 indx maskopd

    %load IFopd/G_dm1_lim3_m2_355nm_30nm_bw20_gap08 km1 km2 
    %load IFopd/G_dm1_dm2_lim3_m3_1500nm_30nm_bw20_gap08 km1 km2 indx

    M = [dwddm1];

    A  = M'*M; 
    sv = diag(A);
    jmax = max(sv); 
    sv = sv / jmax;

    sn = sv;

    lim = 1e-3;
    %lim = 0.02;
    km1 = find(sv > lim);

    ma = 64;
    mx = [1:length(sv)];

    figure(1), clf, semilogy(mx, sort(sv, 'descend'), 'ro-', mx, mx*0+lim,'b-'), grid
        parms = fun_newaxes(14, 2, 4); 
        xlabel('act number index')
        ylabel('eigen-values of M-matrix')
        title('dm1: phase only')
        legend('data','threshold','Location','West')
        print -dpng ../Figs/fig_dm1

    ns = length(km1);
    pp = [ns ma*ma]
    
    dv = reshape(sv, ma,ma);
    mask = (dv > lim);
    nx = [1:ma];
%end
    figure(2), imagesc(nx,nx,log10(dv), [-2 0]), axis image, colormap(jet), 
    %figure(2), imagesc(nx,nx,dv, [0.6 0.8]), axis image, colormap(jet), 
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM1 Act-Strength: N(active) = %i', ns))
        hold on, contour(nx,nx,mask,1,'b','Linewidth',2), hold off
    print -dpng ../Figs/fig_dm1a
    %return

%end % main - if    

if 0

    M = dwddm1(:,km1);

    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;
    
    gi  = -2; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

end

%save -v7.3 ~esidick/Afalco/falco20200916/macos/IFopd/G_matrix G
%save -v7.3 IFopd/G_dm1_lim3_gray_gap08 G km1 indx
%save -v7.3 IFopd/G_dm1_lim5_m3_355nm_30nm_bw20_gap08 G km1 km2 indx
%save -v7.3 IFopd/G_dm1_lim3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx

%return
%end % -----------------------------------------------
%%
%if 0

%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_1500nm_30nm_bw20_gap08 dwddm daddm %indx maskopd
load ~sab/Projects/6MST/Afalco/falco20200916/macos/IFhex_opd/dwddm_dm2_nm_over_nm_1500nm_opd1_27nm_bw20_gap08 dwddm daddm %indx maskopd

M = [dwddm];

    A  = M'*M; 
    sv = diag(A);
    jmax = max(sv); 
    sv = sv / jmax;

    sn = sv;

    km2 = find(sv > lim);

    mx = [1:length(sv)];

    figure(3), clf, semilogy(mx, sort(sv, 'descend'), 'ro-', mx, mx*0+lim,'b-'), grid
        parms = fun_newaxes(14, 2, 4); 
        xlabel('act number index')
        ylabel('eigen-values of M-matrix')
        title('dm2: phase only')
        legend('data','threshold','Location','West')
        print -dpng ../Figs/fig_dm2

    ns = length(km2);
    pp = [ns ma*ma]
    
    dv = reshape(sv, ma,ma);
    mask = (dv > lim);
    nx = [1:ma];
%end
    figure(4), imagesc(nx,nx,log10(dv), [-2 0]), axis image, colormap(jet), 
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM2 Act-Strength: N(active) = %i', ns))
        hold on, contour(nx,nx,mask,1,'b','Linewidth',2), hold off
    print -dpng ../Figs/fig_dm2a

%end % main - if    
%%

%return

    aa = dwddm1(:,km1);
    M = [aa; aa*0];         % nm/nm

    aa = [dwddm(:,km2); daddm(:,km2)];
    M = [M aa];  % nm/nm
        
    A = M'*M;

    sv = diag(A);
    jmax = max(sv);
    sn = sv / jmax;
    
    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;

    end

idd = 2

    gi  = -2; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

save -v7.3 ~sab/Projects/6MST/Afalco/falco20200916/macos/IFhex_opd/G_dm1_dm2_lim1e3_m2_1500nm_30nm_bw20_gap08 G km1 km2 indx

%clear aa M dwddm1 dwddm daddm

return
%end % -----------------------------------------------
%%

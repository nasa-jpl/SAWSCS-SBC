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

if 0
    %load IFopd/G_dm1_dm2_lim5_m2_1064nm G km1 km2 indx

%load ~esidick/Afalco/falco20200916/macos/IFopd/dat_dwddm1_small km1 dwddm1 indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full dwddm1 indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_full_1064nm dwddm daddm indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full_1550nm dwddm1 indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_full_1550nm dwddm daddm indx maskopd
load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_500nm_00nm_bw00_gray_gap08 dwddm1 indx maskopd

%M = [dwddm1(:,km1); daddm];
%M = [dwddm1];
%dw = daddm(:,64*10+32);
%opd = v2m(dw,indx);
%figure(2), clf, imagesc(opd), axis xy image, colormap(jet), colorbar, %return
%figure(1), clf, imagesc(dwddm), colormap(jet), colorbar, return

%end

    A  = M'*M; 
    sv = diag(A);
    jmax = max(sv); 
    sv = sv / jmax;

    sn = sv;

%end
ma = 64;
mx = [1:length(sv)];

    figure(2), clf, semilogy(mx, sort(sv, 'descend'), 'ro-', mx, mx*0+1e-3,'b-'), grid
        parms = fun_newaxes(14, 2, 4); 
        xlabel('act number index')
        ylabel('eigen-values of M-matrix')
        title('dm1: phase only')
        legend('data','threshold','Location','West')
        print -dpng ../Figs/fig_dm1
%return

    lim = 1e-3;

    km1 = find(sv > lim);
    %km2 = km1;

    ns = length(km1);
    pp = [ns ma*ma]
    
    dv = reshape(sv, ma,ma);
    mask = (dv > lim);
    nx = [1:ma];
%end
    figure(1), imagesc(nx,nx,log10(dv), [-8 0]), axis image, colormap(jet), 
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title(sprintf('DM1 Act-Strength: N(active) = %i', ns))
        hold on, contour(nx,nx,mask,1,'b','Linewidth',2), hold off
    print -dpng ../Figs/fig_dm1a

    M = dwddm1(:,km1);

    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;
    
    gi  = -3; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

end

%save -v7.3 ~esidick/Afalco/falco20200916/macos/IFopd/G_matrix G
save -v7.3 IFopd/G_dm1_lim3_gray_gap08 G km1 indx


return
%end % -----------------------------------------------
%%
%load ~esidick/Afalco/falco20200916/macos/IFopd/dat_dwddm1_small km1 dwddm1 indx maskopd
M = dwddm1(:,km2);  % nm/nm
%clear dwddm1
        
    A = M'*M;

    sv = diag(A);
    jmax = max(sv);
    sn = sv / jmax;
    
    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;

%end

    gi  = -1; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

%save -v7.3 ~esidick/Afalco/falco20200916/macos/IFopd/G_matrix G
save -v7.3 IFopd/G_dm2_lim5_m1 G km2

%clear G

return
%end % -----------------------------------------------
%%
flag = 2

load IFopd/G_dm1_dm2_lim5_m2_1064nm km1 km2

%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full_1064nm dwddm1 indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full dwddm1 indx maskopd
load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full_355nm_30nm_bw20_flat dwddm1 indx maskopd
aa = dwddm1(:,km1);
%M = [aa];  % nm/nm
M = [aa; aa*0];  % nm/nm

%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_full_1064nm dwddm daddm indx maskopd
load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_full_355nm_30nm_bw20_flat dwddm daddm indx maskopd
aa = [dwddm(:,km2); daddm(:,km2)];
M = [M aa];  % nm/nm
        
    A = M'*M;

    sv = diag(A);
    jmax = max(sv);
    sn = sv / jmax;
    
    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;

%end

    gi  = -2; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

%save -v7.3 ~esidick/Afalco/falco20200916/macos/IFopd/G_matrix G
%save -v7.3 IFopd/G_dm1_dm2_lim5_m1_1064nm G km1 km2 indx
%save -v7.3 IFopd/G_dm1_lim5_m1_500nm G km1 indx
save -v7.3 IFopd/G_dm1_dm2_lim5_m2_355nm_30nm_bw20_flat G km1 km2 indx

clear aa M dwddm1 dwddm daddm

return
%end % -----------------------------------------------
%%
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_full dwddm1 indx maskopd
%load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm2_nm_over_nm_full dwddm daddm indx maskopd
load ~esidick/Afalco/falco20200916/macos/IFopd/dwddm_dm1_nm_over_nm_500nm_00nm_bw00_gray_gap08 dwddm1 indx maskopd

M = [dwddm1];
        
    A = M'*M;

    sv = diag(A);
    jmax = max(sv);
    sn = sv / jmax;
    
    [u, s, v] = svd(M,0);
    s0 = diag(s);
    sx = s0.^2 / jmax;

%end

    gi  = -3; %[-6 -5 -4 -3 -2 -1];
    gg  = 10 .^ gi;
    G = v * diag( sx ./ ( sx + gg ) ./ s0 ) * u'; 

    load ~esidick/Aoomao/matprog/dat_20230227/opdx_3d_nm_256pix opdx
    opd = opdx(:,:,2);
    dw = m2v(opd, indx);
    du = G * dw;
    V1_nm = reshape(du, 64,64);

    save ~esidick/Aoomao/matprog/dat_20230227/dat_V1_nm V1_nm 
    %save ~esidick/Aoomao/matprog/dat_20230227/dat_V2_nm V2_nm 

    figure(1), clf, imagesc(V1_nm), axis xy image, colormap(jet), colorbar
    

return
%end % -----------------------------------------------
%%



%% EXAMPLE MATLAB RANDOM GENERATION

   addpath ~esidick/matlab_erkin/
   addpath /proj/exospin/users/joegreen/matlab/

   mp.Fend.res = 4;

    fs = 14;

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;
flag_rms = 0;

if 0   
       
pup = mp.P1.full.mask;
[m1 m2 n1 n2] = find_boundary(pup>0, 0.5);
wid1 = m2-m1+1;

lyo = mp.P4.full.croppedMask;
[m1 m2 n1 n2] = find_boundary(lyo>0, 0.5);
wid2 = m2-m1+1;

figure(1), clf, imagesc(pup), axis xy image, colormap(jet), %niceaxes
    parmx = fun_newaxes(14, 0, 0);  
    colorbar, parmx = fun_newaxes(14, 0, 0);  
    title(sprintf('full pupil: dia = %i', wid1))
print -dtiff ~esidick/2022gomap_coro/Figs/fig_pup

figure(2), clf, imagesc(lyo), axis xy image, colormap(jet), %niceaxes
    parmx = fun_newaxes(14, 0, 0);  
    colorbar, parmx = fun_newaxes(14, 0, 0);  
    title(sprintf('full lyo: dia = %i', wid2))
print -dtiff ~esidick/2022gomap_coro/Figs/fig_lyo

return
%end % --------------------------------------------------------
%
x = [5:5:20];
%y1 = [8.68 8.8 8.88 8.94]*1e-8; % pre-control
y1 = [8.67 8.87 8.95 9.02]*1e-7; % pre-control
%y2 = [4.46 9.61 13.1 24.4]*1e-11;  % individual control
y2 = [1.98 4.65 8.84 16.1]*1e-11;  % individual control
%y3 = [6.47 9.64 15.0 24.4]*1e-11; % control bw = 20%
y3 = [3.25 7.19 9.92 16.1]*1e-11; % control bw = 20%

pa = polyfit(x,y1,2);
ya = polyval(pa, x);
pb = polyfit(x,y2,2);
yb = polyval(pb, x);
pc = polyfit(x,y3,2);
yc = polyval(pc, x);

fs = 14;

figure(1), clf, plot(x, y1, 'ro-', x, ya, 'b--'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Scoring Badwidth [$\%$]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Pre-Control: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pa), 'Location', 'Southeast')
   %print -dpng ../Figs/fig_cb

figure(2), clf, plot(x, y2, 'ro-', x, yb, 'b--'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Control and Scoring Badwidth [$\%$]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Post-Control: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pb), 'Location', 'Northwest')
   %print -dpng ../Figs/fig_cb1
    
figure(3), clf, plot(x, y3, 'ro-', x, yc, 'b--'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Scoring Badwidth [$\%$]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Post-Control: Control BW = 20\%, $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pb), 'Location', 'Northwest')
  print -dpng ../Figs/fig_cb1a
return
%end % --------------------------------------------------------
%
load ../dat_bb_wfc2/mean_cb_vs_angle_10mm_05per rx ynom
    y1 = ynom;
load ../dat_bb_wfc2/mean_cb_vs_angle_10mm_10per rx ynom
    y2 = ynom;
load ../dat_bb_wfc2/mean_cb_vs_angle_10mm_15per rx ynom
    y3 = ynom;
load ../dat_bb_wfc2/mean_cb_vs_angle_10mm_20per rx ynom
    y4 = ynom;
x = rx;

figure(1), clf, semilogy(x, y1, 'r-', x, y2, 'b-',x, y3, 'g-',x, y4, 'm-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Mean Contrast vs Field-Angle','Fontsize',fs,'Interpreter', 'Latex');
        axis([2 25 1e-11 1e-8])
        legend('BW = 5%', 'BW = 10%', 'BW = 15%', 'BW = 20%', 'Location', 'North')
   print -dpng ../Figs/fig_cb
return
%end % --------------------------------------------------------
%

load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_05per_3lam_run1_cbx_itr36 dm1 dm2 tx cbx fom Im beta_value
    y1 = cbx; x1 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_10per_5lam_run1_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    y2 = cbx; x2 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_15per_7lam_run1_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    y3 = cbx; x3 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 cbx 
    y4 = cbx; x4 = [1:length(cbx)] - 1;

    s1 = ['BW =  5%: ' sprintf('Last Cb = %0.3g', y1(end))];
    s2 = ['BW = 10%: ' sprintf('Last Cb = %0.3g', y2(end))];
    s3 = ['BW = 15%: ' sprintf('Last Cb = %0.3g', y3(end))];
    s4 = ['BW = 20%: ' sprintf('Last Cb = %0.3g', y4(end))];

    fs = 14;

figure(1), clf, semilogy(x1, y1, 'r-', x2, y2, 'b-',x3, y3, 'g-',x4, y4, 'm-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb (2.4 - 24$\lambda$/D)', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Mean Contrast vs Control Iterations','Fontsize',fs,'Interpreter', 'Latex');
        %axis([2 25 1e-11 1e-8])
        legend(s1,s2,s3,s4, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb

return
%end % --------------------------------------------------------
%%
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    y1 = cbx; x1 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    y2 = cbx; x2 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    y3 = cbx; x3 = [1:length(cbx)] - 1;


    s1 = ['gap = 10mm: ' sprintf('Last Cb = %0.3g', y1(end))];
    s2 = ['gap = 30mm: ' sprintf('Last Cb = %0.3g', y2(end))];
    s3 = ['gap = 50mm: ' sprintf('Last Cb = %0.3g', y3(end))];

    fs = 14;

figure(1), clf, semilogy(x1, y1, 'r-', x2, y2, 'b-',x3, y3, 'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast, Cb (2.4 - 24$\lambda$/D)', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Mean Contrast vs Control Iterations','Fontsize',fs,'Interpreter', 'Latex');
        %axis([2 25 1e-11 1e-8])
        legend(s1,s2,s3, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb

x = [10 30 50];
y1 = [1.62 4.33 6.84]*1e-10;

pa = polyfit(x,y1,1);

x1 = [5:50];
ya = polyval(pa, x1);

figure(2), clf, plot(x, y1, 'ro-', x1, ya, 'b--'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('PM Segment Gap [mm]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Contrast Floor: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        %legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        legend('Data', sprintf('Fit: y = (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
   print -dpng ../Figs/fig_cb1


return
%end % --------------------------------------------------------
%%
% wfc5: changing owa limit in both control and scoring

load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_3lamD_7lam_run5_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
  c1 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_7lamD_7lam_run4_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
  c2 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_12lamD_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
  c3 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_16lamD_7lam_run3_cbx_itr60 dm1 dm2 tx cbx fom Im beta_value
  c4 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run5_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
  c5 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
  c6 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_24lamD_v2_7lam_run6_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
  c7 = cbx(end);

x = [3 6 12 16 20 22 24];
y1 = [c1 c2 c3 c4 c5 c6 c7]
y2 = [3.3 1.17 0.676 0.645 0.776 0.897 c7*1e10] * 1e-10; % 24lamD map evaluated at different OWA

fs = 14;

figure(2), clf, semilogy(x, y1, 'ro-',x,y2,'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('OWA [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('IWA = 2, OWA = varied: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        %legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        %legend('Data', sprintf('Fit: y = (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        legend('varying control/score owa', 'from control-owa = 24\lambda/D', 'Location', 'North')
   %print -dpng ../Figs/fig_cb_change_owa1

%return
%end % --------------------------------------------------------
%%
% wfc6: changing owa, with iwa = 1.5 lam/D, 0.98D-circular lyo, 10mm-gap
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_03lamD_v3_7lam_run3_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
    d1 = cbx(end);
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v2_7lam_run5_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v3_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
  figure(3), clf, semilogy(cbx,'ro-'),grid
    d2 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_12lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
    d3 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_16lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d4 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_20lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d5 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d6 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    d7 = cbx(end);

y3 = [d1 d2 d3 d4 d5 d6 d7];

x3 = [3 7 12 16 20 22 24];

figure(1), clf, semilogy(x, y1, 'ro-',x3,y3,'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('OWA [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('IWA = 2, OWA = varied: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        %legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        %legend('Data', sprintf('Fit: y = (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        legend('iwa = 2.0\lambda/D', 'iwa = 1.5\lambda/D', 'Location', 'North')
   print -dpng ../Figs/fig_cb_change_owa1

return
%end % --------------------------------------------------------
%%
x1 = [3 7 12 16 20 22 24];
y1 = [34.99 35 35.02 34.99 34.63 34.29 34.0];  % iwa = 2 lam/D

x2 = [3 6 12 16 20 22 24];
y2 = [35.03 35.01 35.02 34.99 34.72 34.44 34.08]; % iwa = 1.5 lam/D


figure(1), clf, plot(x1, y1, 'ro-',x2,y2,'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('OWA [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Throughput [%]', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Throughput vs IWA: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        legend('iwa = 2.0\lambda/D', 'iwa = 1.5\lambda/D', 'Location', 'North')
        axis([0 25 33 37])
   print -dpng ../Figs/fig_cb_change_tput


return
%end % --------------------------------------------------------
%%
load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_512_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc7_2048/dat_20mm_cir95_24lamD_7lam_1024_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1;

s1 = sprintf(' 512pix: last Cb = %0.3g', c1(end));
s2 = sprintf('1024pix: last Cb = %0.3g', c2(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'b-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('Effects of PM Diameter Pixel Number: Cb-Ratio = %0.3f', c1(end)/c2(end)),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_pmD


return
%end % --------------------------------------------------------
%%
% wfc6: changing owa, with iwa = 1.5 lam/D, 0.98D-circular lyo, 10mm-gap
% control owa = 24 vs owa = 25, score owa = 24 in both

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_06lamD_v4_7lam_run3_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1;
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_25lamD_7lam_run5_cbx_itr110 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v3_7lam_run5_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1;

s1 = sprintf('c-owa = 6\\lambda/D: last Cb = %0.3g', c1(end));
s2 = sprintf('c-owa = 7\\lambda/D: last Cb = %0.3g', c2(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'b-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 1.5 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('effects of larger control owa on contrast'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_owa


return
%end % --------------------------------------------------------
%%
% wfc6: changing pm aperture, with iwa = 1.5 lam/D, 0.98D-circular lyo, 10mm-gap
% control owa = 24 , score owa = 24 in both
dlamD = 0.25; % lam/D per pix
lim1 = 1;
lim2 = 25;
'g--',
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1;
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_1ring_v3_10mm_cir98_24lamD_7lam_512_run2_cbx_itr170 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx(100:end); x2 = [1:length(c2)] - 1;
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1;
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

s1 = sprintf('monolithic pm: last Cb = %0.3g', c1(end));
s2 = sprintf('1-ring pm: last Cb = %0.3g', c2(end));
s3 = sprintf('2-ring pm: last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'b-', x3, c3, 'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 1.5 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Effects of pm shape on contrast. 20%-BW'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2, s3, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_pm

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'b-', v3, y3, 'g-'), grid
        xline(1.5, 'k--', 'x = 1.5');
        xline(24,  'k--', 'x = 24');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend('monolithic pm', '1-ring pm', '2-ring pm', 'Location', 'North')
        axis([0 25 1e-13 1e-7])
   print -dpng ../Figs/fig_cb_rms

return
%end % --------------------------------------------------------
%%
% wfc8: chg=4, DH R = 1 to 12
dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1; 
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_run3_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1;
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_v2_run2_cbx_itr93 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1;
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
    c4 = cbx; x4 = [1:length(cbx)] - 1;
    [v4, y4] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_2p5to12lamD_chg4_7lam_run2_cbx_itr98 dm1 dm2 tx cbx fom Im beta_value
    c5 = cbx; x5 = [1:length(cbx)] - 1;
    [v5, y5] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_3to12lamD_chg4_7lam_run1_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
    c6 = cbx; x6 = [1:length(cbx)] - 1;
    [v6, y6] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

s1 = sprintf('iwa = 1.0: last Cb = %0.3g', c1(end));
s2 = sprintf('iwa = 1.5-v1: last Cb = %0.3g', c2(end));
s3 = sprintf('iwa = 1.5-v2: last Cb = %0.3g', c3(end));
s4 = sprintf('iwa = 2.0: last Cb = %0.3g', c4(end));
s5 = sprintf('iwa = 2.5: last Cb = %0.3g', c5(end));
s6 = sprintf('iwa = 3.0: last Cb = %0.3g', c6(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'c-',x3,c3,'b-',x4,c4,'m-',x5,c5,'g-',x6,c6,'k'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('Charge = 4, 10mm, cir-lyot-0.98, owa = 12\\lambda/D', c1(end)/c2(end)),'Fontsize',fs);
        legend(s1,s2,s3,s4,s5,s6, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'c-',v3,y3,'b-',v4,y4,'m-',v5,y5,'g-',v6,y6,'k'), grid
        xline(1, 'r--', 'x = 1');
        xline(12,  'r--', 'x = 12');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2,s3,s4,s5,s6, 'Location', 'North')
        axis([0 13 1e-12 1e-7])
   print -dpng ../Figs/fig_cb_rms

   xx = [1:0.5:3];
   yy = [c1(end) c3(end) c4(end) c5(end) c6(end) ];

   figure(3), clf, semilogy(xx, yy, 'ro-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('IWA [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Final Mean Contrast, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('Chrg=4, 10mm, cir-lyot-0.98, owa=12\\lambda/D'), 'Fontsize',fs);
        axis([0.5 3.5 2e-12 6e-12])
   print -dpng ../Figs/fig_cb_vs_iwa1


   
return
%end % --------------------------------------------------------
%% 
% wfc9_psd: with PSD aberrations, chg=4, DH R = 1 to 12

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1; 
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_1to12lamD_chg4_7lam_run5_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; 
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_1to12lamD_chg4_7lam_v2_run5_cbx_itr165 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; 
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
    
s1 = sprintf('(a) aberr-free: last Cb = %0.3g', c1(end));
s2 = sprintf('(b) start from 0: last Cb = %0.3g', c2(end));
s3 = sprintf('(c) start from (a): last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'b-',x3,c3,'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('Charge = 4, 10mm, cir-lyot-0.98, owa = 12\\lambda/D', c1(end)/c2(end)),'Fontsize',fs);
        legend(s1,s2,s3, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'b-',v3,y3,'g-'), grid
        xline(1, 'r--', 'x = 1');
        xline(12,  'r--', 'x = 12');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2,s3, 'Location', 'North')
        axis([0 13 1e-12 1e-7])
   print -dpng ../Figs/fig_cb_rms
    
return
%end % --------------------------------------------------------
%% 
% wfc9_psd: with PSD aberrations, chg=4, DH R = 1 to 12

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
cb = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc8/dat_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cb cbx]; x1 = [1:length(c1)] - 1; 
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg4_7lam_run5_cbx_itr200 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; 
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg6_7lam_run2_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; 
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
    
s1 = sprintf('(a) chrg=4, psd-free: last Cb = %0.3g', c1(end));
s2 = sprintf('(b) chrg=4, with psd: last Cb = %0.3g', c2(end));
s3 = sprintf('(c) chrg=6, with psd: last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x1, c1, 'r-',x2,c2,'b-',x3,c3,'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('10mm, cir-lyot-0.98, dark-hole r = 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s1,s2,s3, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'b-',v3,y3,'g-'), grid
        xline(2, 'm--', 'x = 1');
        xline(12,  'm--', 'x = 12');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s1,s2,s3, 'Location', 'North')
        axis([0 13 1e-12 1e-7])
   print -dpng ../Figs/fig_cb_rms

return
%end % --------------------------------------------------------
%% 

% wfc10_psd: with PSD aberrations, chg=4, DH-R = 2 to 12, gap = 30mm

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_mono_pm_cir98_2to12lamD_chg4_5lam_run1_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c0 = [cbx]; x0 = [1:length(c0)] - 1; t1 = [fom(end)*10];
    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg4_5lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cbx]; x1 = [1:length(c1)] - 1; t1 = [t1 fom(end)*10];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo1000_2to12lamD_chg4_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

cb1 = [c0(end) c1(end) c2(end) c3(end)];    
ta = t1
    
s0 = sprintf('(a) mono-pm:   last Cb = %0.3g', c0(end));
s1 = sprintf('(b) no-sfe:        last Cb = %0.3g', c1(end));
s2 = sprintf('(c) sfe-\\lambda/1000: last Cb = %0.3g', c2(end));
s3 = sprintf('(d) sfe-\\lambda/100:  last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-',x3,c3,'m-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('30mm, cir-lyot-0.98, dark-hole r = 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,s3, 'Location', 'Northeast')
   %print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v0,y0,'r-',v1, y1, 'b-',v2,y2,'g-',v3,y3,'m-'), grid
        xline(2, 'c--', 'x = 1');
        xline(12,  'c--', 'x = 12');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s0,s1,s2,s3, 'Location', 'North')
        axis([0 13 5e-13 1e-7])
   %print -dpng ../Figs/fig_cb_rms

return
%end % --------------------------------------------------------
%% 

% wfc10_psd: with PSD aberrations, chg=4, DH-R = 1 to 12, gap = 30mm

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_mono_pm_cir98_2to12lamD_chg4_5lam_run1_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c0 = [cbx]; x0 = [1:length(c0)] - 1; t1 = [fom(end)*10];
    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_no_psd_2to12lamD_chg4_5lam_run1_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cbx]; x1 = [1:length(c1)] - 1; t1 = [t1 fom(end)*10];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo1000_2to12lamD_chg4_5lam_run1_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo100_2to12lamD_chg4_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

cb2 = [c0(end) c1(end) c2(end) c3(end)];    
    
s0 = sprintf('(a) mono-pm:   last Cb = %0.3g', c0(end));
s1 = sprintf('(b) no-sfe:        last Cb = %0.3g', c1(end));
s2 = sprintf('(c) sfe-\\lambda/1000: last Cb = %0.3g', c2(end));
s3 = sprintf('(d) sfe-\\lambda/100:  last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-',x3,c3,'m-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('10mm, cir-lyot-0.98D, dark-hole r = 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,s3, 'Location', 'Northeast')
   %print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v0,y0,'r-',v1, y1, 'b-',v2,y2,'g-',v3,y3,'m-'), grid
        xline(2, 'c--', 'x = 1');
        xline(12,  'c--', 'x = 12');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s0,s1,s2,s3, 'Location', 'North')
        axis([0 13 5e-13 1e-7])
   %print -dpng ../Figs/fig_cb_rms

   x = [1:4];
   fs = 16;

figure(3), clf, semilogy(x,cb2,'ro-', x, cb1, 'bs-'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('Case Number', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2 - 12$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Cases: mono-pm, no-psd, \lambda/1000-psd, \lambda/100-psd','Fontsize',fs);
        legend('gap = 10mm', 'gap = 30mm', 'Location', 'Northwest')
   print -dpng ../Figs/fig_cb_case

figure(4), clf, plot(x,t1,'ro-', x, ta, 'bs-'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('Case Number', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Absolute Core-Throughput, Tc [%]', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Cases: mono-pm, no-psd, \lambda/1000-psd, \lambda/100-psd','Fontsize',fs);
        legend('gap = 10mm', 'gap = 30mm', 'Location', 'Northeast')
   print -dpng ../Figs/fig_tc_case

return
%end % --------------------------------------------------------
%% gap = 30mm, 5lam-efc vs 7lam-efc, psd = lam/1000, lam/100
cb1 = [3.45 3.12 3.43 3.58] * 1e-11; % psd = lam/1000
cb2 = [8.06 7.27 8.07 8.6] * 1e-11; % psd = lam/100
x = [1:4];

figure(3), clf, semilogy(x,cb1,'ro-', x, cb2, 'bs-'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('Case Number', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2 - 12$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('Cases: 5c/5s, 5c/7s, 7c/7s, 7c/5s','Fontsize',fs);
        legend('\lambda/1000-psd', '\lambda/100-psd', 'Location', 'west')
        text(1.5, 9e-11, 'In title: c = control, s = scoring', 'FontSize', 14)
        axis([1 4 2e-11 1e-10])
   print -dpng ../Figs/fig_cb_case_lam

return
%end % --------------------------------------------------------
%% 






if 1 % wa = [1.5 x] lam/D
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_03lamD_v3_7lam_run3_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v2_7lam_run5_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v3_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_12lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_16lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_20lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
else % wa = [2 x] lam/D
    load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_3lamD_7lam_run5_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_7lamD_7lam_run4_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_12lamD_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_16lamD_7lam_run3_cbx_itr60 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_24lamD_v2_7lam_run6_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value

end

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

lim = 7;
lim1= 1.5;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > lim1) & (rm <= 24);
    maskc1 = (rm >= 1.5) & (rm <= 6);
    maskc2 = (rm >= lim1) & (rm <= lim);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, lim);
    [xc2, yc2] = fun_circle(xv, 24);

    psfn = Im .* maskc2; aa = psfn;
    cb = mean(nonzeros(psfn(:)))
    psfn1 = Im .* maskc1;
    cb1 = mean(nonzeros(psfn(:)));
    
    fs = 14;
    %lim = 6;

    figure(1), clf
        imagesc(xv,xv,log10(Im.*maskc2), [-12 -8]); axis xy image, colormap(jet); 
        %imagesc(xv,xv,log10(aa), [-12 -8]); axis xy image, colormap(jet); 
        %hold on, plot(xc, yc, 'w', xc1, yc1, 'w',  'Linewidth', 2), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        title([sprintf('%0.1f-%i\\lambda/D: Cb = %0.3g', lim1,lim, cb)]);
        axis(25*[-1 1 -1 1])
   print -dpng ../Figs/fig_psfa
return
%end % --------------------------------------------------------
%%

% wfc4: changing pm segment gap, using circular lyot-mask with Dlyo = 0.95D
% Also, for gap = 10mm, change Lyot-mask diameter
load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
  c1 = cbx(end);
  load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_10mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
  c2 = cbx(end);
  load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_30mm_7lam_run4_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
  c3 = cbx(end);
  load /home/esidick/Afalco/falco20200916/dat_bb_wfc4/dat_50mm_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value
  c4 = cbx(end);
  
  x = [0 10 30 50];
  y2 = [c1 c2 c3 c4]
%return
%end % --------------------------------------------------------
%%
% wfc3: changing pm gaps. pm-like Lyot-mask, with gaps 10mm wider than pm at each side. Dlyo = 0.95D
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
  b2 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
  b3 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
  b4 = cbx(end); 
  
  x = [0 10 30 50];
  y1 = [c1 b2 b3 b4]

  fs = 14;

figure(1), clf, semilogy(x, y1, 'ro-', x, y2, 'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('PM Segment Gap [mm]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('PM-Like vs Circular Lyot-Mask: $\lambda$0 = 500nm','Fontsize',fs,'Interpreter', 'Latex');
        %legend('Data', sprintf('Fit: y = (%0.3g)*x^2 + (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        %legend('Data', sprintf('Fit: y = (%0.3g)*x + %0.3g', pa), 'Location', 'Northwest')
        legend('pm-like lyot', 'circular lyot', 'Location', 'Northwest')
   print -dpng ../Figs/fig_cb_cir_pm_like_lyo
return
%end % --------------------------------------------------------
%% dat_bb_wfc11_psd: segment-wise psd, rms = 12 or 6nm 

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;
flag_rms = 0;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/digit_psd_12nm_30mm_cir98_2to12lamD_chg6_9lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    c0 = [cbx]; x0 = [1:length(c0)] - 1; t1 = [fom(end)*10];
    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_digit_rms_12nm_2to12lamD_chg6_4lam_run1_cbx_itr185 dm1 dm2 tx cbx fom Im beta_value
cbx1 = cbx;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_digit_rms_12nm_2to12lamD_chg6_9lam_run1_cbx_itr185 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cbx1]; x1 = [1:length(c1)] - 1; t1 = [t1 fom(end)*10];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_rms_12nm_2to12lamD_chg6_9lam_run1_cbx_itr150 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

cb2 = [c0(end) c1(end) c2(end)];    
    
s0 = sprintf('(a) digit from 0:      last Cb = %0.3g', c0(end));
s1 = sprintf('(b) digit from nom: last Cb = %0.3g', c1(end));
s2 = sprintf('(c) nominal:            last Cb = %0.3g', c2(end));

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('psd-12nm, 30mm, cir-0.98D, 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v0,y0,'r-',v1, y1, 'b-',v2,y2,'g-'), grid
        xline(2, 'c--', 'x = 1');
        xline(12,  'c--', 'x = 12');
        yline(1e-10,  'm--', 'Linewidth', 4);
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s0,s1,s2,'Location', 'North')
        axis([0 13 5e-12 1e-7])
   print -dpng ../Figs/fig_cb_rms
return
%end % --------------------------------------------------------
%%
load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_rms_06nm_start0_2to12lamD_chg6_9lam_run1_cbx_itr158 dm1 dm2 tx cbx fom Im beta_value
    c0 = [cbx]; x0 = [1:length(c0)] - 1; t1 = [fom(end)*10];
    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_rms_06nm_2to12lamD_chg6_9lam_run1_cbx_itr240 dm1 dm2 tx cbx fom Im beta_value
    ns = length(cbx); kk = 152:ns; c1 = [cbx(kk)]; x1 = [1:length(c1)] - 1; t1 = [t1 fom(end)*10];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_rms_12nm_2to12lamD_chg6_9lam_run1_cbx_itr150 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

    a1 = 1.19e-11;
    a2 = 1.40e-11;
    a3 = 2.02e-11;
    c0(end) = a1;
    c1(end) = a2;
    c2(end) = a3;

s0 = sprintf('06nm-rms, from 0:        last Cb = %0.3g', c0(end));
s1 = sprintf('06nm-rms, from 12nm: last Cb = %0.3g', c1(end));
s2 = sprintf('12nm-rms:                    last Cb = %0.3g', c2(end));
s3 = sprintf('aberration-free: last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('psd-06nm, 30mm, cir-0.98D, 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,s3,'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v0,y0,'r-',v1, y1, 'b-',v2,y2,'g-'), grid
        xline(2, 'c--', 'x = 1');
        xline(12,  'c--', 'x = 12');
        yline(1e-10,  'm--', 'Linewidth', 4);
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s0,s1,s2,s3,'Location', 'North')
        axis([0 13 1e-12 1e-7])
   print -dpng ../Figs/fig_cb_rms
return
end % --------------------------------------------------------
%%
   

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_mono_pm_cir98_2to12lamD_chg6_5lam_run1_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c0 = [cbx]; x0 = [1:length(c0)] - 1; t1 = [fom(end)*10];
    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/dm1_psd_30mm_cir98_rms_12nm_2to12lamD_chg6_9lam_run1_cbx_itr150 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
    sub_plot_psf, return

load /home/esidick/Afalco/falco20200916/dat_bb_wfc11_psd/opdz_30mm_cir98_rms_30nm_2to12lamD_chg6_4lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; t1 = [t1 fom(end)*10];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

%end

s0 = sprintf('monolithic: last Cb = %0.3g', c0(end));
s1 = sprintf('no aberr: last Cb = %0.3g', c1(end));
s2 = sprintf('12nm-rms-psd: last Cb = %0.3g', c2(end));
s3 = sprintf('30nm-rms_zern: last Cb = %0.3g', c3(end));

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-',x3,c3,'m-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('35mm-gap, cir-0.98D, 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,s3,'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa

figure(2), clf, semilogy(v0,y0,'r-',v1, y1, 'b-',v2,y2,'g-',v3,y3,'m-'), grid
        xline(2, 'c--', 'x = 1');
        xline(12,  'c--', 'x = 12');
        yline(1e-10,  'm--', 'Linewidth', 4);
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [$\lambda$/D]', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Azimuthal Mean Contrast $ $', 'Fontsize',fs,'Interpreter', 'Latex');
        title(('Azimuthal Mean Contrast vs Field-Angle $ $'),'Fontsize',fs,'Interpreter', 'Latex');
        legend(s0,s1,s2,s3,'Location', 'North')
        axis([0 13 5e-13 1e-7])
   print -dpng ../Figs/fig_cb_rms
return
%end % --------------------------------------------------------
%%
    
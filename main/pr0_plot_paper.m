%% EXAMPLE MATLAB RANDOM GENERATION

   addpath ~esidick/matlab_erkin/
   addpath ~esidick/matlab_joegreen/
   %addpath /proj/exospin/users/joegreen/matlab/

   mp.Fend.res = 4;

    fs = 14;

    flag_rms = 0;

if 0   
    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
        cbb = [cbx(end)]; y1 = cbx; x1 = [1:length(cbx)] - 1;
    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
        cbb = [cbb cbx(end)]; y2 = cbx; x2 = [1:length(cbx)] - 1;
    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
        cbb = [cbb cbx(end)]; y3 = cbx; x3 = [1:length(cbx)] - 1;
    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
        cbb = [cbb cbx(end)]; y4 = cbx; x4 = [1:length(cbx)] - 1;

    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

lim = 24;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2.4) & (rm <= 24);

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));


    thput = fom(end)*10;
    
    fs = 14;

    figure(10), clf
        imagesc(xv,xv,log10(Im.*maskc), [-13 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        ylabel('Field-angle [\lambda/D]');
        title([sprintf('(a) 0pix: Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'gap = 0pix', 'FontSize', 14, 'Color', 'w')
        axis(25*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig3a_psf 

return
%end

    s1 = ['gap = 0pix: ' sprintf('Last Cb = %0.3g', y1(end))];
    s2 = ['gap = 1pix: ' sprintf('Last Cb = %0.3g', y2(end))];
    s3 = ['gap = 3pix: ' sprintf('Last Cb = %0.3g', y3(end))];
    s4 = ['gap = 5pix: ' sprintf('Last Cb = %0.3g', y4(end))];

    fs = 14;

figure(1), clf, semilogy(x1, y1, 'r-', x2, y2, 'b-',x3, y3, 'g-',x4, y4, 'm-'), grid
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('EFC iteration number', 'Fontsize',fs);
        ylabel('Mean Contrast, Cb', 'Fontsize',fs);
        title('Mean Contrast vs Control Iterations','Fontsize',fs);
        axis([0 70 5e-13 1e-4])
        legend(s1,s2,s3, s4,'Location', 'Northeast')
   print -dpng ../Figs_paper/fig2_cb


   return

%end
    %load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    %load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value

    %figure(10), clf, semilogy(cbx, 'ro-'), grid

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));
    thput = fom(end)*10;
    
    figure(1), clf
        imagesc(xv,xv,log10(Im.*maskc), [-13 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        ylabel('Field-angle [\lambda/D]');
        title([sprintf('(d) 5pix: Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'gap = 3pix', 'FontSize', 14, 'Color', 'w')
        axis(25*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig3d_psf    

   return
%end

   x = [0 1 3 5] * 6000/512;

   px = polyfit(x, cbb, 1);
   y = polyval(px, x);

   yy = [35.6 34.3 32.5 30.9];

figure(2), clf, 
[ax, h1, h2] = plotyy(x, cbb, 'rs-', x, y, 'b--', x, y1, 'gd');
grid
hold on, plotyy(x, y1, 'rs-'), hold off
    %parmx = fun_newaxes(14, 2, 7); 
        xlabel('Segmen gap [mm]');
        ylabel('Mean contrast inside 2.4 - 24\lambda/D, Cb');
        title('(b) Mean contrast vs gap size');
        legend(ax(1), 'Data', sprintf('Fit: y = (%0.3g)*x + %0.3g', px), 'Location', 'Southeast')
   %print -dpng ../Figs_paper/fig3b_gap      

%end

fs = 16;

figure(2), clf
[AX,H1,H2] = plotyy(x,cbb,x,yy,'plot');
%You can use the handles returned by plotyy to label the axes and set the line styles used for plotting. With the axes handles you can specify the YLabel properties of the left- and right-side y-axis:
%axis(AX(1))
%set(H1(1), 'Color', 'r')
%ylabel('Mean contrast, Cb','Color','r', 'FontSize', 14)
%set(AX(1), 'YTick', 'Color', 'r')
%if 0
set(get(AX(1),'Ylabel'),'String','Mean contrast, Cb','FontSize', fs)
set(get(AX(2),'Ylabel'),'String','Throughput, Tc [%]','FontSize', fs)

%set(AX(1),'XTick',[212:218]);
set(AX(1),'FontSize',fs);
set(AX(2),'FontSize',fs);

%set(geset(AX(1),'XTick',[212:218]);

%if 0
%set(H1,'Ylabel','String','Mean contrast inside 2.4 - 24\lambda/D, Cb','FontSize', 14)
%set(H2,'Ylabel','String','Core-throughput [%]','FontSize', 14)
%Use the xlabel and title commands to label the x-axis and add a title:


grid
xlabel('Segmen gap [mm]', 'FontSize', fs)
title('Cb and Tc versus gap','FontSize', fs)
%Use the line handles to set the LineStyle properties of the left- and right-side plots:

%if 0
set(H1(1),'LineStyle','-','Linewidth', 2, 'Marker','o','MarkerSize',6);
set(H2(1),'LineStyle','-','Linewidth', 2, 'Marker','s','MarkerSize',6);

legend([H1;H2],'Cb (2.4 - 24\lambda/D)','Tc', 'FontSize', 14, 'Location', 'North')
%legend(AX(2),'thput','FontSize', 14, 'Location', 'North')

%set(H1(2),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerEdgeColor','b');  
%set(H2(1),'LineStyle','none','Marker','*','MarkerSize',6,'MarkerEdgeColor','k');

%set(H1,'LineStyle','rs-')
%set(H2,'LineStyle','bd-')
print -dpng ../Figs_paper/fig3bb_gap      

return
%end % --------------------------------------------------------
% Fig. 6: vs bandwidth
%close all

dlamD = 0.25; % lam/D per pix
lim1 = 1;
lim2 = 25;

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc1/dat_10mm_1lam_run1_cbx_itr20 dm1 dm2 tx cbx fom Im beta_value
    y0 = cbx; x0 = [1:length(cbx)] - 1; %figure(1), clf, semilogy(x0, y0, 'ro-'), grid, return
    [v0, z0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_05per_3lam_run1_cbx_itr36 dm1 dm2 tx cbx fom Im beta_value
    y1 = cbx; x1 = [1:length(cbx)] - 1;
    [v1, z1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_10per_5lam_run1_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    y2 = cbx; x2 = [1:length(cbx)] - 1;
    [v2, z2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_15per_7lam_run1_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
    y3 = cbx; x3 = [1:length(cbx)] - 1;
    [v3, z3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 cbx dm1 dm2 tx cbx fom Im beta_value
    y4 = cbx; x4 = [1:length(cbx)] - 1;
    [v4, z4] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

if 0    
    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    lim = 24;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2.4) & (rm <= 24);

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));
    thput = fom(end)*10;
    fs = 14;

    figure(1), clf
        imagesc(xv,xv,log10(Im.*maskc), [-11 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        ylabel('Field-angle [\lambda/D]');
        title(['(a) 5%: ' sprintf('Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        axis(25*[-1 1 -1 1])
   %print -dpng ../Figs_paper/fig6a_psf  

end

x = [0:5:20];
y = [y0(end) y1(end) y2(end) y3(end) y4(end)];

pb = polyfit(x,y,2);
yb = polyval(pb, x);
pb = pb * 1e11;

figure(2), clf, plot(x, y, 'ro-', x, yb, 'b--'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('Control and scoring badwidth [%]');
        ylabel('Mean contrast, Cb');
        title('(a) Cb vs bandwidth');
        legend('Data', sprintf('Fit: y = [(%0.3g)*x^2 + (%0.3g)*x + %0.3g]*1e-11', pb), 'Location', 'Northwest')
        axis([0 20 0 2e-10])
   print -dpng ../Figs_paper/fig7_bw_new

figure(3), clf, semilogy(v0,z0,'r',v1, z1, 'b-',v2,z2,'g-', v3, z3, 'm-', v4, z4, 'k-'), grid
        xline(2.4, 'k--', 'x = 2.4', 'FontSize', 14);
        xline(24,  'k--', 'x = 24', 'FontSize', 14);
        parmx = fun_newaxes(16, 2, 0); 
        xlabel('Field-angle [\lambda/D]');
        ylabel('Azimuthal mean contrast');
        title(('(b) Cb vs field-angle'));
        legend('BW = 0%', 'BW = 5%', 'BW = 10%', 'BW = 15%', 'BW = 20%', 'Location', 'North')
        axis([0 26 1e-13 1e-6])
   print -dpng ../Figs_paper/fig7b_cb_rms_new
   
return
%end % --------------------------------------------------------
%
% wfc4: changing pm segment gap, using circular lyot-mask with Dlyo = 0.95D
% Also, for gap = 10mm, change Lyot-mask diameter
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_00mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
  c1 = cbx(end); ta1 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
  c2 = cbx(end); ta2 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_30mm_7lam_run4_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
  c3 = cbx(end); ta3 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_50mm_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value
  c4 = cbx(end); ta4 = fom(end)*10;
  
  x = [0 10 30 50];
  y2 = [c1 c2 c3 c4]
  ta = [ta1 ta2 ta3 ta4];

%
% wfc3: changing pm gaps. pm-like Lyot-mask, with gaps 10mm wider than pm at each side. Dlyo = 0.95D
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_10mm_v2_7lam_run3_cbx_itr41 dm1 dm2 tx cbx fom Im beta_value
  b2 = cbx(end); tb2 = fom(end)*10;
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_30mm_v2_7lam_run3_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
  b3 = cbx(end); tb3 = fom(end)*10;
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc3/dat_50mm_7lam_run4_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
  b4 = cbx(end); tb4 = fom(end)*10;
  
  x = [0 1 3 5] * 6000/512;
  y1 = [c1 b2 b3 b4]
  tb = [ta1 tb2 tb3 tb4];

  fs = 14;

figure(1), clf, plot(x, y1, 'ro-', x, y2, 'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Segment gap [mm]');
        title('(a) Cb versus segment gap');
        legend('PM-like Lyot', 'Circular Lyot', 'Location', 'Northwest')
   print -dpng ../Figs_paper/fig8a_cb_cir_pm_like_lyo

figure(2), clf, plot(x, tb, 'ro-', x, ta, 'bs-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('Segment gap [mm]');
        title('(b) Tc versus segment gap');
        legend('PM-like Lyot', 'Circular Lyot', 'Location', 'Northeast')
   print -dpng ../Figs_paper/fig8b_tc_cir_pm_like_lyo

return
%end % --------------------------------------------------------
%%
% Also, for gap = 10mm, change Lyot-mask diameter
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir95_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
  c1 = cbx(end); ta1 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir96_7lam_run3_cbx_itr46 dm1 dm2 tx cbx fom Im beta_value
  c2 = cbx(end); ta2 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir97_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
  c3 = cbx(end); ta3 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir98_7lam_run3_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
  c4 = cbx(end); ta4 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir99_7lam_run4_cbx_itr44 dm1 dm2 tx cbx fom Im beta_value
  c5 = cbx(end); ta5 = fom(end)*10;
  load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc4/dat_10mm_cir100_7lam_run4_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
  c6 = cbx(end); ta6 = fom(end)*10;
  
  x = [0.95:0.01:1];
  y2 = [c1 c2 c3 c4 c5 c6];
  ta = [ta1 ta2 ta3 ta4 ta5 ta6];


fs = 16;

figure(2), clf
[AX,H1,H2] = plotyy(x,y2,x,ta,'plot');
set(get(AX(1),'Ylabel'),'String','Mean contrast, Cb','FontSize', fs)
set(get(AX(2),'Ylabel'),'String','Throughput, Tc [%]','FontSize', fs)
ylim(AX(1), [0 7e-10])
set(AX(1),'YTick',[0:1:7]*1e-10);
ylim(AX(2), [30 36])
set(AX(2),'YTick',[30:1:36]);
set(AX(1),'FontSize',fs);
set(AX(2),'FontSize',fs);
grid
xlabel('Lyot mask diameter [D]', 'FontSize', fs)
title('Cb and Tc versus Lyot diameter','FontSize', fs)
set(H1(1),'LineStyle','-','Linewidth', 2, 'Marker','o','MarkerSize',6);
set(H2(1),'LineStyle','-','Linewidth', 2, 'Marker','s','MarkerSize',6);
legend([H1;H2],'Cb (2.4 - 24\lambda/D)','Tc', 'FontSize', 14, 'Location', 'North')
print -dpng ../Figs_paper/fig9_lyoD      

return
%end % --------------------------------------------------------
%%


% wfc5: changing owa limit in both control and scoring

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_3lamD_7lam_run5_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
  c1 = cbx(end);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_7lamD_7lam_run4_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
  c2 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_12lamD_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
  c3 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_16lamD_7lam_run3_cbx_itr60 dm1 dm2 tx cbx fom Im beta_value
  c4 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run5_cbx_itr80 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc5/dat_10mm_cir98_20lamD_7lam_run5_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
  c5 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
  c6 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_24lamD_v2_7lam_run6_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
  c7 = cbx(end);

x = [3 6 12 16 20 22 24];
y1 = [c1 c2 c3 c4 c5 c6 c7];

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;
    lim = 24;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);

    [xc, yc] = fun_circle(xv, 2);

    maskc = (rm > 2) & (rm <= 24);
    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)))

    figure(1), clf
        imagesc(xv,xv,log10(Im), [-12 -8]); axis image, colormap(jet); niceaxes
        hold on, plot(xc, yc, 'b'), hold off
%return


% wfc6: changing owa, with iwa = 1.5 lam/D, 0.98D-circular lyo, 10mm-gap
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_03lamD_v3_7lam_run3_cbx_itr55 dm1 dm2 tx cbx fom Im beta_value
    d1 = cbx(end); maskc = (rm > 1.5) & (rm <= 3); a = [Im.*maskc];
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc6/dat_10mm_cir98_07lamD_v2_7lam_run5_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_07lamD_v3_7lam_run5_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
  figure(3), clf, semilogy(cbx,'ro-'),grid
    d2 = cbx(end); maskc = (rm > 1.5) & (rm <= 6); a = [a Im.*maskc];
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_12lamD_7lam_run2_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
    d3 = cbx(end); maskc = (rm > 1.5) & (rm <= 12); a = [a Im.*maskc];
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_16lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d4 = cbx(end); maskc = (rm > 1.5) & (rm <= 16); b = [Im.*maskc];
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_20lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d5 = cbx(end); maskc = (rm > 1.5) & (rm <= 20); b = [b Im.*maskc];
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_22lamD_7lam_run4_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
    d6 = cbx(end); 
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    d7 = cbx(end); maskc = (rm > 1.5) & (rm <= 24); b = [b Im.*maskc];

    psf = [a;b];


    maskc = (rm > 1.5) & (rm <= 24);
    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));
    thput = fom(end)*10;
    fs = 14;

    %[xc, yc] = fun_circle(xv, 2);

    close all
   

    figure(1), clf
        imagesc(log10(psf), [-12 -8]); axis image, colormap(jet); niceaxes
        %imagesc(log10(psf), [-11 -8]); colormap(jet); %niceaxes
        %hold on, plot(xc, yc, 'b'), hold off
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        %ylabel('Field-angle [\lambda/D]');
        title(['Contrast map when control OWA is varied']);
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        %axis(25*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig10_new_psf  

return
%end

y3 = [d1 d2 d3 d4 d5 d6 d7];

x3 = [3 7 12 16 20 22 24];

figure(1), clf, semilogy(x, y1, 'ro-',x3,y3,'bs-'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('OWA [\lambda/D]');
        %ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title('(a) Cb versus OWA');
        legend('IWA = 2.0\lambda/D', 'IWA = 1.5\lambda/D', 'Location', 'North')
   print -dpng ../Figs_paper/fig9a_cb_change_owa

x1 = [3 7 12 16 20 22 24];
y1 = [34.99 35 35.02 34.99 34.63 34.29 34.0];  % iwa = 2 lam/D

x2 = [3 6 12 16 20 22 24];
y2 = [35.03 35.01 35.02 34.99 34.72 34.44 34.08]; % iwa = 1.5 lam/D

figure(2), clf, plot(x1, y1, 'ro-',x2,y2,'bs-'), grid
    parmx = fun_newaxes(16, 2, 5); 
        xlabel('OWA [\lambda/D]');
        %ylabel('Throughput [%]', 'Fontsize',fs,'Interpreter', 'Latex');
        title('(b) Tc versus OWA');
        legend('IWA = 2.0\lambda/D', 'IWA = 1.5\lambda/D', 'Location', 'North')
        axis([0 25 33 37])
   print -dpng ../Figs_paper/fig9b_tc_change_owa



return
%end % --------------------------------------------------------
%
% wfc6: changing pm aperture, with iwa = 1.5 lam/D, 0.98D-circular lyo, 10mm-gap
% control owa = 24 , score owa = 24 in both
dlamD = 0.25; % lam/D per pix
lim1 = 1;
lim2 = 25;

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_00mm_cir98_24lamD_7lam_run2_cbx_itr50 dm1 dm2 tx cbx fom Im beta_value
    c1 = cbx; x1 = [1:length(cbx)] - 1; thput = fom(end)*10;
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_1ring_v3_10mm_cir98_24lamD_7lam_512_run2_cbx_itr170 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx(100:end); x2 = [1:length(c2)] - 1; thput = fom(end)*10;
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc6/dat_10mm_cir98_24lamD_7lam_run6_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; thput = fom(end)*10;
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;
    lim = 24;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);

    maskc = (rm > 1.5) & (rm <= 24);
    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));

    figure(2), clf
        imagesc(xv,xv,log10(Im.*maskc), [-12 -8]); axis xy image, colormap(jet); %niceaxes
        %imagesc(xv,xv,log10(Im.*maskc), [-11 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        ylabel('Field-angle [\lambda/D]');
        title([sprintf('(c) Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        axis(25*[-1 1 -1 1])
   %print -dpng ../Figs_paper/fig12_v2_psf3

pupil = fitsread('/home/esidick/2022habex/dat/pupil_10mm_512pix_soft.fits');
load ~esidick/Afalco/falco20200916/erkin_iris/masks_6mst/lyo_512pix_098D lyo
aa = 2*pupil - lyo;

close all

    figure(3), clf
        imagesc(aa); axis xy image, colormap(gray); niceaxes
   print -dpng ../Figs_paper/fig12_v2_mask3


   return

   close all

s1 = sprintf('monolithic pm: last Cb = %0.3g', c1(end));
s2 = sprintf('1-ring pm: last Cb = %0.3g', c2(end));
s3 = sprintf('2-ring pm: last Cb = %0.3g', c3(end));

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'b-', v3, y3, 'g-'), grid
        xline(1.5, 'k--', 'x = 1.5', 'FontSize', 14);
        xline(24,  'k--', 'x = 24', 'FontSize', 14);
        parmx = fun_newaxes(16, 2, 0); 
        xlabel('Field-angle [\lambda/D]');
        ylabel('Azimuthal mean contrast');
        title(('Azimuthal mean contrast vs field-angle'));
        legend('monolithic PM', '1-ring PM', '2-ring PM', 'Location', 'North')
        axis([0 26 1e-13 1e-6])
   print -dpng ../Figs_paper/fig13_cb_rms

return
%end % --------------------------------------------------------
%% Changing IWA
% wfc8: chg=4, DH R = 1 to 12
dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;
    lim = 24;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);

    pn = 100;

    flag_rms = 0;


load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_1to12lamD_chg4_7lam_run2_cbx_itr90 dm1 dm2 tx cbx fom Im beta_value
c1 = cbx; x1 = [1:length(cbx)] - 1; thput = [fom(end)*10]; maskc = (rm > 1) & (rm <= 12); psf1 = [pad(Im.*maskc, pn)];
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_run3_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1;
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_1p5to12lamD_chg4_7lam_v2_run2_cbx_itr93 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; thput = [thput fom(end)*10]; maskc = (rm > 1.5) & (rm <= 12); psf1 = [psf1 pad(Im.*maskc, pn)];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value
    c4 = cbx; x4 = [1:length(cbx)] - 1; thput = [thput fom(end)*10]; maskc = (rm > 2) & (rm <= 12); psf2 = [pad(Im.*maskc, pn)];
    [v4, y4] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_2p5to12lamD_chg4_7lam_run2_cbx_itr98 dm1 dm2 tx cbx fom Im beta_value
    c5 = cbx; x5 = [1:length(cbx)] - 1; thput = [thput fom(end)*10]; 
    [v5, y5] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_3to12lamD_chg4_7lam_run1_cbx_itr45 dm1 dm2 tx cbx fom Im beta_value
    c6 = cbx; x6 = [1:length(cbx)] - 1; thput = [thput fom(end)*10]; maskc = (rm > 3) & (rm <= 12); psf2 = [psf2 pad(Im.*maskc, pn)];
    [v6, y6] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

psf = [psf1; psf2];


    figure(1), clf
        imagesc(log10(psf), [-12 -8]); axis image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        %ylabel('Field-angle [\lambda/D]');
        title('IWA = 1, 1.5, 2, 3\lambda/D; OWA = 12\lambda/D');
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        %axis(25*[-1 1 -1 1])
   %print -dpng ../Figs_paper/fig12aa_4psf  
%return

   x = [1:0.5:3];
   y2 = [c1(end) c3(end) c4(end) c5(end) c6(end) ];

   close all

fs = 16;

figure(2), clf
[AX,H1,H2] = plotyy(x,y2,x,thput,'plot');
set(get(AX(1),'Ylabel'),'String','Mean contrast, Cb','FontSize', fs)
set(get(AX(2),'Ylabel'),'String','Throughput, Tc [%]','FontSize', fs)
ylim(AX(1), [2e-12 6e-12])
set(AX(1),'YTick',[2:1:6]*1e-12);
ylim(AX(2), [38 42])
set(AX(2),'YTick',[38:1:42]);
set(AX(1),'FontSize',fs);
set(AX(2),'FontSize',fs);
grid
xlabel('IWA [\lambda/D]', 'FontSize', fs)
title('Cb and Tc versus IWA','FontSize', fs)
set(H1(1),'LineStyle','-','Linewidth', 2, 'Marker','o','MarkerSize',6);
set(H2(1),'LineStyle','-','Linewidth', 2, 'Marker','s','MarkerSize',6);
legend([H1;H2],'Cb','Tc', 'FontSize', 14, 'Location', 'North')
print -dpng ../Figs_paper/fig13_iwa      

return
%end % --------------------------------------------------------
% wfc9_psd: dm digitization, chg=4,6, DH R = 1 to 12
%close all

% 9-surface psd contronl:
%load ~esidick/2022habex/dat/dat_wfc_residual_psd psd_nm dm1  


dlamD = 0.25; % lam/D per pix
lim1 = 1;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dm1_psd_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr205 dm1 dm2 tx cbx fom Im beta_value
    y1 = cbx; x1 = [1:length(cbx)] - 1; tc = [fom(end)*10];
    [v1, z1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_digit_v2_chg4_7lam_run1_cbx_itr225 dm1 dm2 tx cbx fom Im beta_value
    y2 = cbx; x2 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [v2, z2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg6_7lam_run2_cbx_itr260 dm1 dm2 tx cbx fom Im beta_value
    y3 = cbx; x3 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [v3, z3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_digit_v2_chg6_7lam_run1_cbx_itr290 cbx dm1 dm2 tx cbx fom Im beta_value
    y4 = cbx; x4 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [v4, z4] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);

if 0    
    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    lim = 12;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = Im .* maskc; aa = psfn;
    cb = mean(nonzeros(psfn(:)));
    thput = fom(end)*10;
    fs = 14;

    figure(1), clf
        imagesc(xv,xv,log10(Im.*maskc), [-11 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        ylabel('Field-angle [\lambda/D]');
        title([sprintf('(d) Charge-6, DM digitized\nCb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %xlabel([sprintf('Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        axis(13*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig17d_psf  
return
end

x = [1:4];
y = [y1(end) y2(end) y3(end) y4(end)];

figure(3), clf, semilogy(v1, z1, 'r-',v2,z2,'b-', v3, z3, 'g-', v4, z4, 'm-'), grid
        xline(2, 'k--', 'x = 2', 'FontSize', 14);
        xline(12,  'k--', 'x = 12', 'FontSize', 14);
        parmx = fun_newaxes(16, 2, 0); 
        xlabel('Field-angle [\lambda/D]');
        title(('(b) Azimuthal Cb vs field-angle'));
        legend('chrg-4, not-ditized', 'chrg-4, ditized', 'chrg-6, not-ditized', 'chrg-6, ditized', 'Location', 'North')
        axis([1 14 8e-11 1e-6])
   %print -dpng ../Figs_paper/fig18b_cb_rms

   close all

fs = 16;

figure(1), clf
[AX,H1,H2] = plotyy(x,y,x,tc,'plot');
set(get(AX(1),'Ylabel'),'String','Mean contrast, Cb','FontSize', fs)
set(get(AX(2),'Ylabel'),'String','Throughput, Tc [%]','FontSize', fs)
if 0
    ylim(AX(1), [2e-12 6e-12])
    set(AX(1),'YTick',[2:1:6]*1e-12);
    ylim(AX(2), [38 42])
    set(AX(2),'YTick',[38:1:42]);
end
set(AX(1),'FontSize',fs);
set(AX(2),'FontSize',fs);
grid
xlabel('Case Number', 'FontSize', fs)
title('(a) Cb and Tc of four cases','FontSize', fs)
set(H1(1),'LineStyle','-','Linewidth', 2, 'Marker','o','MarkerSize',6);
set(H2(1),'LineStyle','-','Linewidth', 2, 'Marker','s','MarkerSize',6);
legend([H1;H2],'Cb (2 - 12\lambda/D)','Tc', 'FontSize', 14, 'Location', 'North')
print -dpng ../Figs_paper/fig18a_cb_tc   


return
%end % --------------------------------------------------------
%%
% Fig. 22: DM commands

load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
dm0 = dm1*1e9;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dm1_psd_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr205 dm1 dm2 tx cbx fom Im beta_value
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value

    aa = dm1; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1]
    figure(1), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        %ylabel('DM1', 'FontSize', 20)
        title(sprintf('dm1: not-flattened, rms = %0.1fnm', rms0(1)));
   print -dpng ../Figs_paper/fig22_v2_dma1 

    aa = dm2; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1]
    figure(2), clf
        imagesc(aa,cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        %ylabel('DM2', 'FontSize', 20)
        title(sprintf('dm2: not-flattened, rms = %0.1fnm', rms0(1)));
   print -dpng ../Figs_paper/fig22_v2_dma2 

return
%end % --------------------------------------------------------
%

load  ~esidick/2022habex/dat/dat_wfc_residual_psd_11srf_lamo100 del_psd psd_nm indx_opd dm1
%aa = psd_nm; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
aa = del_psd; rms0 = stat2d(aa); cx = rms0(1)*3*[-1 1];
    figure(1), clf
        imagesc(aa, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title([sprintf('(b) post-flattening wfe: rms = %0.1fnm', rms0(1))]);
   print -dpng ../Figs_paper/fig19b_pos_wfe  

return
%end % --------------------------------------------------------
%% 

% wfc10_psd: with PSD aberrations, chg=4, DH-R = 2 to 12, gap = 30mm

dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

flag_rms = 0;

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dat_mono_pm_cir98_2to12lamD_chg4_5lam_run1_cbx_itr65 dm1 dm2 tx cbx fom Im beta_value
%    c0 = [cbx]; x0 = [1:length(c0)] - 1; tc = fom(end)*10;
%    [v0, y0] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2);
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr140 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cbx]; x1 = [1:length(c1)] - 1; tc = fom(end)*10;
    [v1, y1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo1000_2to12lamD_chg6_5lam_run1_cbx_itr95 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo1000_2to12lamD_chg6_5lam_run1_cbx_itr125 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [v2, y2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [v3, y3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

cb = [c1(end) c2(end) c3(end)];

if 0

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
    c1 = [cbx]; x1 = [1:length(c1)] - 1; tc = [fom(end)*10];
    [w1, z1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo1000_2to12lamD_chg6_5lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo1000_2to12lamD_chg6_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value
    c2 = cbx; x2 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [w2, z2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value
    c3 = cbx; x3 = [1:length(cbx)] - 1; tc = [tc fom(end)*10];
    [w3, z3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2,flag_rms);

cb2 = [c1(end) c2(end) c3(end)];

end

if 0    

    cb1 = cbx(end);

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    lim = 12;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = Im .* maskc; aa = psfn;
    cb1 = mean(nonzeros(psfn(:)));
    cb = cbx(end);
    thput = fom(end)*10;
    fs = 14;

    pp = [cb1 cb]

    figure(3), clf
        imagesc(xv,xv,log10(Im.*maskc), [-12 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        ylabel('Field-angle [\lambda/D]');
        title([sprintf('(c) \\lambda/100: Cb=%0.3g, Tc=%0.1f', cb, thput) '%']);
        %title([sprintf('(a) no sfe: Cb=%0.3g, Tc=%0.1f', cb, thput) '%']);
        %xlabel([sprintf('Cb = %0.3g, Tc = %0.1f', cb, thput) '%']);
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        axis(13*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig20_v2_psf3  

figure(2), clf, semilogy(cbx, 'r-'), grid   

return 
end

x = [1 2 3];

fs = 16;

figure(1), clf
[AX,H1,H2] = plotyy(x,cb,x,tc,'plot');
set(get(AX(1),'Ylabel'),'String','Mean contrast, Cb','FontSize', fs)
set(get(AX(2),'Ylabel'),'String','Throughput, Tc [%]','FontSize', fs)
if 1
    ylim(AX(1), [0 10e-11])
    set(AX(1),'YTick',[0:2:10]*1e-11);
    ylim(AX(2), [29 35])
    set(AX(2),'YTick',[29:1:35]);
end
set(AX(1),'FontSize',fs);
set(AX(2),'FontSize',fs);
grid
xlabel('Case Number', 'FontSize', fs)
title('Cb and Tc of three cases in Fig.23','FontSize', fs)
set(H1(1),'LineStyle','-','Linewidth', 2, 'Marker','o','MarkerSize',6);
set(H2(1),'LineStyle','-','Linewidth', 2, 'Marker','s','MarkerSize',6);
legend([H1;H2],'Cb (2 - 12\lambda/D)','Tc', 'FontSize', 14, 'Location', 'North')
print -dpng ../Figs_paper/fig24_cb_tc   

return

    
%s0 = sprintf('(a) mono-pm:   last Cb = %0.3g', c0(end));
s1 = sprintf('(a) no-sfe:        Cb = %0.3g', c1(end));
s2 = sprintf('(b) sfe-\\lambda/1000: Cb = %0.3g', c2(end));
s3 = sprintf('(c) sfe-\\lambda/100:   Cb = %0.3g', c3(end));

figure(2), clf, semilogy(v1, y1, 'r-',v2,y2,'b-',v3,y3,'g-'), grid
hold on, semilogy(w1,z1,'r--',w2, z2, 'b--',w3,z3,'g--'), hold off
        xline(2, 'k--');
        xline(12,  'k--');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [\lambda/D]');
        ylabel('Azimuthal Mean Contrast');
        title(('Solid: gap = 11.7mm. Dashed: gap = 35.2mm'));
        legend(s1,s2,s3, 'Location', 'North')
        axis([0 13 2e-12 1e-8])
print -dpng ../Figs_paper/fig25_cb_radial   

   return

figure(1), clf, semilogy(x0,c0,'r-', x1, c1, 'b-',x2,c2,'g-',x3,c3,'m-'), grid
    parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations', 'Fontsize',fs,'Interpreter', 'Latex');
        ylabel('Mean Contrast Inside 2.4 - 24$\lambda$/D, Cb', 'Fontsize',fs,'Interpreter', 'Latex');
        title(sprintf('10mm, cir-lyot-0.98D, dark-hole r = 2 - 12\\lambda/D'),'Fontsize',fs);
        legend(s0,s1,s2,s3, 'Location', 'Northeast')
   print -dpng ../Figs/fig_cb_vs_iwa


return
%end % --------------------------------------------------------
%% 
% wfc9_psd: with PSD aberrations, chg=4, DH-R = 2 to 12, gap = 10mm
% wfc10_psd: with PSD aberrations, chg=4, DH-R = 2 to 12, gap = 10mm
load /home/esidick/Afalco/falco20200916/dat_bb_wfc9_psd/dat_10mm_cir98_2to12lamD_chg4_7lam_run5_cbx_itr200 dm1 dm2 tx cbx fom Im beta_value
aa = dm1;
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_10mm_cir98_lamo100_2to12lamD_chg4_5lam_run1_cbx_itr130 dm1 dm2 tx cbx fom Im beta_value

bb = dm1; rms0 = stat2d(bb); cx = rms0(1)*3*[-1 1];
    figure(2), clf
        imagesc(bb, cx); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        title([sprintf('(b) with wfe-flattening: rms = %0.1fnm', rms0(1))]);
   print -dpng ../Figs_paper/fig22b_dm1  

return
%end % --------------------------------------------------------
%% 
dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

flag_rms = 0;

%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_5lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    cb3 = cbx(end); t3 = fom(end)*10; bb1 = Im; cb3 = cbx(end);
    [v1, z1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_7lam_run1_cbx_itr70 dm1 dm2 tx cbx fom Im beta_value
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_7lam_run1_cbx_itr110 dm1 dm2 tx cbx fom Im beta_value
    cb4 = cbx(end); t4 = fom(end)*10; bb2 = Im; cb4 = cbx(end),
%load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_lamo100_2to12lamD_chg6_9lam_run1_cbx_itr30 dm1 dm2 tx cbx fom Im beta_value
    cb4 = cbx(end); t5 = fom(end)*10; bb3 = Im; cb5 = cbx(end),
    
    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;

    lim = 12;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = abs(bb3) .* maskc; 
    rms0 = stat2d(psfn);

if 0    
    figure(1), clf
        imagesc(xv,xv,log10(psfn), [-12 -8]); axis xy image, colormap(jet); niceaxes
        %imagesc(log10([Im bb1])); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        ylabel('Control', 'FontSize', 20);
        %title([sprintf('|C(9\\lambda) - C(5\\lambda)|: rms = %0.3g', rms0(1))]);
        title([sprintf('7\\lambda: Cb = %0.3g', cb4)]);
        axis(13*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig26a_psf  
   return
end

if 1

    flag_rms = 1;

load ../dat_paper/fig26_Im_5c_7s_lamo100_5lam_30mm Im xv cb maskc 
    aa = abs(Im - bb1); Im0 = Im;
    [v2, z2] = fun_rms_azimuthal_mean_contrast(aa, dlamD, lim1, lim2, flag_rms);
    aa = aa .* maskc;
    rmsa = stat2d(aa);

load ../dat_paper/fig26_Im_5c_9s_lamo100_5lam_30mm_v2 Im xv cb maskc 
    bb = abs(Im - bb1);
    [v3, z3] = fun_rms_azimuthal_mean_contrast(bb, dlamD, lim1, lim2, flag_rms);
    bb = bb .* maskc;
    rmsb = stat2d(bb);

end


    figure(1), clf
        imagesc(xv,xv,log10(bb), [-13 -9]); axis xy image, colormap(jet); niceaxes
        %imagesc(log10([Im0 bb1])); axis xy image, colormap(jet); %niceaxes
        %imagesc([abs(Im0- bb1).*maskc abs(Im- bb1).*maskc]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        ylabel('Field-angle [\lambda/D]');
        %ylabel('Scoring', 'FontSize', 20);
        title([sprintf('(b) |C(9\\lambda) - C(5\\lambda)|: rms = %0.2g', rmsb(1))]);
        %xlabel([sprintf('9\\lambda: Cb = %0.3g', cb)]);
        axis(13*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig27b_psf  

   return

   close all

figure(3), clf, semilogy(v2, z2, 'r-',v3,z3,'b-', v1, z1, 'g-'), grid
        xline(2, 'k--', 'x = 2', 'FontSize', 14);
        xline(12,  'k--', 'x = 12', 'FontSize', 14);
        parmx = fun_newaxes(16, 2, 0); 
        xlabel('Field-angle [\lambda/D]');
        ylabel('Azimuthal rms of \Deltacontrast');
        %title(('(b) Azimuthal Cb vs field-angle'));
        legend('rms, C(7\lambda) - C(5\lambda)', 'rms, C(9\lambda) - C(5\lambda)', 'contrast, C(5\lambda)', 'Location', 'North')
        axis([1 13 1e-12 1e-7])
   print -dpng ../Figs_paper/fig28_rms_cb

return
%end % --------------------------------------------------------
%% 
load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg4_5lam_run1_cbx_itr140 dm1 dm2 tx cbx fom Im beta_value
cba = [cbx(end)];
load ../dat_paper/dat_fig27_01lamD_chrg4_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cba = [cba cb]; 
load ../dat_paper/dat_fig27_05lamD_chrg4_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cba = [cba cb];
load ../dat_paper/dat_fig27_10lamD_chrg4_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cba = [cba cb];

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
cbb = [cbx(end)];
load ../dat_paper/dat_fig27_01lamD_chrg6_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cbb = [cbb cb];
load ../dat_paper/dat_fig27_05lamD_chrg6_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cbb = [cbb cb]; 
load ../dat_paper/dat_fig27_10lamD_chrg6_Im Im cb maskc xv
psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)));
cbb = [cbb cb]; 

load ../dat_paper/chrg4_cbxx_vs_lamD mx cbxx
x1 = mx; y1 = cbxx;
load ../dat_paper/chrg6_cbxx_vs_lamD mx cbxx
x2 = mx; y2 = cbxx;
load ../dat_paper/z1z2z3_chrg4_cbxx_vs_lamD mx cbxx psfm
x3 = mx; y3 = cbxx;
load ../dat_paper/z1z2z3_chrg6_cbxx_vs_lamD mx cbxx psfm
x4 = mx; y4 = cbxx;

x = [0 0.01 0.05 0.1];

figure(1), clf, semilogy(x1, y1, 'ro-', x3, y3, 'rs--',x2,y2,'bd-',x4,y4,'b^--'), grid
        parms = fun_newaxes(16, 2, 7);  
        xlabel('Stellar angular diameter [\lambda/D]');
        ylabel('Mean contrast, Cb');
        legend('Charge 4, no z2/z3 control', 'Charge 4, with z2/z3 control','Charge 6, no z2/z3 control','Charge 6, with z2/z3 control', 'Location', 'northwest')
        axis([0 0.1 1e-11 1e-6])
   print -dpng ../Figs_paper/fig32_cb_vs_star
return

cb

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / mp.Fend.res;

lim = 12;

    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);
    
    [xc, yc]  = fun_circle(xv, 2);
    [xc1, yc1]  = fun_circle(xv, lim);
    [xc2, yc2] = fun_circle(xv, 24);

    psfn = Im .* maskc; aa = psfn;
    psfn = Im .* maskc; cb = mean(nonzeros(psfn(:)))
    cb = cbx(end);
    
    fs = 16;

    figure(2), clf
        imagesc(xv,xv,log10(Im.*maskc), [-10 -6]); axis xy image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        ylabel('Charge 6', 'FontSize', 20);
        title([sprintf('0.10\\lambda/D')], 'FontSize', 20);
        axis(13*[-1 1 -1 1])
%print -dpng ../Figs_paper/fig27f_chrg6_10lamD

return
%end % --------------------------------------------------------
%% 
% wfc10_psd: with PSD aberrations = lam/100, chg=6, DH-R = 2 to 12, gap = 30mm

y3 = [3.54 3.55 3.6 3.6 3.58 3.57 3.56]*1e-11; x3  =[3:9];
y4 = [3.59 3.2 3.14 3.09 3.07 3.05]*1e-11; x4 = [4:9];
y5 = [4.15 3.67 3.59 3.53 3.34 ]*1e-11; x5 = [5:9];

figure(1), clf, plot(x3, y3, 'ro-', x4, y4, 'bs-', x5, y5, 'gd-'), grid
        parms = fun_newaxes(16, 2, 7);  
        xlabel('Number of scoring passbands');
        ylabel('Mean contrast, Cb');
        legend('3\lambda-control', '4\lambda-control', '5\lambda-control', 'Location', 'northeast')
        %axis([1 13 1e-12 1e-7])
   print -dpng ../Figs_paper/fig29_cb_vs_passband_number

return
%end % --------------------------------------------------------
%% 
load /home/esidick/Afalco/falco20200916/dat_paper/jeff_apodization_mask.mat
mask = ap_mask;

load /home/esidick/Afalco/falco20200916/dat_paper/jeff_OnAxis_and_7LambdaoverD.mat
Im = QContrast;

[nr, nc] = size(Im);
xv = [-nr/2:nr/2-1] / 4;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);
    maskc = (rm > 2) & (rm <= 12);

    psfn = Im .* maskc; 
    cb = mean(nonzeros(psfn(:)))

    figure(1), clf
        imagesc(mask); axis xy image, colormap(gray); niceaxes
        parms = fun_newaxes(16, 0, 0);  
        colorbar, parms = fun_newaxes(16, 0, 0);  
        %ylabel('Charge 6', 'FontSize', 20);
        title('(a) Apodizer', 'FontSize', 20);
        %axis(13*[-1 1 -1 1])
print -dpng ../Figs_paper/fig32a_apod


figure(2), clf
        imagesc(xv,xv,log10(Im), [-12 -8]); axis xy image, colormap(jet); %niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        ylabel('Field-Angle [\lambda/D]', 'FontSize', 20);
        title(sprintf('(b) Narrowband Contrast: Cb = %0.3g', cb), 'FontSize', 20);
        axis(13*[-1 1 -1 1])
print -dpng ../Figs_paper/fig32b_psf

    
return
%end % --------------------------------------------------------
%% 
dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

flag_rms = 0;

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg4_5lam_run1_cbx_itr140 dm1 dm2 tx cbx fom Im beta_value
    cb1 = cbx(end); t1 = fom(end)*10; bb1 = Im; cb3 = cbx(end);
    [v1, z1] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/dm1_psd_30mm_cir98_no_psd_2to12lamD_chg6_5lam_run1_cbx_itr160 dm1 dm2 tx cbx fom Im beta_value
    cb2 = cbx(end); t2 = fom(end)*10; 
    [v2, z2] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/z1z2z3_30mm_cir98_no_psd_2to12lamD_chg4_4lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    cb3 = cbx(end); t3 = fom(end)*10; 
    [v3, z3] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

load /home/esidick/Afalco/falco20200916/dat_bb_wfc10_psd/z1z2z3_30mm_cir98_no_psd_2to12lamD_chg6_4lam_run1_cbx_itr100 dm1 dm2 tx cbx fom Im beta_value
    cb4 = cbx(end); t4 = fom(end)*10; 
    [v4, z4] = fun_rms_azimuthal_mean_contrast(Im, dlamD, lim1, lim2, flag_rms);

    cb = [cb1 cb2; cb3 cb4]
    tc = [t1 t2; t3 t4]

figure(2), clf, semilogy(v1, z1, 'r-',v3,z3,'r--',v2,z2,'b-', v4,z4,'b--'), grid
%hold on, semilogy(w1,z1,'r--',w2, z2, 'b--',w3,z3,'g--'), hold off
        xline(2, 'k--');
        xline(12,  'k--');
        parmx = fun_newaxes(14, 2, 0); 
        xlabel('Field-Angle [\lambda/D]');
        ylabel('Azimuthal Mean Contrast');
        title('Stellar angular diameter = 0\lambda/D');
        legend('Charge 4, no z2/z3 control', 'Charge 4, with z2/z3 control','Charge 6, no z2/z3 control','Charge 6, with z2/z3 control', 'Location', 'north')
        axis([0 13 2e-12 1e-8])
print -dpng ../Figs_paper/fig33_cb_radial   

    
return
%end % --------------------------------------------------------
%% 
x = [1 2 3];
y = [4.32e-13 1.07e-10 1.36e-10];
figure(1), clf, semilogy(x, y, 'ro-'), grid
        parmx = fun_newaxes(14, 2, 5); 
        xlabel({'(a)                                      (b)                                      (c)',  'Case Number'});
        %xlabel({'\fontsize{15}Population','\fontsize{20}(in thousands)'})
        ylabel('Mean Contrast, Cb');
        axis([1 3 1e-13 1e-9])
print -dpng ../Figs_paper/fig40_cb_case3   

    
return
end % --------------------------------------------------------
%% 
dlamD = 0.25; % lam/D per pix
lim1 = 0.5;
lim2 = 13;

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc8/dat_10mm_cir98_2to12lamD_chg4_7lam_run1_cbx_itr75 dm1 dm2 tx cbx fom Im beta_value

    [nr, nc] = size(Im);
    xv = [-nr/2:nr/2-1] / mp.Fend.res;
    lim = 24;
    [xm, ym] = meshgrid(xv);
    rm = abs(xm + 1i * ym);

    pn = 100;

    flag_rms = 0;

    c4 = cbx; x4 = [1:length(cbx)] - 1; thput = [thput fom(end)*10]; maskc = (rm > 2) & (rm <= 12); psf1 = [pad(Im.*maskc, pn)];

load /home/esidick/Afalco/falco20200916/dat_wfc/dat_bb_wfc5/dat_10mm_cir98_12lamD_7lam_run3_cbx_itr35 dm1 dm2 tx cbx fom Im beta_value
  x1 = [1:length(c4)] - 1;
  x2 = [1:length(cbx)] - 1;
figure(2), clf, semilogy(x1, c4, 'r-',x2,cbx,'b-'), grid
        parmx = fun_newaxes(14, 2, 5); 
        xlabel('EFC Iterations');
        ylabel('Mean Contrast, Cb');
        %axis([1 3 1e-13 1e-9])
        legend('IWA: Further EFC from 1-12\lambda/D dark-hole case, chg-4', 'OWA: From start, EFC not exausted, chg-6', 'Location', 'Northeast')
print -dpng ../Figs_paper/fig40_cb_case2   
  

psf = [psf1 pad(Im,pn).*pad(maskc,pn)];

    figure(1), clf
        imagesc(log10(psf), [-12 -8]); axis image, colormap(jet); niceaxes
        parms = fun_newaxes(14, 0, 0);  
        colorbar, parms = fun_newaxes(14, 0, 0);  
        %xlabel('Field-angle [\lambda/D]');
        %ylabel('Field-angle [\lambda/D]');
        title('vs IWA                  vs OWA');
        %text(-24, 23, 'BW = 5%', 'FontSize', 14, 'Color', 'w')
        %axis(25*[-1 1 -1 1])
   print -dpng ../Figs_paper/fig40_psf  
%return

return
%end % --------------------------------------------------------
%% 


close all

% figure position: [left bottom width height]

fig_size = [10 10 750 1200];
handles.master = figure('Color',[1 1 1],'Position',fig_size);
set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')

flag_dm = 0;

            handles.textbox1 = annotation('textbox', [0.05 0.87 0.95 0.1], ...
                'String', ['1500nm-WFC with DM1 + DM2: dm-drift = 1.5nm, Itr = ' sprintf('%i/%i', ii, ns)],'Fontsize',20,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');

    aa = ampi; rms0 = stat2d(aa); 
    if ii == 1, cx1 = rms0(1)*5*[-1 1]; end

            fst = 14; %--Font size for titles in the subplots

% OuterPosition: [left bottom width height]

            h_psf = subplot(2,2,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.46, [1 1]*0.45])
            imagesc(aa, cx_amp); niceaxes
            %hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            title(sprintf('\\Deltaamp: rms = %0.2g', rms0(1)),'FontSize',14);
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            %title(sprintf('amp: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

            aa = opdi; rms0 = stat2d(aa); 

            h_psf = subplot(2,2,3); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.02, [1 1]*0.45])
            imagesc(aa, cx_opd); niceaxes
            %hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('\\Deltaphase: rms = %0.1fpm', rms0(1)),'FontSize',fst);
            %title(sprintf('phase: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

            h_psf = subplot(2,2,2); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.46, [1 1]*0.45])
            %set(h_psf, 'OuterPosition', [0.35, 0.46, [1 1]*0.45])
            imagesc(xv,xv,log10(psfn0), [-11 -8]);  
            axis xy equal tight; colorbar(h_psf); colormap(jet);
            axis(12*[-1 1 -1 1])
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('nominal NI: Cb = %0.3g', cb0))
            
            h_psf = subplot(2,2,4); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.02, [1 1]*0.45])
            %set(h_psf, 'OuterPosition', [0.35, 0.02, [1 1]*0.45])
            imagesc(xv,xv,log10(psfn), [-11 -8]);  
            axis xy equal tight; colorbar(h_psf); colormap(jet);
            axis(12*[-1 1 -1 1])
            ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('post-wfc NI: Cb = %0.3g', cbi))

return
close all

% figure position: [left bottom width height]

fig_size = [10 10 750 1200];
handles.master = figure('Color',[1 1 1],'Position',fig_size);
set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')

flag_dm = 0;

if flag_dm

            handles.textbox1 = annotation('textbox', [0.1 0.87 0.8 0.1], ...
                'String', sprintf('Poked DM2 Act(%i,%i) with h = %0.1fnm [%i/%i]',ri,ci,dval, ii, ns),'Fontsize',24,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');
else
            handles.textbox1 = annotation('textbox', [0.05 0.87 0.95 0.1], ...
                'String', sprintf('Perturbed Seg(%i) with Tilt RMS = %ipm [%i/%i]',ii, round(dval*1e3), ii, ns),'Fontsize',24,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');
end

    amp = (amp_target - ampnom) / max(ampnom(:));
    aa = amp; rms0 = stat2d(aa); 
    if ii == 1, cx1 = rms0(1)*5*[-1 1]; end

            fst = 14; %--Font size for titles in the subplots

% OuterPosition: [left bottom width height]

            h_psf = subplot(2,2,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.46, [1 1]*0.45])
            imagesc(aa); niceaxes
            hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            title(sprintf('\\Deltaamp: rms = %0.2g, pv = %0.1g', rms0),'FontSize',14);
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            %title(sprintf('amp: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

    if flag_dm        
        aa = opd_target; rms0 = stat2d(aa); 
    else
        aa = opdz; rms0 = stat2d(aa); 
    end

    if ii == 1, cx2 = [min(aa(:)) max(aa(:))]; end

            h_psf = subplot(2,2,3); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.02, [1 1]*0.45])
            imagesc(aa); niceaxes
            hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis xy equal tight; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('\\Deltaphase: rms = %0.2f, pv = %0.2fnm', rms0),'FontSize',fst);
            %title(sprintf('phase: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

            h_psf = subplot(2,2,2); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.46, [1 1]*0.45])
            %set(h_psf, 'OuterPosition', [0.35, 0.46, [1 1]*0.45])
            imagesc(xv,xv,log10(psf0), [-11 -8]);  
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
            title(sprintf('distorted NI: Cb = %0.3g', cb1))

return
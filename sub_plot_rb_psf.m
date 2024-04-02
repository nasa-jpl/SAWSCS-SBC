% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

close all

% figure position: [left bottom width height]

%fig_size = [500 10 750 1200];
fig_size = [500 10 1700 1210];
handles.master = figure('Color',[1 1 1],'Position',fig_size);
set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')

            handles.textbox1 = annotation('textbox', [0.05 0.87 0.95 0.1], ...
                'String', [sprintf('DM Drift Error RMS = 1000pm: Case-%i', ii)],'Fontsize',20,...
                'HorizontalAlignment','center','LineStyle','none','Interpreter','latex');

            fst = 14; %--Font size for titles in the subplots

% OuterPosition: [left bottom width height]

            aa = opdi*1e3; rms0 = stat2d(aa);

            h_psf = subplot(2,2,1); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.46, [1 1]*0.45])
            imagesc(aa,rms0(1)*3*[-1 1]); axis image, niceaxes
            %hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis image; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            title(sprintf('before-control wfe: rms = %0.1fpm', rms0(1)),'FontSize',14);
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            %title(sprintf('amp: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');

            aa = opdx(:,:,end)*1e3; rms0 = stat2d(aa);

            h_psf = subplot(2,2,2); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.46, [1 1]*0.45])
            %set(h_psf, 'OuterPosition', [0.35, 0.46, [1 1]*0.45])
            imagesc(aa, rms0(1)*3*[-1 1]); axis image, niceaxes 
            axis image; colorbar(h_psf); colormap(jet);
            %axis(12*[-1 1 -1 1])
            %ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('after-control wfe: rms = %0.1fpm', rms0(1)),'FontSize',14);
            
            h_psf = subplot(2,2,3); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.05, 0.02, [1 1]*0.45])
            imagesc(xv,xv,log10(psfi), [-11 -6]); axis image, niceaxes
            %hold on, contour(maskopd*max(aa(:)), 1, 'w'), hold off
            axis image; colorbar(h_psf); colormap(jet);
        %     xlabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('pre-wfc NI: Cb = %0.3g', cbi),'FontSize',fst);
            %title(sprintf('phase: poke act(%i,%i) = 0.5nm', ri, ci),'Fontsize',fst);%,'Fontweight','Bold');
            axis(12*[-1 1 -1 1])

            h_psf = subplot(2,2,4); % Save the handle of the subplot
            set(h_psf, 'OuterPosition', [0.55, 0.02, [1 1]*0.45])
            %set(h_psf, 'OuterPosition', [0.35, 0.02, [1 1]*0.45])
            imagesc(xv,xv,log10(psfo), [-11 -6]); axis image, niceaxes
            axis image; colorbar(h_psf); colormap(jet);
            %axis(12*[-1 1 -1 1])
            %ylabel('$\lambda_0$/D','FontSize',16,'Interpreter','LaTeX'); 
            set(gca,'FontSize',14,'FontName','Times','FontWeight','Normal')
            title(sprintf('post-wfc NI: Cb = %0.3g', cbo),'FontSize',fst);
            axis(12*[-1 1 -1 1])

        %set(gca,'nextplot','replacechildren');
        %frame = getframe(gcf);
        %writeVideo(video,frame);
        %close(video); 

return

    
ii = 1;

mp.wfc.fn = sprintf('~esidick/Afalco/falco20200916/macos/dat_dm_err/wfc_dm1_dm2_gap08_30nm_bw20_1500nm_dm_err_1000pm_run%i_v2', ii);
    load(mp.wfc.fn, 'dm1x', 'dm2x', 'opdx', 'opdi', 'rms_opd', 'rms_amp', 'dmrms')

        video = VideoWriter('dat_dm_err/movie_opd_1000pm.avi'); %create the video object
        video.FrameRate = 1;
        open(video); %open the file for writing

        ns = size(opdx,3)

figure(2), clf

fs = 14;

for jj = 1:ns
    opd = opdx(:,:,jj) * 1e3;
    rms0 = stat2d(opd);
    cx = rms0(1)*3*[-1 1];

        imagesc(opd, 100*[-1 1]); axis xy image, colormap(jet); 
        parms = fun_newaxes(fs, 0, 0);  
        colorbar, parms = fun_newaxes(fs, 0, 0);  
        title([sprintf('iter = %i: rms = %0.1fpm',  jj, rms0(1)) ]), drawnow

        set(gca,'nextplot','replacechildren');
        frame = getframe(gcf);
        writeVideo(video,frame);
end % ii - loop
close(video); 

return
%% ----------------------------------------------------


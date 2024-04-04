%% EXAMPLE MATLAB RANDOM GENERATION

   addpath ~esidick/matlab_erkin/
   %addpath ~esidick/matlab_joegreen/
   addpath /proj/exospin/users/joegreen/matlab/
       
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

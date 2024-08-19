% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
%SHOW_OPD - Display an OPD.
%
%fignum = show_opd(opd, my_title, units, clim)

function fignum=show_opd(opd, my_title, units, clim)

%% Check input arguments
if nargin < 1
   return
end
if nargin < 2 
   my_title = [];
end
if nargin < 3 
   units=' ';
end
if nargin < 4
   clim = [min(opd(:)) max(opd(:))];
end

%% Display the OPD
colormap('jet')
if nnz(opd)
   imagesc(opd,clim), axis image %xy
else
   imagesc(opd), axis image %xy
end
title(my_title)
xlabel(sprintf('RMS = %0.3f, PV = %0.3f',stat2d(opd)))
set(get(colorbar,'ylabel'),'String',units,'Rotation',90);

drawnow

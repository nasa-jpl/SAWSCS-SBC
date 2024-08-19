% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
function [thestats] = stat2d (opd)

 indx=find(isnan(opd));
 opd(indx)=0;
 
 therms = rms2(nonzeros(opd));
 thepv  = max(opd(:))-min(opd(:));

 thestats = [therms thepv];
 

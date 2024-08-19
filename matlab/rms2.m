% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% function [root] = rms2 (vec)
%
% Computes the root mean squared of a vector
%
% written by Scott Basinger, JPL, Caltech
%
function [root]=rms2(vec)

if size(vec,1)==0,
   root = 0;
else
   root=sqrt(sum((vec-mean(vec)).^2)/size(vec,1));
end

% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% Utility to vectorize matrix data
%   [mat]=v2m(vec,indx);

	function [mat]=vec2mat(vec,indx);
	m=indx.size(1);
	n=indx.size(2);
	mat = full(sparse(indx.i,indx.j,vec,m,n));

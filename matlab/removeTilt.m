% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
	function phase=removeTilt(phase,mask);

% 	function phase=removeTilt(phase,mask);

% Pad to square array
		[m n]=size(phase);
		mx=max([m,n]);
		if m<mx,
		   ma=floor((n-m)/2);
		   mb=mx-m-ma;
		   phase=[zeros(ma,n);phase;zeros(mb,n)];
		   mask=[zeros(ma,n);mask;zeros(mb,n)];
		elseif n<mx,
		   na=floor((m-n)/2);
		   nb=mx-n-na;
		   phase=[zeros(m,na) phase zeros(m,nb)];
		   mask=[zeros(m,na) mask zeros(m,nb)];
		end;

% Remove tilt
		[X,Y]=meshgrid(-(mx-1)/2:(mx-1)/2);
		X=mask.*X;
		Y=mask.*Y;
		tiltx=X/norm(nonzeros(X));
		tilty=Y/norm(nonzeros(Y));

		phase=mask.*(phase-mean(nonzeros(phase)));

		fx=phase(:)'*tiltx(:);
		phase=mask.*(phase-fx*tiltx);
		fy=phase(:)'*tilty(:);
		phase=mask.*(phase-fy*tilty);

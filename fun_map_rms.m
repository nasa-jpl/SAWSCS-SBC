% Copyright 2018-2020 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%


function rmso = fun_map_rms(map, nlength)

    rmsx = [];
    for ii = 1:nlength
        mapi = map(:,:,ii);
        rms0 = stat2d(mapi); 
        rmsx = [rmsx; rms0];
    end

    rmso = rmsx(:,1); 

return
%end % -----------------------------------------------
%

function [EP,segAp,m_CntrlObs,m_Struts,m_segGaps] = fun_make_pupil(N,gap_size_m,diam_m,obs_m,strut_size_m,lambda_m,tilt_set)
 
    
% Original Source from Dave Redding, 28 May 2019
%
% Edited - J. Jewell to work with other code,
% and to provide aberrations if needed
%
% EDited - 2 July 2019 by J. Jewell to generate
% an entrance pupil structure with the discrete
% wavelength set for broadband optimization, and possibly
% with a set of tilted modes for finite stellar disc (to be added
% later)

% Note - without primary segment tip tilt errors, there is
% redundancy in the discrete wavelength EP...
% With tip/tilt errors, there will be a "lambda_0 / lambda"
% wavelength dependence...

% Easier for coding to just assume the EP is a discrete
% wavelength and mode structure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the central wavelength and on-axis EP mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
radius = N/2;
samp_m = diam_m / (N - 1);
gap = gap_size_m / samp_m;
go2 = gap/2;
rad_obs = obs_m / 2 / samp_m;
rad_mid = rad_obs + (radius - rad_obs)/2;

m = zeros(radius*2);
[xx,yy] = meshgrid(-radius:1:radius-1,-radius:1:radius-1);
r2 = sqrt(xx.^2 + yy.^2);

% Zero values outside the aperture
m(find(r2 <= radius)) = 1;

% Zero values in the radial gaps -- outer ring
for ith = 1:6
    th = ith/6 * pi;    % Define a rectangular gap
    cth = cos(th);
    sth = sin(th);
    pu = radius*[cth; sth];
    C = pu + go2*[-sth; cth];
    D = pu + go2*[sth; -cth];
    pl = -pu;
    A = pl + go2*[sth; -cth];
    B = pl + go2*[-sth; cth];
    AB = B-A;
    BC = C-B;
    ABdotAB = AB'*AB;
    BCdotBC = BC'*BC;
    
    for i = -radius:1:radius-1  % Set values in each gap to zero
        for j = -radius:1:radius-1
            P = [i; j];
            AP = P-A;
            BP = P-B;
            ABdotAP = AB'*AP;
            BCdotBP = BC'*BP;
            if (ABdotAP>=0) && (ABdotAP<=ABdotAB) && (BCdotBP>=0) && (BCdotBP<=BCdotBC)
                m(i+radius+1,j+radius+1) = 0;
            end
        end
    end
end
            
% figure(2), clf
% plot(A(1),A(2),'x',B(1),B(2),'o',...
% C(1),C(2),'x',D(1),D(2),'o'), axis equal

% Zero values in the ring gap
m(find(r2 < rad_mid + go2)) = 0;
m(find(r2 < rad_mid - go2)) = 1;

% Zero values in the radial gaps -- inner ring
for ith = 1:3
    th = ith/3 * pi;    % Define a rectangular gap
    cth = cos(th);
    sth = sin(th);
    pu = rad_mid*[cth; sth];
    C = pu + go2*[-sth; cth];
    D = pu + go2*[sth; -cth];
    pl = -pu;
    A = pl + go2*[sth; -cth];
    B = pl + go2*[-sth; cth];
    AB = B-A;
    BC = C-B;
    ABdotAB = AB'*AB;
    BCdotBC = BC'*BC;
    
    for i = -radius:1:radius-1  % Set values in each gap to zero
        for j = -radius:1:radius-1
            P = [i; j];
            AP = P-A;
            BP = P-B;
            ABdotAP = AB'*AP;
            BCdotBP = BC'*BP;
            if (ABdotAP>=0) && (ABdotAP<=ABdotAB) && (BCdotBP>=0) && (BCdotBP<=BCdotBC)
                m(i+radius+1,j+radius+1) = 0;
            end
        end
    end
end

m_segGaps = m;
m_segGaps(find(r2 < rad_obs)) = 1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add in the struts...
%
% At angles (deg) [30,120,270]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = zeros(radius*2);
[xx,yy] = meshgrid(-radius:1:radius-1,-radius:1:radius-1);
r2 = sqrt(xx.^2 + yy.^2);

% Zero values outside the aperture
m(find(r2 <= radius)) = 1;

m_CntrlObs = m;
m_CntrlObs(find(r2 < rad_obs)) = 0;

strut_angles = [0,120,240]*pi/180;
strut = strut_size_m / samp_m;
st2 = strut/2;

for ith = 1:size(strut_angles,2)
    th = strut_angles(1,ith);    % Define a rectangular gap
    cth = cos(th);
    sth = sin(th);
    %pu = rad_mid*[cth; sth];
    pu = radius*[cth; sth];
    C = pu + st2*[-sth; cth];
    D = pu + st2*[sth; -cth];
    %pl = -pu;
    %A = pl + st2*[sth; -cth];
    %B = pl + st2*[-sth; cth];
    A = st2*[sth; -cth];
    B = st2*[-sth; cth];
    AB = B-A;
    BC = C-B;
    ABdotAB = AB'*AB;
    BCdotBC = BC'*BC;
    
    for i = -radius:1:radius-1  % Set values in each gap to zero
        for j = -radius:1:radius-1
            %for i = 0.0:1:radius-1  % Set values in each gap to zero
            %for j = 0.0:1:radius-1
            P = [i; j];
            AP = P-A;
            BP = P-B;
            ABdotAP = AB'*AP;
            BCdotBP = BC'*BP;
            if (ABdotAP>=0) && (ABdotAP<=ABdotAB) && (BCdotBP>=0) && (BCdotBP<=BCdotBC)
                m(i+radius+1,j+radius+1) = 0;
            end
        end
    end
end

%m(find(r2 < rad_obs)) = 0;

m_Struts = m;
m_Struts(find(r2 < rad_obs)) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the EP(nmode,nlambda).field structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlam = size(lambda_m,2);
nmode = size(tilt_set,2);

centerWL_index = round(nlam/2);
centerWL_m = lambda_m(1,centerWL_index);

EP='';
for jmode=1:nmode
    for jlam=1:nlam
        EP(jmode,jlam).field = m_CntrlObs.*m_Struts.*m_segGaps;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Aperture structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ap0 = m_CntrlObs.*m_Struts.*m_segGaps;

[Xgrid,Ygrid] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);
[THETA,RHO]=cart2pol(Xgrid,Ygrid);

ap_outer_supp = zeros(N,N);
[JR]=RHO(ap0 == 1);
Rmax = max(max(JR));
ap_outer_supp(RHO <= Rmax)=1.0;

segAp = '';
segAp.test_aperture = ap0;
segAp.aperture = ap0;
segAp.ap_outer_supp = ap_outer_supp;

segAp.ap_inner_supp = ap_outer_supp;
segAp.crop_CP = [1:N];
segAp.ap_inner_supp_crop = ap_outer_supp;
%segAp.ap_LyotStop_crop = ap_LyotStop_crop;

%segAp.ap_LyotStop = ap_LyotStop;

%segAp.diam_outer_m = diam_outer_LB;
%segAp.diam_inner_m = diam_inner_LB;
%segAp.pcnt_seg_gap = 100*gap_size_m/diam_outer_LB;

%segAp.EP_stop = EP_stop;
%segAp.LPM_InnerR = 0.0;
%segAp.LPM_OuterR = LP_stop;
%segAp.Lyot_fractionInnerR = Lyot_fractionInnerR;


%hexAp.nrings = nrings;
%hexAp.seg_flat_side2side_m = seg_flat_side2side_m;
%hexAp.gap_size_m = gap_size_m;
%hexAp.seg_center2center = seg_center2center;
%hexAp.seg_diam_flat2flat = seg_diam_flat2flat;
%hexAp.center2vertex = seg_center2vertex;
%hexAp.pcnt_central_obscr = 0.0;
%segAp.apDiam = apDiam;
%hexAp.lambdaOverD = npix/apDiam;
segAp.lambda_m = lambda_m;
segAp.tilt_set = tilt_set;
%segAp.ZV = ZV;

segAp.gap_size_m = gap_size_m;
segAp.strut_size_m = strut_size_m;
segAp.primary_outer_diam_m = diam_m;
segAp.secondary_diam_m = obs_m;
%ap0 = prop_hex_aperture(bmi,nr,hr,hs);

%crop = npix/2-vcg.N/2+1:npix/2+vcg.N/2;
%ap_crop = ap0(crop,crop);

segAp
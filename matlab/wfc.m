% Copyright 2024, by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------
%
% The following loads in surface parallel actuator influence functions
% computer by a FEM.  There are 42 actuators.
% This script removes tip/tilt from the inf fcns, computes a simple control 
% matrix, then generates a random wavefront and corrects it using the
% mirror actuators
load influence-functions-20220309.mat

% Set Flags
Display = 0; useZernike = 1; zernike_fcn_available = 0;
randn('seed',211660);
    
% There are 43 influence functions in the matrix D, but  only 42 actuators.
% The first column is all zeros, so remove it

D1 = D(:,2:end);
clear D;
nacts = size(D1,2);

% Create a binary mask for visible part of the hexagonal mirror
maskz = v2m(D1(:,1),indx) ~= 0;

% Remove first three zernike terms from all influence functions (piston, tip, tilt)
% Put results back into matrix D
D = [];
for ii = 1:nacts
   if zernike_fcn_available
      [opdi, zk] = zernike_remove(v2m(D1(:,ii),indx), [1 2 3]);
   else
      opdi = removeTilt(v2m(D1(:,ii),indx), maskz);
   end
   wi = m2v(opdi, indx);
   D = [D wi(:)];
   if Display
      figure(1),show_opd(v2m(1e9*D(:,ii),indx),sprintf('Inf Fcn %d, ptt removed',ii),'nm'),drawnow
   end
end

% Compute simple control (or "Gain") matrix
G = pinv(D);

if useZernike
   if zernike_fcn_available
      zk = 533e-9*randn(1,15);
      zk(1:3) = 0;
      opdi = zernike_compose(maskz, zk);
      wi = m2v(opdi,indx);
   else
      load ZernikeWF % contains opdi & wi
   end
else
   pert = 2*randn(42,1);
   wi = D*pert; % Perturbed wavefront
end
figure(2),show_opd(1e9*v2m(wi,indx),'Random actuator disturbance','nm')

du = -G * wi; % Actuator commands to fix wavefront
wo = D * du + wi;  % corrected vector wavefront
opdo = v2m(wo, indx); % corrected matrix wavefront
figure(3),show_opd(removeTilt(1e9*opdo,maskz),'Corrected wavefront','nm')

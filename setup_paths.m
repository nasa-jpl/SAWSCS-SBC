% Copyright 2024 by the California Institute of Technology. ALL RIGHTS
% RESERVED. United States Government Sponsorship acknowledged. Any
% commercial use must be negotiated with the Office of Technology Transfer
% at the California Institute of Technology.
% -------------------------------------------------------------------------

codedir = '/proj/jwst2/Afalco/'; % This is where to set paths to read code from

datadir = '/proj/jwst2/Afalco/'; %'/home/esidick/Afalco/';
falcodir = [datadir 'falco20200916/']; % This is where to load and save data files from

if ~exist('falco_init_storage_arrays')
   addpath([codedir 'matlab_erkin/']);
   addpath(genpath([codedir 'falco20200916/'])); %savepath;
   addpath([codedir 'proper_v3.0.1_matlab_22aug17/']); %savepath;
   addpath([codedir 'zernike/']); % Need to point to your favorite zernike util directory
end

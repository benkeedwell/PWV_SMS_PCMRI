function [ksp, data] = ReadTWIX(data_name, data_type)
%-------------------------------------------------------------------------------%
% READ TWIX DATA FROM MEAS.DAT FILE USING MAPVBVD
% By Benjamin Keedwell (2026)
% MapVBVD by Philipp Ehses with modifications by Yang Ji
%-------------------------------------------------------------------------------%
% Read raw data from meas.dat file using MapVBVD
% This modified MapVBVD is designed for scans using asymmetric echo.
%-------------------------------------------------------------------------------%
% INPUT
% data_name - Path to meas.dat file
% data_type - 'img' (SMS PC-MRI) or 'acs' (calibration)
%-------------------------------------------------------------------------------%
% OUTPUT
% ksp - TWIX k-space data as MATLAB array
% data - Extracted TWIX data including headers
%-------------------------------------------------------------------------------%

%Read meas.dat file
data = mapVBVD_ae(data_name);

%Set data as the second nested structure (i.e., discard noise scan)
data = data{length(data)};

if strcmp(data_type, 'img')
    %Ignore unused segmentation dimension
    data.image.flagIgnoreSeg = true;

    %Remove oversampling (factor 2) in read direction
    data.image.flagRemoveOS  = true;

    %Extract k-space data
    ksp = data.image{''};

elseif strcmp(data_type, 'acs')
    %Ignore unused segmentation dimension
    data.refscan.flagIgnoreSeg = true;

    %Remove oversampling (factor 2) in read direction
    data.refscan.flagRemoveOS  = true;

    %Extract k-space data
    ksp = data.refscan{''};

    %Extract second dataset within 'ksp'
    ksp = ksp(:,:,:,:,2);
end


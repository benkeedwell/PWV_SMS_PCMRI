function msk = IdentifyWraps(phs, Niqr, fig)

%-------------------------------------------------------------------------------%
% IDENTIFY VOXELS WITH WRAPPED WAVEFORMS
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Identify wrapped waveforms as voxels with outlier differences in
% consecutive time points.
%-------------------------------------------------------------------------------%
% INPUTS
% phs - Phase Image (NFrq NPhs NCph)
% Niqr - Number of IQRs from first and third quartiles
% fig - Show figures (bool)
%-------------------------------------------------------------------------------%
% OUTPUTS
% msk - Mask of wrapped voxels (NFrq NPhs)
%-------------------------------------------------------------------------------%

% Difference in consecutive time points
phs_dff = diff(phs, 1, 3);                                 %(NFrq NPhs NCph-1)

% Interquartile range and quartiles (1st and 3rd)
[phs_dff_iqr, phs_dff_quar] = iqr(phs_dff, 3);             %(NFrq NPhs) (NFrq NPhs 2)

% Calculate upper and lower fences
phs_lf = phs_dff_quar(:,:,1) - Niqr*phs_dff_iqr;           %(NFrq NPhs)
phs_uf = phs_dff_quar(:,:,2) + Niqr*phs_dff_iqr;           %(NFrq NPhs)

% Identify wrapped waveforms
msk = any(phs_dff > phs_uf, 3) | any(phs_dff < phs_lf, 3); %(NFrq NPhs)

% Plot figure
if fig==1
    f = figure(color = 'white');
    imshow(msk)
    title(f.Children, 'Wrapped Voxels Mask')
end
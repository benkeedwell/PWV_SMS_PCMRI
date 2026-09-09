function VoxelWaveforms(phs, rec, msk, msk_sub, phslim, clr, label)

%-------------------------------------------------------------------------------%
% DISPLAY GRID OF INDIVIDUAL VOXEL WAVEFORMS
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Display grid of waveforms for input region of image. Masks can also be
% input to colour-code waveforms.
%-------------------------------------------------------------------------------%
% INPUTS
% phs - Phase image (NFrq NPhs NCph)
% rec - Rectangular region [x1, y1, x2, y2]
% msk - Mask [] OR (NFrq NPhs)
% msk_sub - Subset of mask [] or (NFrq NPhs)
% phs_lim - Phase axis limits [upper, lower]
% clr - Waveform colour [R, G, B]
% label - Figure title
%-------------------------------------------------------------------------------%

% Select image region
phs_rec = phs(rec(2):rec(4), rec(1):rec(3), :);
[NFrq, NPhs, NCph] = size(phs_rec);

% If requested, select mask and sub-mask region
if ~isempty(msk)
    msk_rec = msk(rec(2):rec(4), rec(1):rec(3));
    if ~isempty(msk)
        msk_sub_rec = msk_sub(rec(2):rec(4), rec(1):rec(3));
        msk_rec = msk_rec - msk_sub_rec;
    else
        msk_sub_rec = zeros(NFrq, NPhs);
    end
else
    msk_rec = zeros(NFrq, NPhs);
    msk_sub_rec = zeros(NFrq, NPhs);
end

% Plot figure
figure('Units', 'normalized', 'Position', [0,0,1,1], 'Color', 'white');
t = tiledlayout(NFrq, NPhs, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, label)
for i = 1:NFrq
    for j = 1:NPhs
        nexttile
        wf_ij = squeeze(phs_rec(i,j,:));
        plot(wf_ij, 'Color', clr, 'LineWidth', 1);
        xlim([0, NCph+1])
        ylim(phslim)
        if msk_rec(i,j) == 1
            clr_ax = [1,0.95,0.95];
        elseif msk_sub_rec(i,j) == 1
            clr_ax = [0.95,1,0.95];
        else
            clr_ax = [1,1,1];
        end
        if msk_rec(i,j) == 1 || msk_sub_rec(i,j) == 1
            %pat_ij = PATFilterTangents(wf_ij.', 4, 10, 100, 0.2, 0, '', [], 1);
            pat_ij = PATFilterTangents(wf_ij.', 1, 20, 80, 30, 1, 0);
            text(3, phslim(2)-0.5, num2str(pat_ij, 3), 'FontSize', 12);
        end
        set(gca,'Color', clr_ax, 'LineWidth', 1)
        if (i == NFrq && j == 1)
            xlabel('Cardiac Phase')
            ylabel('Phase')
        else
            set(gca,'XTick',[], 'YTick', [])
        end
        pbaspect([1,1,1])
    end
end
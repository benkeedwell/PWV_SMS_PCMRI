function msk = SmoothnessMask(phs, mode, fig)
%-------------------------------------------------------------------------------%
% GENERATE MASK THRESHOLDED BY SMOOTHNESS
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Compute phase-shift smoothness and plot histogram of values. Either
% select threshold using this histogram or override.
%-------------------------------------------------------------------------------%
% INPUTS
% phs - Phase-shift image (NFrq NPhs NCph)
% mode - 'histogram'/thr controls what threshold is selected
% fig - Show figures (boolean)
%-------------------------------------------------------------------------------%
% OUTPUTS
% msk - Smoothness mask (NFrq NPhs)
%-------------------------------------------------------------------------------%

% Generate phase smoothness map
phs_smo = 1 - (var(diff(phs, 1, 3), [] ,3) ./ var(phs, [], 3));  %(NFrq NPhs)

% Estimate phase smoothness threshold
[h_cnt, h_edg] = histcounts(phs_smo(:)); %(1 NBin) (1 NBin+1)
h_bin = h_edg(1:end-1) + (h_edg(2) - h_edg(1))/2; %(1 NBin)
wth = 1.5;
[~, max_ind] = max(h_cnt);
bin_max = h_bin(max_ind);
bin_bool = (h_bin >= bin_max - wth/2) & (h_bin <= bin_max + wth/2);
gauss = fit(h_bin(bin_bool).', h_cnt(bin_bool).', 'gauss1');

if strcmp(mode, 'histogram')
    thr = gauss.b1 + 3*gauss.c1;
else
    thr = mode;
end

% Generate phase smoothness mask
msk = phs_smo > thr; %(NFrq NPhs)

% Plot figures
if fig==1
    figure(color = 'white')
    t = tiledlayout('TileSpacing', 'compact');
    title(t, 'PHASE SMOOTHNESS')
    nexttile
    ImgScale(phs_smo, 2, 98)
    title('Phase Smoothness Map')
    nexttile
    histogram(phs_smo(:))
    hold on
    plot(gauss)
    xline(thr)
    xline([bin_max-wth/2, bin_max+wth/2], '--')
    xlabel('Phase Smoothness')
    ylabel('Count')
    legend('Histogram', 'Gaussian', 'Threshold')
    title('Phase Smoothness Histogram')
    nexttile
    ImgScale(msk, 0, 100)
    title('Phase Smoothness Mask')
end
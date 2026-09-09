function msk = MaxPhaseMask(phs, bidirec, mode, fig)
%-------------------------------------------------------------------------------%
% GENERATE MASK THRESHOLDED BY MAXIMUM PHASE
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Compute maximum phase-shift and plot histogram of values. Either
% select threshold using this histogram or override.
%-------------------------------------------------------------------------------%
% INPUTS
% phs - Phase-shift image (NFrq NPhs NCph)
% bidirec - Include maximum velocities in both directions (boolean)
% mode - 'static'/'vessel'/thr controls what threshold is selected
% fig - Show figures (boolean)
%-------------------------------------------------------------------------------%
% OUTPUTS
% msk - Maximum phase mask (NFrq NPhs)
%-------------------------------------------------------------------------------%

% Generate maximum phase map
if bidirec == 1
    phs_max = max(abs(phs), [], 3); %(NFrq NPhs)
elseif bidirec == 0
    phs_max = max(phs, [], 3); %(NFrq NPhs)
end

% Estimate maximum phase threshold
[h_cnt, h_edg] = histcounts(phs_max(:)); %(1 NBin) (1 NBin+1)
h_bin = h_edg(1:end-1) + (h_edg(2) - h_edg(1))/2; %(1 NBin)
start = [max(h_cnt(h_bin < pi/2)),...                          % a1
         h_bin(find(h_cnt == max(h_cnt(h_bin < pi/2)), 1)),... % b1
         0.2,...                                               % c1
         max(h_cnt(h_bin > pi/2)),...                          % a2
         h_bin(find(h_cnt == max(h_cnt(h_bin > pi/2)), 1)),... % b2
         0.2];                                                 % c2
gauss = fit(h_bin.', h_cnt.', 'gauss2', 'StartPoint', start);
coeff = reshape(coeffvalues(gauss), 3, 2); %(3 2)
[~, order] = sort(coeff(2,:));
coeff = coeff(:,order); %(3 2)


if strcmp(mode, 'static')
    thr = coeff(2,1) + 3*coeff(3,1);
elseif strcmp(mode, 'vessel')
    thr = (coeff(2,2) + coeff(2,1))/2;
else
    thr = mode;
end


% Generate maximum phase mask
msk = phs_max > thr; %(NFrq NPhs)

% Plot figures
if fig==1
    figure(color = 'white')
    t = tiledlayout('TileSpacing', 'compact');
    title(t, 'MAXIMUM PHASE')
    nexttile
    ImgScale(phs_max, 1, 99)
    title('Maximum Phase Map')
    nexttile
    histogram(phs_max(:))
    hold on
    plot(gauss)
    xline(thr)
    xlabel('Maximum Phase')
    ylabel('Count')
    legend('Histogram', 'Double Gaussian', 'Threshold')
    title(['Maximum Phase Histogram (Mode: ' mode ')'])
    nexttile
    ImgScale(msk, 0, 100)
    title('Maximum Phase Mask')
end
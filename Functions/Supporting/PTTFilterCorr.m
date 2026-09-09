function ptt = PTTFilterCorr(wf, pp, lb, ub, type, fig)

%-------------------------------------------------------------------------------%
% FILTER WAVEFORMS & ESTIMATE PTT USING SLIDING CORRELATION
% By Benjamin Keedwell (2026)
% Filter by M Schmid, D Rath and U Diebold
% PCC method based on a similar implementation by E Peper, B Coolen et al
%-------------------------------------------------------------------------------%
% Filter waveforms using a modified sinc kernel and estimate pulse transit
% time through cross-correlation, mean squared error or Pearson's correlation 
% coefficient (PCC).
%-------------------------------------------------------------------------------%
% INPUTS
% wf - Waveforms (NWav NCph)
% pp - Preprocessing: 0 = Interpolation, 1 = Filtering
% lb - Lower bound (%) of upstroke
% ub - Upper bound (%) of upstroke
% type - 'XCorr', 'MSE' or 'Pearsons'
% fig - No figures (0), basic figure (1) or Technical Note figure (2)
%-------------------------------------------------------------------------------%
% OUTPUTS
% ptt - Pulse transit time relative to first waveform (1 NWav)
%-------------------------------------------------------------------------------%

% Parameters
n_MS = 4;              % MS filter order
framehalf_MS = 10;     % MS filter frame halflength
interp_factor = 100;   % Interpolation factor
ExWav = 2;             % Example waveforms for figure (e.g. 2:3)

% Setup
[NWav, NCph] = size(wf);
ptt = zeros(1, NWav);
lb = lb/100;
ub = ub/100;
if NWav < 2
    error('At least two waveforms required to compute PTT')
end
if ~any(strcmp(type, {'XCorr', 'MSE', 'Pearsons'}))
    error('type must be XCorr, MSE or Pearsons')
end
if fig == 1
    figure(Color = 'white');
    tiledlayout(length(ExWav), 2, 'TileSpacing', 'tight');
end

% Interpolate input time series
x = 1:NCph;                                   %(1 NCph)
x_int = 1:(1/interp_factor):NCph;             %(1 NCph_int)
y_int = interp1(x, wf.', x_int, 'spline').';  %(NWav NCph_int)

for i = 1:NWav

    % If requested, filter interpolated waveform
    if pp == 1
        y = smoothMS(y_int(i,:), n_MS, framehalf_MS * interp_factor).'; %(NCph_int 1)
    elseif pp == 0
        y = y_int(i,:).'; %(NCph_int 1)
    end

    % Normalise waveform
    [maxy, maxx] = max(y);
    [miny, ~] = min(y(1:maxx));
    y = (y - miny)./(maxy - miny);

    % Select template region
    if i == 1
        
        % Isolate upstroke
        yup = y(1:maxx);
        xup = (1:maxx).';

        % Limit upstroke
        tem_bool = (yup >= lb) & (yup <= ub);
        grp_start = find(diff([0; tem_bool(:); 0]) == 1);
        grp_end = find(diff([0; tem_bool(:); 0]) == -1) - 1;
        [~, grp_idx] = max(grp_end - grp_start + 1);
        tem_bool(:) = false;
        tem_bool(grp_start(grp_idx):grp_end(grp_idx)) = true;
        ytem = yup(tem_bool);
        xtem = xup(tem_bool);
        tem_bool = logical([tem_bool; zeros(length(y) - length(yup), 1)]);

    % Compare to target waveform
    else

        % Compute selected correlation metric
        ytar = y;
        xtar = 1:length(ytar);
        shift = (xtem(end) - xtar(end)):(xtem(1) - xtar(1));
        r = zeros(length(shift), 1);
        for j = 1:length(shift)
            ytar_shift = circshift(ytar, shift(j));
            if strcmp(type, 'XCorr')
                r(j) = sum(ytem .* ytar_shift(tem_bool));
            elseif strcmp(type, 'MSE')
                r(j) = mean((ytem - ytar_shift(tem_bool)).^2);
            elseif strcmp(type, 'Pearsons')
                r(j) = corr(ytem, ytar_shift(tem_bool));
            end
        end
        
        % Select PTT
        if strcmp(type, 'XCorr')
            [~, shift_idx] = max(r);
        elseif strcmp(type, 'MSE')
            [~, shift_idx] = min(r);
        elseif strcmp(type, 'Pearsons')
            [~, shift_idx] = max(r);
        end
        shift_max = shift(shift_idx);
        ptt(i) = -shift_max / interp_factor;

    end

    % Figures

    if (i > 1) && any(i == ExWav)
        if fig == 1

            nexttile
            if pp == 1
                title(['Waveform Filtering (Waveform ' num2str(i) ')'])
            elseif pp == 0
                title(['Waveform Interpolation (Waveform ' num2str(i) ')'])
            end
            hold on
            plot(x,  wf(i,:), '*', 'MarkerSize', 3, 'Color', 'k');
            plot(x_int, y * (maxy - miny) + miny, 'r');
            if pp == 1
                legend({'Input', 'Filtered'})
            elseif pp == 0
                legend({'Input', 'Interpolated'})
            end
            xlabel('Cardiac Phase')
            ylabel('Waveform')
            box on;

            nexttile
            title(['Sliding ' type ' (Waveform ' num2str(i) ')'])
            hold on
            plot(xup, yup, 'k-', 'LineWidth', 1)
            plot(xtem, ytem, 'k-', 'LineWidth', 3)
            plot(xtar, ytar, 'r-', 'LineWidth', 1)
            plot(xtar(tem_bool) - shift_max, ytar(circshift(tem_bool, -shift_max)), 'r-', 'LineWidth', 3)
            plot(xtar(tem_bool), ytar(circshift(tem_bool, -shift_max)), 'r:', 'LineWidth', 3)
            xlabel('Interpolated Cardiac Phase')
            ylabel('Normalised Waveform')
            ylim([0,1])
            legend({'', 'Slice-1', '', ['Slice-' num2str(i)], ['Slice-' num2str(i) '-Shifted']})
            box on

        end

        if fig == 2

            TR = 6.38;
            TE = 3.92;
            xup_cphs = 1 + (xup - 1)/interp_factor;
            xup_ms = (xup_cphs * TR) + TE;
            xtem_cphs = 1 + (xtem - 1)/interp_factor;
            xtem_ms = (xtem_cphs * TR) + TE;
            xtar_cphs = 1 + (xtar - 1)/interp_factor;
            xtar_ms = (xtar_cphs * TR) + TE;
            shift_max_ms = shift_max / interp_factor * TR;

            figure;
            set(gcf, 'WindowState', 'maximized');
            hold on
            plot(xup_ms, yup, 'k-', 'LineWidth', 2)
            plot(xtem_ms, ytem, 'k-', 'LineWidth', 6)
            plot(xtar_ms, ytar, 'm-', 'LineWidth', 2)
            plot(xtar_ms(tem_bool) - shift_max_ms, ytar(circshift(tem_bool, -shift_max)), 'm-', 'LineWidth', 6)
            plot(xtar_ms(tem_bool), ytar(circshift(tem_bool, -shift_max)), 'm:', 'LineWidth', 6)
            xlabel('Time (ms)')
            ylabel('Normalised Waveform')
            ylim([0,1])
            xlim([0, xup_ms(end)])
            xticks(0:20:xup_ms(end))
            %legend({'', 'Slice 1', '', ['Slice ' num2str(i)], ['Slice ' num2str(i) ' Shifted']}, 'Location', 'northwest');
            legend({'', 'Template', '', 'Target', 'Target (Shifted)'}, 'Location', 'northwest');
            legend('boxoff')
            pbaspect([1.5,1,1])
            fontsize(40, 'points')
            ax = gca;
            ax.LineWidth = 4.0;
            box on

        end
    end
end


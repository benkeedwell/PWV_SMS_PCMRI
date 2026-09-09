function pat = PATFilterTangents(wf, pp, lb, ub, win, bl, fig)

%-------------------------------------------------------------------------------%
% FILTER WAVEFORMS & ESTIMATE PAT BY TANGENT INTERSECTION
% By Benjamin Keedwell (2026)
% Filter by M Schmid, D Rath and U Diebold
%-------------------------------------------------------------------------------%
% Filter waveforms using a modified sinc kernel and estimate pulse arrival
% time as the intersection of tangents fitted to the baseline and upstroke.
%-------------------------------------------------------------------------------%
% INPUTS
% wf - Waveforms (NWav NCph)
% pp - Preprocessing: 0 = Interpolation, 1 = Filtering
% lb - Lower bound (%) of upstroke to fit
% ub - Upper bound (%) of upstroke to fit
% win - Sliding window (%) of upstroke between bounds to fit ([] to ignore this step)
% bl - Baseline zero (0), median (1) or minimum (2)
% fig - No figures (0), basic figure (1) or Technical Note figure (2 or 3)
%-------------------------------------------------------------------------------%
% OUTPUTS
% pat - Pulse arrival time for each waveform (1 NWav)
%-------------------------------------------------------------------------------%

% Parameters
n_MS = 4;              % MS filter order
framehalf_MS = 10;     % MS filter frame halflength
interp_factor = 100;   % Interpolation factor
ExWav = 1;             % Example waveforms for figure (e.g 2:3)
NWin = 100;            % Number of windows to trial

% Setup
[NWav, NCph] = size(wf);
pat = zeros(1, NWav);
lb = lb/100;
ub = ub/100;
win = win/100;
if fig == 1
    figure(Color = 'white');
    tiledlayout(length(ExWav), 2, 'TileSpacing', 'tight');
end

% Interpolate input time series
x = 1:NCph;                                       %(1 NCph)
x_int = 1:(1/interp_factor):NCph;                 %(1 NCph_int)
if NWav>1
    y_int = interp1(x, wf.', x_int, 'spline').';  %(NWav NCph_int)
elseif NWav==1
    y_int = interp1(x, wf.', x_int, 'spline');    %(NWav NCph_int)
end

for i = 1:NWav

    % If requested, filter interpolated waveform
    if pp == 1
        y = smoothMS(y_int(i,:), n_MS, framehalf_MS * interp_factor).'; %(NCph_int 1)
    elseif pp == 0
        y = y_int(i,:).'; %(NCph_int 1)
    end

    % Isolate upstroke
    [maxy, maxx] = max(y);
    [miny, ~] = min(y(1:maxx));
    yup = y(1:maxx);
    xup = (1:maxx).';
  
    % Select region to fit tangent
    ylb = miny + lb * (maxy-miny);
    yub = miny + ub * (maxy-miny);
    upline_bool = (yup >= ylb) & (yup <= yub);
    yupline = yup(upline_bool);
    xupline = xup(upline_bool);

    % Calculate most linear window
    if ~isempty(win)
        ywin = win * (maxy-miny);
        starts = linspace(ylb, yub - ywin, NWin);
        r2 = zeros(NWin, 1);
        for j = 1:NWin
            bool = yupline >= starts(j) & yupline <=  (starts(j) + ywin);
            ytrial = yupline(bool);
            xtrial = xupline(bool);
            Sx  = sum(xtrial);
            Sy  = sum(ytrial);
            Sxx = sum(xtrial.^2);
            Syy = sum(ytrial.^2);
            Sxy = sum(xtrial.*ytrial);
            L = length(ytrial);
            num = (L .* Sxy - Sx .* Sy).^2;
            den = (L .* Sxx - Sx.^2) .* (L .* Syy - Sy.^2);
            r2(j) = num ./ den;
            if r2(j) == max(r2)
                ylinear = ytrial;
                xlinear = xtrial;
            end
        end
    else
        ylinear = yupline;
        xlinear = xupline;
    end

    % Fit tangent
    upfit = fit(xlinear, double(ylinear), 'poly1');

    % Calculate intersection
    if any(bl == [1,2])
        x_tangent = (miny - upfit.p2) / upfit.p1;
        y_tangent = miny;

        if bl == 1
            % Recalculate baseline using median of points before intersection
            baseline_median = median(y(1:floor(x_tangent)));
            x_tangent = (baseline_median - upfit.p2) / upfit.p1;
            y_tangent = baseline_median;
        end

    elseif bl == 0
        x_tangent = (0 - upfit.p2) / upfit.p1;
        y_tangent = 0;

    end

    % Convert PAT to uninterpolated cardiac phases
    pat(i) = 1 + (x_tangent - 1)/interp_factor;
    if (pat(i) < 1) || (pat(i) > NCph)
        pat(i) = nan;
    end

    % Figures
    if any(i == ExWav) || (i==5 && fig==3)
        if fig == 1

            nexttile
            if pp == 1
                title(['Waveform Filtering (Waveform ' num2str(i) ')'])
            elseif pp == 0
                title(['Waveform Interpolation (Waveform ' num2str(i) ')'])
            end
            hold on
            plot(x,  wf(i,:), '*', 'MarkerSize', 3, 'Color', 'k');
            plot(x_int, y, 'r');
            xline(pat(i), 'k', 'PAT')
            if pp == 1
                legend({'Input', 'Filtered'})
            elseif pp == 0
                legend({'Input', 'Interpolated'})
            end
            xlabel('Cardiac Phase')
            ylabel('Waveform')
            box on;

            nexttile
            title(['Tangent Intersection (Waveform ' num2str(i) ')'])
            hold on
            plot(1:length(y), y, 'k', 'LineWidth', 0.5)
            plot(xup, yup, 'r', 'LineWidth', 1)
            if ~isempty(win)
                plot(xupline, yupline, 'b', 'LineWidth', 2)
            end
            plot(xlinear, ylinear, 'm', 'LineWidth', 3)
            plot(1:xupline(end), upfit(1:xupline(end)), 'k', 'LineWidth', 1)
            yline(y_tangent, 'k', 'LineWidth', 1)
            scatter(x_tangent, y_tangent, 50, 'k')
            if ~isempty(win)
                legend({'Filtered', 'Upstroke', 'Initial Fit Region', 'Final Fit Region'});
            else
                legend({'Filtered', 'Upstroke', 'Fit Region'});
            end
            xlabel('Interpolated Cardiac Phase')
            ylabel('Waveform')
            ylim([-0.1*maxy, 1.1*maxy])
            box on;

        end

        if fig == 2 || (i==5 && fig==3)

            TR = 6.38;
            TE = 3.92;
            VENC = 50;
            xup_cphs = 1 + (xup - 1)/interp_factor;
            xup_ms = (xup_cphs * TR) + TE;
            xupline_cphs = 1 + (xupline - 1)/interp_factor;
            xupline_ms = (xupline_cphs * TR) + TE;
            xlinear_cphs = 1 + (xlinear - 1)/interp_factor;
            xlinear_ms = (xlinear_cphs * TR) + TE;
            xupfit = 1:xupline(end);
            xupfit_cphs = 1 + (xupfit - 1)/interp_factor;
            xupfit_ms = (xupfit_cphs * TR) + TE;
            x_tangent_cphs = 1 + (x_tangent - 1)/interp_factor;
            x_tangent_ms = (x_tangent_cphs * TR) + TE;

            figure;
            set(gcf, 'WindowState', 'maximized');
            hold on
            plot(xup_ms, yup/pi*VENC, 'k', 'LineWidth', 6)
            if ~isempty(win)
                plot(xupline_ms, yupline/pi*VENC, 'm', 'LineWidth', 9)
                plot(xlinear_ms, ylinear/pi*VENC, 'b', 'LineWidth', 12)
                plot(xupfit_ms, upfit(1:xupline(end))/pi*VENC, 'k', 'LineWidth', 2)
            else
                plot(xlinear_ms, ylinear/pi*VENC, 'm', 'LineWidth', 9)
                plot(xupfit_ms, upfit(1:xupline(end))/pi*VENC, 'k', 'LineWidth', 2)
            end
            yline(y_tangent/pi*VENC, 'k', 'LineWidth', 2, 'Alpha', 1)
            xline(x_tangent_ms, 'k:', 'LineWidth', 4)
            if ~isempty(win)
                legend({'Pulse Wave', 'Initial Fit Region', 'Final Fit Region', '', '', 'Intersection'}, 'Location', 'northwest');
            else
                legend({'Pulse Wave', 'Fit Region', '', '', 'Intersection'}, 'Location', 'northwest');
            end
            legend('boxoff')
            xlabel('Time (ms)')
            ylabel('Velocity (cm/s)')
            ylim([10,70])
            xlim([0, xup_ms(end)])
            xticks(0:20:xup_ms(end))
            pbaspect([1.5,1,1])
            fontsize(40, 'points')
            ax = gca;
            ax.LineWidth = 4.0;
            box on;

        end
    end
end


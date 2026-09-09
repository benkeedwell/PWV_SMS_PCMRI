function pat = PATFilterUpstroke(wf, pp, fig)

%-------------------------------------------------------------------------------%
% FILTER WAVEFORMS & ESTIMATE PAT AS MAXIMUM FIRST DERIVATIVE
% By Benjamin Keedwell (2026)
% Filter by M Schmid, D Rath and U Diebold
%-------------------------------------------------------------------------------%
% Filter waveforms using a modified sinc kernel and estimate pulse arrival
% time as the point of maximum first derivative.
%-------------------------------------------------------------------------------%
% INPUTS
% wf - Waveforms (NWav NCph)
% pp - Preprocessing: 0 = Interpolation, 1 = Filtering
% fig - No figures (0), basic figure (1) or Technical Note figure (2)
%-------------------------------------------------------------------------------%
% OUTPUTS
% pat - Pulse arrival time for each waveform (1 NWav)
%-------------------------------------------------------------------------------%

% Parameters
n_MS = 4;              % MS filter order
framehalf_MS = 10;     % MS filter frame halflength
interp_factor = 100;   % Interpolation factor
ExWav = 1;             % Example waveforms for figure (e.g 2:3)

% Setup
[NWav, NCph] = size(wf);
pat = zeros(1, NWav);
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
    [~, maxx] = max(y);
    yup = y(1:maxx);
    xup = (1:maxx).';

    % Differentiate
    dydx = diff(yup);
    [~, x_max_dydx] = max(dydx);
  
    % Convert PAT to uninterpolated cardiac phases
    pat(i) = 1 + (x_max_dydx - 1)/interp_factor;

    % Figures
    if any(i == ExWav)
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
            title(['First Derivative (Waveform ' num2str(i) ')'])
            plot(xup, yup ./ max(yup), 'k', 'LineWidth', 0.5)
            hold on
            plot(xup(1:end-1), dydx ./ max(dydx), 'r', 'LineWidth', 0.5)
            xline(x_max_dydx, 'k', 'PAT')
            xlabel('Interpolated Cardiac Phase')
            ylabel('Normalised Waveform')
            legend({'Upstroke', 'Derivative'})
            box on;

        end

        if fig == 2

            TR = 6.38;
            TE = 3.92;
            xup_cphs = 1 + (xup - 1)/interp_factor;
            xup_ms = (xup_cphs * TR) + TE;
            x_max_dydx_cphs = 1 + (x_max_dydx - 1)/interp_factor;
            x_max_dydx_ms = (x_max_dydx_cphs * TR) + TE;

            figure;
            set(gcf, 'WindowState', 'maximized');
            plot(xup_ms, rescale(yup), 'k', 'LineWidth', 6)
            hold on
            plot(xup_ms(1:end-1), rescale(dydx), 'm', 'LineWidth', 6)
            xline(x_max_dydx_ms, 'm:', 'LineWidth', 4)
            xlabel('Time (ms)')
            ylabel('Normalised Waveform')
            ylim([0,1])
            xlim([0, xup_ms(end)])
            xticks(0:20:xup_ms(end))
            legend({'Pulse Wave', 'First Derivative', 'Max. First Derivative'}, 'Location', 'northwest')
            legend('boxoff')
            pbaspect([1.5,1,1])
            fontsize(40, 'points')
            ax = gca;
            ax.LineWidth = 4.0;
            box on;

        end
    end
end


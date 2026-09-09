%-------------------------------------------------------------------------------%
% GENERATE FIGURES
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Generate requested figure for paper.
%-------------------------------------------------------------------------------%
% INPUTS
% fig - Requested figure ('TOF', 'Waveforms', 'PTT', 'Perc', 'TTF') (string)
%-------------------------------------------------------------------------------%

function GenerateFigures(fig)

% Figure 1a: TOF MIP
if strcmp(fig, 'TOF')

    addpath(genpath('/Users/benjaminkeedwell/Data_Storage/Physiological_Variation_Study/20260429_J/TOF'))

    sag = dicomread('MR.1.3.12.2.1107.5.2.43.66050.2026042911041395934504424');
    cor = dicomread('MR.1.3.12.2.1107.5.2.43.66050.2026042911041396229104425');

    sli = [33, 80, 127, 175, 222];

    clr = [0.0902 0.4667 0.7019  % Blue
        0.9608 0.4980 0.1373  % Orange
        0.2706 0.6275 0.3333  % Green
        0.9412 0.2314 0.1255  % Red
        0.6000 0.3020 0.6392  % Purple
        ];

    figure
    set(gcf, 'WindowState', 'maximized');
    imshow(rot90(sag, 2), [])
    for i = 1:5
        yline(sli(i)+1, 'Color', clr(i,:), 'LineWidth', 5, 'Alpha', 1)
    end

    figure
    set(gcf, 'WindowState', 'maximized');
    imshow(rot90(cor, 2), [])
    for i = 1:5
        yline(sli(i)+1, 'Color', clr(i,:), 'LineWidth', 5, 'Alpha', 1)
    end

end

% Figure 1d: Waveforms
if strcmp(fig, 'Waveforms')

    load('20260429_J/TS_J_1_LICA.mat', 'ts_phs', 'TR', 'TE')

    VENC = 50;

    [NSli, NCph] = size(ts_phs);

    ts_vel = ts_phs ./ pi * VENC;

    t = ((1:NCph) * TR) + TE;

    clr = [0.0902 0.4667 0.7019  % Blue
           0.9608 0.4980 0.1373  % Orange
           0.2706 0.6275 0.3333  % Green
           0.9412 0.2314 0.1255  % Red
           0.6000 0.3020 0.6392  % Purple
          ];

    figure
    set(gcf, 'WindowState', 'maximized');
    legpd = cell(NSli, 1);
    hold on
    for i = 1:NSli
        plot(t, ts_vel(i,:), 'Color', clr(i,:), 'LineStyle', '-', 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 18);
        legpd{i} = append('Slice ', string(i));
    end
    legend(legpd, 'Location', 'northwest');
    legend('boxoff')
    xlabel('Time (ms)');
    ylabel('Velocity (cm/s)');
    xlim([0, 220])
    pbaspect([1.1,1,1])
    fontsize(36, 'points')
    ax = gca;
    ax.LineWidth = 4.0;
    box on

end

% Figure 2: PTT Methods
if strcmp(fig, 'PTT')

    load('20260429_J/TS_J_1_LICA.mat', 'ts_phs')

    PTTFilterCorr(ts_phs, 1, 0, 80, 'Pearsons', 2);

    PATFilterUpstroke(ts_phs, 1, 2);

    PATFilterTangents(ts_phs, 1, 20, 80, [], 1, 2);

    PATFilterTangents(ts_phs, 1, 20, 80, 30, 1, 2);

end

% Figure S3: TTF Methods Example
if strcmp(fig, 'TTF')

    load('20260317_B/TS_B_3_LICA.mat', 'ts_phs')

    PATFilterTangents(ts_phs, 1, 20, 80, [], 1, 3);

    PATFilterTangents(ts_phs, 1, 20, 80, 30, 1, 3);

end

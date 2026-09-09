%-------------------------------------------------------------------------------%
% COMPARE PULSE WAVE VELOCITY ESTIMATES ACROSS PARTICIPANTS
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Load PWV across participants and compare across methods
%-------------------------------------------------------------------------------%
% INPUTS
% direc - PWV directories (cell array of strings)
% age_file - Age data file name (string)
% pl_file - Pathlength data file name (string)
%-------------------------------------------------------------------------------%

function PWV4_CompareParticipants(direc, age_file, pl_file)

% SETUP
vessel = {'LICA', 'RICA'};

clr_par = [                   % (NPar, 3)
    0.0000, 0.4470, 0.7410;   % Blue
    0.8500, 0.3250, 0.0980;   % Orange
    0.9290, 0.6940, 0.1250;   % Yellow
    0.4940, 0.1840, 0.5560;   % Purple
    0.4660, 0.6740, 0.1880;   % Green
    0.3010, 0.7450, 0.9330;   % Cyan
    0.6350, 0.0780, 0.1840;   % Dark red
    0.0000, 0.0000, 0.0000;   % Black
    1.0000, 0.0000, 1.0000;   % Magenta
    0.5000, 0.5000, 0.5000];  % Grey

sig = 0.05;

% Load PWV
NPar = length(direc);
NVes = length(vessel);
for i = 1:NPar
    for j = 1:NVes
        load(fullfile(direc{i}, ['PWV_' direc{i}(end) '_' vessel{j} '.mat']), 'pwv', 'Method', 'PrePro')
        if i == 1 && j == 1
            [NRep, NMet, NPre] = size(pwv);
            pwv_load = zeros(NPar, NVes, NRep, NMet, NPre);
        end
        pwv_load(i,j,:,:,:) = pwv;
    end
end
pwv = pwv_load; %(NPar NVes NRep NMet NPre)

% COMPARE PWV METHODS
pwv_mnVesRep = squeeze(mean(mean(pwv, 2), 3)); %(NPar NMet NPre)
pwv_mnVes_sdRep = squeeze(std(mean(pwv, 2), [], 3)); %(NPar NMet NPre)
pwv_mnVes_cvRep = pwv_mnVes_sdRep ./ pwv_mnVesRep * 100; %(NPar NMet NPre)
disp('KEY STATISTICS:')

% Mean PWV - Main method comparison
compare_pwv(pwv_mnVesRep, Method, PrePro, {'XCorr 0-80%', 'Filtered'; 'TTU', 'Filtered'; 'TTF 20-80%', 'Filtered'; 'TTF Sliding', 'Filtered'}, 'Mean PWV (m/s)', [5,14], [0.18,0.07], clr_par, 0, sig)
xticklabels({'XCorr', 'TTU', 'TTF_{20-80%}', 'TTF_{Sliding}'})
ax = gca; ax.TickLabelInterpreter = 'tex';
pbaspect([1.3,1,1]);

% Coefficient of Variation PWV - Main method comparison
compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'XCorr 0-80%', 'Filtered'; 'TTU', 'Filtered'; 'TTF 20-80%', 'Filtered'; 'TTF Sliding', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,20], [0.42,0.07], clr_par, 0, sig)
xticklabels({'XCorr', 'TTU', 'TTF_{20-80%}', 'TTF_{Sliding}'})
ax = gca; ax.TickLabelInterpreter = 'tex';
pbaspect([1.3,1,1]);

% Coefficient of Variation PWV - Interpolated vs filtered
compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'XCorr 0-80%', 'Interpolated'; 'XCorr 0-80%', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,35], [0.7,0.07], clr_par, 1, 0)
xticklabels({'Interpolated', 'Filtered'})
title('XCorr')
pbaspect([1,1.5,1]);

compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'TTU', 'Interpolated'; 'TTU', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,35], [0.7,0.07], clr_par, 1, 0)
xticklabels({'Interpolated', 'Filtered'})
title('TTU')
pbaspect([1,1.5,1]);

compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'TTF 20-80%', 'Interpolated'; 'TTF 20-80%', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,35], [0.7,0.07], clr_par, 1, 0)
xticklabels({'Interpolated', 'Filtered'})
title('TTF_{20-80%}')
pbaspect([1,1.5,1]);

compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'TTF Sliding', 'Interpolated'; 'TTF Sliding', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,35], [0.7,0.07], clr_par, 1, 0)
xticklabels({'Interpolated', 'Filtered'})
title('TTF_{Sliding}')
pbaspect([1,1.5,1]);

% Mean PWV - XCorr region comparison
compare_pwv(pwv_mnVesRep, Method, PrePro, {'XCorr 0-50%', 'Filtered'; 'XCorr 0-80%', 'Filtered'; 'XCorr 0-100%', 'Filtered'; 'XCorr 20-80%', 'Filtered'}, 'Mean PWV (m/s)', [5,14], [0.18,0.07], clr_par, 0, 0)
xticklabels({'0-50%', '0-80%', '0-100%', '20-80%'})
pbaspect([1.3,1,1]);

% Coefficient of Variation PWV - XCorr region comparison
compare_pwv(pwv_mnVes_cvRep, Method, PrePro, {'XCorr 0-50%', 'Filtered'; 'XCorr 0-80%', 'Filtered'; 'XCorr 0-100%', 'Filtered'; 'XCorr 20-80%', 'Filtered'}, 'PWV Coefficient of Variation (%)', [0,20], [0.42,0.07], clr_par, 0, 0)
xticklabels({'0-50%', '0-80%', '0-100%', '20-80%'})
pbaspect([1.3,1,1]);

% PLOT PWV AGAINST AGE
% Load Age
age_table = readtable(age_file);
age = zeros(NPar, 1); %(NPar 1)
for i = 1:NPar
    age(i) = age_table{strcmp(age_table.Participant, direc{i}(end)), 2};
end

% Compare Age and PWV
compare_age(pwv, age, Method, PrePro, {'XCorr 0-80%', 'Filtered'; 'TTU', 'Filtered'; 'TTF 20-80%', 'Filtered'; 'TTF Sliding', 'Filtered'}, [6,13], clr_par)


% CALCULATE AVERAGE PATHLENGTH & PTT
% Pathlength
pl_table = readtable(pl_file);
for i = 1:NPar
    for j = 1:NVes
        pl_ij = pl_table{strcmp(pl_table.Participant, direc{i}(end)) & strcmp(pl_table.Vessel, vessel{j}), 3:7}; %(1 NSli)
        if i == 1 && j == 1
            NSli = length(pl_ij);
            pl = zeros(NPar, NVes, NSli);
        end
        pl(i,j,:) = pl_ij;
    end
end
disp([newline 'MEAN PATHLENGTH & PTT:'])
disp(['Mean left pathlength: ' num2str(mean(pl(:,1,5)), 4) ' +- ' num2str(std(pl(:,1,5)), 4) ' mm'])
disp(['Mean right pathlength: ' num2str(mean(pl(:,2,5)), 4) ' +- ' num2str(std(pl(:,2,5)), 4) ' mm'])
disp(['Mean averaged pathlength: ' num2str(mean(mean(pl(:,:,5), 2)), 4) ' +- ' num2str(std(mean(pl(:,:,5),2)), 4) ' mm'])

% Pulse Transit Time
ptt_all = zeros(NPar, NVes, NSli, NRep, NMet, NPre); %(NPar NVes NSli NRep NMet NPre)
for i = 1:NPar
    for j = 1:NVes
        load(fullfile(direc{i}, ['PWV_' direc{i}(end) '_' vessel{j} '.mat']), 'ptt_ms'); %(NSli NRep NMet NPre)     
        ptt_all(i,j,:,:,:,:) = ptt_ms;
    end
end
ptt_mnVesRep = squeeze(mean(mean(ptt_all, 2), 4)); %(NPar NSli NMet NPre)
method_prepro = {'XCorr 0-80%', 'Filtered';...
                 'TTU', 'Filtered';...
                 'TTF 20-80%', 'Filtered';...
                 'TTF Sliding', 'Filtered'};
disp([newline 'Mean Slice 1-5 PTTs across participants (after left-right averaging):'])
for i = 1:size(method_prepro, 1)
    ptt_i = ptt_mnVesRep(:, 5, strcmp(Method, method_prepro{i,1}), strcmp(PrePro, method_prepro{i,2}));
    disp([method_prepro{i,1} ' (' method_prepro{i,2} '): ' num2str(mean(ptt_i), 4) ' +- ' num2str(std(ptt_i), 4) 'ms'])
end

end %Function end


function compare_pwv(pwv_all, method_all, prepro_all, method_prepro, yax, ybd, jtr, clr, lines, alpha)
    % pwv_all - All vessel-averaged inputs (NPar NMet NPre)
    % method_all - All methods             (1 NMet)
    % prepro_all - All preprocessing       (1 NPre)
    % method_prepro - Selected methods     {Met1, Pre1; Met2, Pre2; ...}
    % yax  - Label for y-axis              (string)
    % ybd  - Limits for y-axis             [ylb, yub]
    % jtr  - Minimum jitter distances      [ry, rx]
    % clr  - Colour of datapoints          (NPar 3)
    % lines - Connect datapoints           (bool)
    % alpha - Significance level           (scalar) (set to 0 to stop statistical tests)

    % Compile PWVs for selected methods
    [NPar, ~, ~] = size(pwv_all);
    LMet = size(method_prepro, 1);
    pwv = zeros(NPar, LMet);
    for i = 1:LMet
        pwv(:,i) = pwv_all(:, strcmp(method_all, method_prepro{i,1}), strcmp(prepro_all, method_prepro{i,2}));
    end
    
    % Plot figure
    x = (1:LMet) - 0.5;
    figure;
    set(gcf, 'WindowState', 'maximized');
    hold on
    for i = 1:LMet
        boxplot(pwv(:,i), 'Position', x(i), 'Colors', 'k', 'Whisker', inf, 'Width', 0.5)
        set(findobj(gca,'Tag','Upper Whisker'),'LineStyle','-')
        set(findobj(gca,'Tag','Lower Whisker'),'LineStyle','-')
        set(findobj(gca,'Type','Line'),'LineWidth', 2)
    end
    x_bs = zeros(NPar, LMet);
    for i = 1:LMet
        x_bs(:,i) = BeeSwarm(max(min(pwv(:,i), ybd(2)), ybd(1)), x(i), jtr(1), jtr(2));
        for j = 1:NPar
            scatter(x_bs(j,i), max(min(pwv(j,i), ybd(2)), ybd(1)), 200, clr(j,:), 'filled')
        end
    end
    if lines == 1
        for j = 1:NPar
            plot(x_bs(j,:), max(min(pwv(j,:), ybd(2)), ybd(1)), 'Color', clr(j,:), 'LineWidth', 1)
        end
    end
    xticks(x)
    LMethods = strcat(method_prepro(:,1), ' (', method_prepro(:,2), ')');
    xticklabels(LMethods)
    xlim([0, LMet])
    ylabel(yax)
    ylim(ybd)
    fontsize(28, 'points')
    ax = gca;
    ax.LineWidth = 3.0;
    box on

    if alpha > 0
        % Run Wilcoxon signed rank tests
        NTst = LMet*(LMet-1)/2;
        p_wsr = zeros(NTst, 1);
        idx = 1;
        for i = 1:LMet
            for j = (i+1):LMet
                p_wsr(idx) = signrank(pwv(:,i), pwv(:,j));
                idx = idx + 1;
            end
        end
        p_wsr_hb = Holm_Bonferroni(p_wsr);
        disp([newline yax ' - Wilcoxon Signed Rank:'])
        idx = 1;
        NSig = 0;
        for i = 1:LMet
            for j = (i+1):LMet
                if (p_wsr_hb(idx)<alpha)
                    disp([LMethods{i} ' & ' LMethods{j} ': p = ' num2str(p_wsr_hb(idx), 4) ' (Uncorrected: '  num2str(p_wsr(idx), 4) ')'])
                    NSig = NSig + 1;
                end
                idx = idx + 1;
            end
        end
        if NSig == 0
            disp('No significant differences')
        end

        % Display median, maximum, mean and standard deviation
        disp([newline yax ' - Summary statistics:'])
        for i = 1:LMet
            disp([LMethods{i} ' Median = ' num2str(median(pwv(:,i)), 4), ', Maximum = ' num2str(max(pwv(:,i)), 4), ', Mean = ' num2str(mean(pwv(:,i)), 4), ', SD = ' num2str(std(pwv(:,i)), 4)])
        end
    end

end


function compare_age(pwv_all, age, method_all, prepro_all, method_prepro, ybd, clr)
    % pwv_all - All PWV values             (NPar NVes NRep NMet NPre)
    % age - Participant age                (NPar 1)
    % method_all - All methods             (1 NMet)
    % prepro_all - All preprocessing       (1 NPre)
    % method_prepro - Selected methods     {Met1, Pre1; Met2, Pre2; ...}
    % ybd  - Limits for y-axis             [ylb, yub]
    % clr  - Colour of datapoints          (NPar 3)

    % Compile PWVs for selected methods
    pwv_all_mn = squeeze(mean(mean(pwv_all, 2), 3));    %(NPar NMet NPre)
    pwv_all_sd = squeeze(std(mean(pwv_all, 2), [], 3)); %(NPar NMet NPre)
    NPar = length(age);
    LMet = size(method_prepro, 1);
    pwv_mn = zeros(NPar, LMet); %(NPar LMet)
    pwv_sd = zeros(NPar, LMet); %(NPar LMet)
    for i = 1:LMet
        pwv_mn(:,i) = pwv_all_mn(:, strcmp(method_all, method_prepro{i,1}), strcmp(prepro_all, method_prepro{i,2}));
        pwv_sd(:,i) = pwv_all_sd(:, strcmp(method_all, method_prepro{i,1}), strcmp(prepro_all, method_prepro{i,2}));
    end
    LMethods = strcat(method_prepro(:,1), ' (', method_prepro(:,2), ')');
    [r_age, p_age] = deal(zeros(LMet, 1)); %(LMet 1)

    % Plot figures
    figure;
    set(gcf, 'WindowState', 'maximized');
    tiledlayout('TileSpacing', 'compact');
    for i = 1:LMet
        [r_age(i), p_age(i)] = corr(age, pwv_mn(:,i));
    end
    p_age_hb = Holm_Bonferroni(p_age);
    for i = 1:LMet
        nexttile
        hold on
        for j = 1:NPar
            scatter(age(j), max(min(pwv_mn(j,i), ybd(2)), ybd(1)), 100 , clr(j,:), 'Marker', 'x', 'LineWidth', 2.5)
            errorbar(age(j), max(min(pwv_mn(j,i), ybd(2)), ybd(1)), pwv_sd(j,i), 'Color', clr(j,:), 'LineStyle', 'none', 'LineWidth', 1)
        end
        if abs(r_age(i)) > 0.5
            lrm = fitlm(age, pwv_mn(:,i));
            lrm_x = (min(age):0.1:max(age)).';
            [lrm_y, lrm_ci] = predict(lrm, lrm_x);
            hold on
            plot(lrm_x, lrm_y, 'k', 'LineWidth', 3)
            plot(lrm_x, lrm_ci(:,1), 'k--', 'LineWidth', 1)
            plot(lrm_x, lrm_ci(:,2), 'k--', 'LineWidth', 1)
            leg = repmat({''}, 1, 2*NPar + 2);
            leg{end-1} = ['r = ' num2str(r_age(i), 2) ', p = ' num2str(p_age(i), 3) ', p (Holm) = ' num2str(p_age_hb(i), 3)];
            leg{end} = '95% CI';
            legend(leg, 'Location', 'northwest')
        else
            plot(nan, nan, 'LineStyle','none', 'Marker','none');
            leg = repmat({''}, 1, 2*NPar + 1);
            leg{end} = ['r = ' num2str(r_age(i), 2) ', p = ' num2str(p_age(i), 3) ', p (Holm) = ' num2str(p_age_hb(i), 3)];
            blank_leg = legend(leg,  'Location', 'northeast');
            blank_leg.ItemTokenSize = [1, 1];
        end
        xlim([min(age)-2, max(age)+2])
        xlabel('Age (years)')
        ylim(ybd)
        ylabel('Mean PWV (m/s)')
        title([LMethods{i}])
        set(gca, 'LineWidth', 1.5)
        fontsize(18, 'points')
        box on
    end
end













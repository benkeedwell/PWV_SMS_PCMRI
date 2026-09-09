%-------------------------------------------------------------------------------%
% ESTIMATE PULSE WAVE VELOCITY & COMPARE ACROSS METHODS
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Load waveforms from all repeats & estimate PWV with various methods
%-------------------------------------------------------------------------------%
% INPUTS
% TS_file - TS file names for all repeats (cell array of strings)
% pl_file - Pathlength data file name (string)
% id - Participant ID (string)
% vessel - Target vessel ('LICA' or 'RICA')
% direc - Directory for TS files (string)
% save_data - Save PWV file (bool)
% fig - Generate figures (0 - No Figures, 1 - All Figures, 
%                         2 - Technical Note Figures)
%-------------------------------------------------------------------------------%
% SAVED
% PWV - Pulse wave velocities
%-------------------------------------------------------------------------------%

function PWV3_EstimatePWV(TS_file, pl_file, id, vessel, direc, save_data, fig)

% SETUP
% Load waveforms
disp('Loading Files...')
NSet = length(TS_file); 
for i = 1:NSet
    load(fullfile(direc, TS_file{i}), 'ts_phs', 'TR', 'TE') %(NSli NCph)
    if i == 1
        [NSli, NCph] = size(ts_phs);
        ts_phs_all = zeros(NCph, NSli, NSet);
    end
    ts_phs_all(:,:,i) = ts_phs.';   %(NCph NSli NSet)
end

clr_sli = [
            0.0902 0.4667 0.7019  % Blue
            0.9608 0.4980 0.1373  % Orange
            0.2706 0.6275 0.3333  % Green
            0.9412 0.2314 0.1255  % Red
            0.6000 0.3020 0.6392  % Purple
            ]; %(NSli, 3)

clr_set = [[(NSet:-2:1), zeros(1,floor(NSet/2))]'/NSet, ... %R
           [(2:2:NSet), (NSet:-2:1)]'/NSet, ...             %G
           [zeros(1,floor(NSet/2)), (1:2:NSet)]'/NSet] ...  %B
           * 0.75; %(NSet, 3)


% ESTIMATE PTT
disp('Estimating PTT...')
Method = {'XCorr 0-50%', 'XCorr 0-80%', 'XCorr 0-100%', 'XCorr 20-80%',...
          'TTU', 'TTF 20-80%', 'TTF Sliding'}; %(1 NMet)
PrePro = {'Interpolated', 'Filtered'}; %(1 NPre)
NMet = length(Method);
NPre = length(PrePro);

ptt = zeros(NSli, NSet, NMet, NPre); %(NSli NSet NMet NPre)
for i = 1:NSet
    for j = 1:NPre

        % SLIDING PEARSON'S
        ptt(:, i, strcmp(Method, 'XCorr 0-50%'), j) =  PTTFilterCorr(ts_phs_all(:,:,i).', j-1, 0, 50, 'Pearsons', 0);
        ptt(:, i, strcmp(Method, 'XCorr 0-80%'), j) =  PTTFilterCorr(ts_phs_all(:,:,i).', j-1, 0, 80, 'Pearsons', 0);
        ptt(:, i, strcmp(Method, 'XCorr 0-100%'), j) = PTTFilterCorr(ts_phs_all(:,:,i).', j-1, 0, 100, 'Pearsons', 0);
        ptt(:, i, strcmp(Method, 'XCorr 20-80%'), j) = PTTFilterCorr(ts_phs_all(:,:,i).', j-1, 20, 80, 'Pearsons', 0);

        % MAXIMUM FIRST DERIVATIVE
        ptt(:, i, strcmp(Method, 'TTU'), j) = PATFilterUpstroke(ts_phs_all(:,:,i).', j-1, 0);

        % INTERSECTION OF TANGENTS
        ptt(:, i, strcmp(Method, 'TTF 20-80%'), j) =  PATFilterTangents(ts_phs_all(:,:,i).', j-1, 20, 80, [], 1, 0);
        ptt(:, i, strcmp(Method, 'TTF Sliding'), j) = PATFilterTangents(ts_phs_all(:,:,i).', j-1, 20, 80, 30, 1, 0);

    end
end

% Convert PTT to ms
ptt = ptt - ptt(1,:,:,:);
ptt_ms = (ptt * TR); %(NSli NSet NMet NPre)


% ESTIMATE PWV
% Load pathlengths
pl_table = readtable(pl_file);
pl = pl_table{strcmp(pl_table.Participant, id) & strcmp(pl_table.Vessel, vessel), 3:7}; %(1 NSli)

% Estimate linear fits
fit_param = zeros(NSet, NMet, NPre, 2);
for i = 1:NSet
    for j = 1:NMet
        for k = 1:NPre
            fit_i = fit(pl.', ptt_ms(:,i,j,k), 'poly1');
            fit_param(i,j,k,:) = [fit_i.p1, fit_i.p2];
        end
    end
end
pwv = (1 ./ fit_param(:,:,:,1)); %(NSet NMet NPre)


% SAVE PWV
save_name = fullfile(direc, ['PWV_' id '_' vessel '.mat']);
if save_data == 1 && ~isfile(save_name)
    save(save_name, 'pwv', 'ptt_ms', 'pl', 'Method', 'PrePro', 'TR', 'TE', '-v7.3');
elseif save_data == 1
    error('File with the requested name already exists.')
end


% FIGURES
if fig == 1

    % WAVEFORM FIGURES
    % All waveforms for the same slice
    figure(color = 'white');
    t = tiledlayout('TileSpacing', 'compact');
    title(t, ['Participant ' id ' ' vessel])
    for i = 1:NSli
        nexttile
        hold on
        for j = 1:NSet
            plot(ts_phs_all(:,i,j), 'Color', clr_set(j,:))
        end
        xlabel('Cardiac Phase')
        ylabel('Mean Phase-shift')
        title(['Slice ', num2str(i)])
        legend(cat(2, repmat('Repeat ', NSet, 1), num2str((1:NSet).')))
        box on
    end
    title(t, 'All waveforms for each slice')

    % All waveforms together
    f = figure(color = 'white');
    hold on
    for i = 1:NSet
        for j = 1:NSli
            plot(ts_phs_all(:,j,i), 'Color', clr_sli(j,:))
        end
    end
    xlabel('Cardiac Phase')
    ylabel('Mean Phase-shift')
    title(f.Children, ['All waveforms - Participant ' id ' ' vessel])
    legend(cat(2, repmat('Slice ', NSli, 1), num2str((1:NSli).')))
    box on;
   
    % PTT FIGURE
    figure(color = 'white');
    t = tiledlayout(NPre*2, ceil(NMet/2), 'TileSpacing', 'tight');
    title(t, ['Participant ' id ' ' vessel])
    for j = 1:NPre
        for i = 1:NMet
            nexttile
            for k = 1:NSet
                scatter(pl, ptt_ms(:,k,i,j), 80 , clr_set(k,:), 'Marker', 'x', 'LineWidth', 2)
                hold on
                plot(pl, fit_param(k,i,j,1)*pl + fit_param(k,i,j,2), 'Color', clr_set(k,:))
            end
            xlabel('Pathlength (mm)')
            xlim([-20, pl(end)+20])
            xline(pl, '--')
            ylabel('PTT (ms)')
            title({[Method{i} ' - ' PrePro{j}], ['PWV = ' num2str(mean(pwv(:,i,j)), 3) ' +- ' num2str(std(pwv(:,i,j)), 2) ...
                ' (' num2str(std(pwv(:,i,j)) ./ mean(pwv(:,i,j)) * 100, 2) '%)']})
            box on; pbaspect([1,1,1]);
        end
        if mod(NMet, 2) == 1
            ax = nexttile;
            axis(ax, 'off')
        end
    end

elseif fig == 2

    % Figure 3: PTT vs Pathlength
    met = strcmp(Method, 'XCorr 0-80%');
    pre = strcmp(PrePro, 'Filtered');

    figure
    set(gcf, 'WindowState', 'maximized');
    for k = 1:NSet
        scatter(pl, ptt_ms(:,k,met,pre), 250 , clr_set(k,:), 'Marker', 'x', 'LineWidth', 3)
        hold on
        plot(pl, fit_param(k,met,pre,1)*pl + fit_param(k,met,pre,2), 'Color', clr_set(k,:), 'LineWidth', 1.5)
    end
    xlabel('Pathlength (mm)')
    xlim([-20, pl(end)+20])
    xticks(0:25:pl(end)+20)
    xline(pl, '--', 'LineWidth', 1, 'Alpha', 1)
    ylabel('Pulse Transit Time (ms)')
    pbaspect([1,1,1]);
    fontsize(28, 'points')
    ax = gca;
    ax.LineWidth = 3.0;
    box on;

 end

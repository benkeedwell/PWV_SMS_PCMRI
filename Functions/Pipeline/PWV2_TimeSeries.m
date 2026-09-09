%-------------------------------------------------------------------------------%
% VESSEL SEGMENTATION AND TIME SERIES EXTRACTION FOR PWV ESTIMATION
% By Benjamin Keedwell (2026)
% Background correction by Aaron Hess
%-------------------------------------------------------------------------------%
% Segment vessel and extract mean phase difference waveform for target vessel.
%-------------------------------------------------------------------------------%
% INPUTS
% IMG_file - IMG file name (string)
% read_direc - Directory for IMG file (string)
% vessel - Target vessel ('LICA' or 'RICA')
% old_mask - Mask setting (0 - Generate new mask, 1 - Reload segmentation mask,
%                          2 - Grow new mask from previous vessel centre)
% fig - Generate figures (0 - No Figures, 1 - All Figures, 
%                         2 - All Figures but voxel-wise, 3 - Only waveforms,
%                         4 - Technical Note Figures )
% save_data - Save TS file (bool)
% save_direc - Directory to save TS file (string)
%-------------------------------------------------------------------------------%
% SAVED
% TS - Mean phase difference time series
% Vessel_Masks - Reusable vessel masks
%-------------------------------------------------------------------------------%

function PWV2_TimeSeries(IMG_file, read_direc, vessel, old_mask, fig, save_data, save_direc)

% OPTIONS/PARAMETERS
vessel_perc = 50;       % Central percentage of segmentation mask
vox_min = 10;           % Minimum number of voxels to maintain in segmentation

% SETUP
% Load image data 
disp('Loading Files...')
load(fullfile(read_direc, IMG_file), 'img', 'TR', 'TE'); %(NFrq NPhs NSli NCph)
[NFrq, NPhs, NSli, NCph] = size(img);
mag = abs(img);
phs = angle(img);

% If requested, load previous mask
if old_mask == 1 || old_mask == 2
    load(fullfile(save_direc, 'Vessel_Masks'), 'masks_saved');
end

% PHASE UNWRAPPING
disp('Unwrapping Phase...')
for i = 1:NSli
    mask_wraps_i = IdentifyWraps(squeeze(phs(:,:,i,:)), 7, 0);
    phs_lin_i = reshape(phs(:,:,i,:), NFrq*NPhs, NCph);
    phs_lin_i(mask_wraps_i(:), :) = unwrap(phs_lin_i(mask_wraps_i(:), :), [], 2);
    phs(:,:,i,:) = reshape(phs_lin_i, NFrq, NPhs, 1, NCph);
end

% PHASE BACKGROUND CORRECTION
disp('Correcting Background...')
bgfit = BackgroundCorrection(phs, pi, 2);
phs = phs - bgfit; %(NFrq NPhs NSli NCph) - (NFrq NPhs NSli)

% SEGMENT VESSEL
% Setup
disp('Segmenting Vessel...')
[masks, masks_seg] = deal(zeros(NFrq, NPhs, NSli));
[ts_phs, ts_mag] = deal(zeros(NSli, NCph));

for i = 1:NSli
    disp(['Slice = ' num2str(i) ' / ' num2str(NSli)])

    % Segment vessels and select vessel
    if old_mask == 0
        masks_seg(:,:,i) = SegmentandSelect(squeeze(phs(:,:,i,:)), squeeze(mag(:,:,i,:)), 0, [], vessel);
    elseif old_mask == 1
        masks_seg(:,:,i) = masks_saved(:,:,i);
    elseif old_mask == 2
       [MY, MX] = find(logical(masks_saved(:,:,i))); 
       masks_seg(:,:,i) = SegmentandSelect(squeeze(phs(:,:,i,:)), squeeze(mag(:,:,i,:)), 0, round([mean(MY), mean(MX)]), vessel);
    end
    
    % Select central percentage of voxels
    vessel_perc_i = min(max(vessel_perc, vox_min * 100 / length(find(masks_seg(:,:,i)))), 100);
    masks(:,:,i) = VesselCentre(masks_seg(:,:,i), vessel_perc_i, 0);

    % Average waveforms across vessel
    ts_phs(i,:) = mean_time_series(squeeze(phs(:,:,i,:)), masks(:,:,i), 0);
    ts_mag(i,:) = mean_time_series(squeeze(mag(:,:,i,:)), masks(:,:,i), 0);
end

% SAVE TIME SERIES
disp('Waveform Extraction Complete.')
SaveName = fullfile(save_direc, ['TS' IMG_file(4:end-4) '_' vessel '.mat']);
if save_data == 1 && ~isfile(SaveName)
    save(SaveName, 'ts_phs', 'masks', 'TR', 'TE', '-v7.3');
elseif save_data == 1
    error('File with the requested name already exists.')
end

% SAVE CURRENT MASKS
if (old_mask == 0) && (save_data == 1)
    masks_saved = masks_seg;
    save(fullfile(save_direc, 'Vessel_Masks'), 'masks_saved');
end

% FIGURES
if fig > 0 

    clr = [
           0.0902 0.4667 0.7019  % Blue
           0.9608 0.4980 0.1373  % Orange
           0.2706 0.6275 0.3333  % Green
           0.9412 0.2314 0.1255  % Red
           0.6000 0.3020 0.6392  % Purple
          ];

    % Voxel-wise time series
    if fig == 1
        for i = 1:NSli
            [row, col] = find(masks_seg(:,:,i));
            rec_i = [min(col)-2, min(row)-2, max(col)+2, max(row)+2];
            if (rec_i(4)-rec_i(2))*(rec_i(3)-rec_i(1)) < 400
                phs_i = squeeze(phs(:,:,i,:));
                VoxelWaveforms(phs_i, rec_i, masks_seg(:,:,i), masks(:,:,i),...
                    [0, max(phs_i(rec_i(2):rec_i(4), rec_i(1):rec_i(3), :), [], 'all') + 0.2], clr(i,:), ['Slice ' num2str(i)])
            end
        end
    end

    if any(fig == [1,2])
        % ROI on cardiac-cycle-averaged images
        figure(Color = 'white');
        t_mag = tiledlayout('TileSpacing', 'tight');
        for i = 1:NSli
            nexttile
            img_mask_scale(mean(squeeze(mag(:,:,i,:)), 3), 'prc', 2, 98, masks(:,:,i))
        end
        title(t_mag, 'ROI on magnitude images (averaged temporally)')

        figure(Color = 'white');
        t_phs = tiledlayout('TileSpacing', 'tight');
        for i = 1:NSli
            nexttile
            img_mask_scale(mean(squeeze(phs(:,:,i,:)), 3), 'raw', -pi/2, pi/2, masks(:,:,i))
        end
        title(t_phs, 'ROI on phase images (averaged temporally)')
    end

    if any(fig == [1,2,3])

        figure(color = 'white');
        tl = tiledlayout('TileSpacing', 'tight');
        title(tl, ['TS' IMG_file(4:end-4) '_' vessel], 'Interpreter', 'none')

        %Plot all phase difference time series
        nexttile
        legpd = cell(NSli, 1);
        hold on
        for i = 1:NSli
            plot(ts_phs(i,:), 'Color', clr(i,:), 'LineStyle', '-', 'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 8);
            legpd{i} = append('Slice-', string(i));
        end
        title('Phase Difference Time Series')
        legend(legpd)
        xlabel('Cardiac Phase');ylabel('Mean phase difference');
        box on

        %Plot normalised phase differece time series
        nexttile
        ts_phs_norm = zeros(size(ts_phs));
        legpdn = cell(NSli, 1);
        hold on
        for i = 1:NSli
            ts_phs_norm(i,:) = rescale(ts_phs(i,:));
            plot(ts_phs_norm(i,:), 'Color', clr(i,:), 'LineStyle', '-', 'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 8);
            legpdn{i} = append('Slice-', string(i));
        end
        title('Normalised Phase Difference Time Series')
        legend(legpdn)
        xlabel('Cardiac Phase');ylabel('Mean phase difference');
        box on

    end

    if any(fig == [1,2])
        %Plot all magnitude time series
        legmag = cell(NSli, 1);
        f_ts_mag = figure(color = 'white');
        hold on
        for i = 1:NSli
            plot(ts_mag(i,:), 'Color', clr(i,:), 'LineStyle', '-', 'LineWidth', 1.2, 'Marker', '.', 'MarkerSize', 8);
            legmag{i} = append('Slice-', string(i));
        end
        title(f_ts_mag.Children, 'Magnitude Time Series')
        legend(legmag)
        xlabel('Cardiac Phase');ylabel('Mean magnitude');
        box on
    end

    if fig == 4
        % Figure 1b: PC-MRI
        ExCph = 20;
        for i = 1:NSli
            f_pcmri = figure;
            set(gcf, 'WindowState', 'maximized');
            imshow(phs(:,:,i,ExCph), [-pi, pi])
            title(f_pcmri.Children, num2str(i))
        end
   
        % Figure 1c: Segmentation
        f_seg = figure;
        set(gcf, 'WindowState', 'maximized');
        ExCph = 20;
        ExSli = 1;
        imshow(phs(:,:,ExSli,ExCph), [-pi, pi])
        title(f_seg.Children, num2str(ExSli))
        clr_cen = cat(3, 0*ones(NFrq, NPhs), 255*ones(NFrq, NPhs), 0*ones(NFrq, NPhs)); %RGB
        mask_clr_cen = clr_cen .* repmat(masks(:,:,ExSli), [1,1,3]);
        hold on
        pause(1)
        m_cen = imshow(mask_clr_cen);
        m_cen.AlphaData = masks(:,:,ExSli)*0.3;
        clr_all = cat(3, 0*ones(NFrq, NPhs), 0*ones(NFrq, NPhs), 255*ones(NFrq, NPhs)); %RGB
        mask_clr_all = clr_all .* repmat(masks_seg(:,:,ExSli)-masks(:,:,ExSli), [1,1,3]);
        pause(1)
        m_all = imshow(mask_clr_all);
        m_all.AlphaData = masks_seg(:,:,ExSli)*0.3;
        axis([155, 177, 85, 107])
    end

end

end

% FUNCTIONS
% Function to show image scaled by percentage (or raw values) with mask overlayed 
function img_mask_scale(img, type, lower, upper, mask)
    R=0; G=255; B=0;
    clr = cat(3, R*ones(size(img)), G*ones(size(img)), B*ones(size(img)));
    mask_clr = clr .* repmat(mask, [1,1,3]);
    if strcmp(type, 'prc')
        c_scale = [prctile(img, lower, 'all'), prctile(img, upper, 'all')];
    elseif strcmp(type, 'raw')
        c_scale = [lower, upper];
    end
    imshow(img, c_scale)
    hold on
    h = imshow(mask_clr);
    set(h, 'AlphaData', 0.3*mask);
end

% Calculate mean time series from masked time series
function ts_mean = mean_time_series(img, mask, figs)    
    [~, ~, NCph] = size(img);
    if figs == 1
        figure, img_mask_scale(mask, 'raw', 0, 1, mask)
    end
    ts_mean = zeros(NCph, 1);
    for j = 1:NCph
        img_j=img(:,:,j);
        ts_mean(j) = mean(img_j(mask>0));
    end
end




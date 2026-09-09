function [msk] = SegmentandSelect(phs, mag, fig, seed, vessel)

%-------------------------------------------------------------------------------%
% SEGMENT VESSELS AND SELECT REGION-OF-INTEREST
% By Benjamin Keedwell (2026)
% Vessel selection UI based on a similar function by Yang Ji
%-------------------------------------------------------------------------------%
% Segment vessels from maximum phase-shift and smoothness before selecting ROI.
%-------------------------------------------------------------------------------%
% INPUTS
% phs - Phase-shift image (NFrq NPhs NCph)
% mag - Magnitude image (NFrq NPhs NCph)
% fig - Show figures (boolean)
% seed - Optional coordinates to grow mask from ([] OR [y, x])
% vessel - Target vessel ('LICA' or 'RICA')
%-------------------------------------------------------------------------------%
% OUTPUTS
% msk - Final vessel mask (NFrq NPhs)
%-------------------------------------------------------------------------------%

% Threshold by maximum phase-shift
msk_max = MaxPhaseMask(phs, 0, pi/3, fig);

% Threshold by smoothness
msk_smo = SmoothnessMask(phs, 0.85, fig);

% Combine masks
msk_seg = msk_max & msk_smo;


% OPTION 1: SELECT VESSEL
if isempty(seed) || ~msk_seg(seed(1), seed(2))
    % Initialise empty final mask
    msk = false(size(msk_seg));

    % Mean magnitude for interactive figure
    mean_mag = mean(mag, 3); %(NFrq NPhs)

    % Create figure with user interface button to break while loop
    f = figure('Color', 'white', 'Units', 'normalized', 'Position', [0,0,1,1]);
    H = uicontrol(f, 'Style', 'PushButton', ...
        'String', 'Selection Complete', ...
        'Callback', 'delete(gcbf)',...
        'Units', 'normalized',...
        'Position', [0.5-0.05, 0.04, 0.16 0.05]);
    ImgScale(mean_mag, 2, 98)
    if strcmp(vessel, 'LICA')
        title(f.Children(2), 'Draw ROI around the left CCA/ICA (right of image)')
    elseif strcmp(vessel, 'RICA')
        title(f.Children(2), 'Draw ROI around the right CCA/ICA (left of image)')
    end
    i=1;

    while true
        % Show blue initial mask over mean magnitude image
        hold on
        ImgScale(mean_mag, 2, 98)
        msk_blue = cat(3, 0*ones(size(msk_seg)), 0*ones(size(msk_seg)), 255*ones(size(msk_seg))) .* repmat(msk_seg, [1,1,3]);
        h = imshow(msk_blue);
        set(h, 'AlphaData', 0.5 .* msk_seg);

        % Show green mask of selected voxels
        msk_green = cat(3, 0*ones(size(msk)), 255*ones(size(msk)), 0*ones(size(msk))) .* repmat(msk, [1,1,3]);
        h = imshow(msk_green);
        set(h, 'AlphaData', 0.5 .* msk);

        % Create interactive polygon tool. 'roipoly' returns a mask of the selected voxels
        msk_roi = roipoly();

        % End loop when the 'break' button is pressed
        if ~ishandle(H)
            break
        end

        % If a ROI has been selected, add ROI to final mask
        if any(msk_roi(:) > 0.5)
            mask_roi_seg = msk_roi .* msk_seg;
            msk_seg = msk_seg - msk_roi;
            msk = msk + mask_roi_seg;
            disp(['ROI ' num2str(i) ' Drawn']);
            i = i + 1;
        end
    end
    fprintf('Draw ROI ended\n');

% OPTION 2: GROW MASK FROM SEED
else
    msk_lab = labelmatrix(bwconncomp(msk_seg, 4));
    msk = (msk_lab == msk_lab(seed(1), seed(2)));
end



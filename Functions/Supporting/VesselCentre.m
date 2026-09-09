function msk_out = VesselCentre(msk_in, perc, fig)
%-------------------------------------------------------------------------------%
% SELECT CENTRAL PERCENTAGE OF INPUT MASK
% By Benjamin Keedwell (2026)
%-------------------------------------------------------------------------------%
% Select a central percentage of voxels from input mask (single vessel
% mask). To account for the elliptical nature of vessel cross-sections, 
% voxels are ordered using a weighted coordinate system based on a rotated
% ellipse.
%-------------------------------------------------------------------------------%
% INPUTS
% msk_in - Mask consisting of a single vessel cross section (NY NX)
% perc - Central percentage of voxels
% fig - Show figures (boolean)
%-------------------------------------------------------------------------------%
% OUTPUTS
% msk_out - Mask consisting of central vessel cross section (NY NX)
%-------------------------------------------------------------------------------%

% Options
NRot = 360;                   % Number of ellipse rotations to trial

% Check inputs
if ~ismatrix(msk_in)
    error('Input mask must have 2 dimensions')
elseif perc < 0 || perc > 100
    error('Input percentage must be between 0 and 100')
elseif fig ~= 1 && fig~= 0
    error('Figure parameter must be 0 or 1')
end

% Voxel coordindates and linear indices
msk_in = logical(msk_in);
[msk_in_y, msk_in_x] = find(msk_in);       % (NIn 1) (NIn 1) y- and x-coordinates of input vessel mask
NIn = length(find(msk_in));                % Number of input mask voxels
NOut = ceil(NIn * perc / 100);             % Number of output mask voxels

% Calculate centroid
CY = mean(msk_in_y);                       % Centroid y-coordinate
CX = mean(msk_in_x);                       % Centroid x-coordinate
msk_in_y_rel = msk_in_y - CY;              % (NIn 1) Relative y-coordinates of input vessel mask
msk_in_x_rel = msk_in_x - CX;              % (NIn 1) Relative x-coordinates of input vessel mask

% Determine ellipse rotation
theta_range = 0 : (pi/2) / NRot : pi/2 - ((pi/2) / NRot);                                       % (1 NRot) Number of ellipse rotations to trial
msk_in_y_allrot = msk_in_x_rel .* sin(theta_range) + msk_in_y_rel .* cos(theta_range);          % (NIn NRot) Rel. x-coordinates of input voxels in rotated coordinate systems
msk_in_x_allrot = msk_in_x_rel .* cos(theta_range) - msk_in_y_rel .* sin(theta_range);          % (NIn NRot) Rel. y-coordinates of input voxels in rotated coordinate systems
ecc_allrot = std(msk_in_y_allrot, 0, 1) ./ std(msk_in_x_allrot, 0, 1);                          % (1 NRot) Measure of eccentricity for each rotation
for i = 1:NRot                                                                                    % Ensure that the minimum eccentricity is 1 (i.e, a circle)
    if ecc_allrot(i) < 1
        ecc_allrot(i) = 1 / ecc_allrot(i);
    end
end
theta_index = find(ecc_allrot == max(ecc_allrot));
theta_index = theta_index(1);                          % For case when vessel is perfectly symmetrical and all theta are selected
theta = theta_range(theta_index);
msk_in_y_rot = msk_in_y_allrot(:, theta_index);        % (NIn) Rel. y-coordinates of input voxels in the rotated coordinate system of maximum eccentricity
msk_in_x_rot = msk_in_x_allrot(:, theta_index);        % (NIn) Rel. x-coordinates of input voxels in the rotated coordinate system of maximum eccentricity
el_b = 2 * std(msk_in_y_rel);                          %  Ellipse Equation (defined relative to centroid here):
el_a = 2 * std(msk_in_x_rel);                          %  (x - h)^2 / a^2 + (y - k)^2 / b^2 = 1
if el_b == 0                                           %  In case of 1D vessel
    el_b = 0.0001;
end
if el_a == 0
    el_a = 0.0001;
end
el_k = 0;
el_h = 0;
el_b_rot = 2 * std(msk_in_y_rot);
el_a_rot = 2 * std(msk_in_x_rot);
if el_b_rot == 0                                       %  In case of 1D vessel
    el_b_rot = 0.0001;
end
if el_a_rot == 0
    el_a_rot = 0.0001;
end
el_k_rot = 0;
el_h_rot = 0;

% Determine distance to centroid of voxels weighted by rotated ellipse axes
msk_in_dist = sqrt((msk_in_y_rot/el_b_rot).^2 + (msk_in_x_rot/el_a_rot).^2);       % (NIn 1) Weighted distance to centroid of vessel mask voxels

% Select voxels by weighted distance from centroid
[~, order] = sort(msk_in_dist);                              % (NIn 1) Reordering of voxels into ascending distance order
msk_out_y = msk_in_y(order(1:NOut));                         % (NOut 1) y-coordinates of output vessel mask
msk_out_x = msk_in_x(order(1:NOut));                         % (NOut 1) x-coordinates of output vessel mask

% Create output mask
msk_out = zeros(size(msk_in));                               % (NY NX) Create empty output mask
for i = 1:NOut
    msk_out(msk_out_y(i), msk_out_x(i)) = 1;                 % Fill-in output mask
end

% Plot figures
if fig==1
    %FIGURE 1
    f_masks = figure(Color = 'white');
    imshow(msk_in)
    hold on
    R=0; G=255; B=0;
    clr = cat(3, R*ones(size(msk_out)), G*ones(size(msk_out)), B*ones(size(msk_out)));
    mask_out_clr = clr .* repmat(msk_out, [1,1,3]);
    i_mask_out = imshow(mask_out_clr);
    set(i_mask_out, 'AlphaData', 0.5 * msk_out);
    scatter(CX, CY, '*', 'r')
    title(f_masks.Children(1), 'Input mask with overlayed output mask (green)')
    
    el_y = el_k + el_b*sin(0:2*pi/100:2*pi);
    el_x = el_h + el_a*cos(0:2*pi/100:2*pi);

    el_y_rot_tmp = el_k_rot + el_b_rot*sin(0:2*pi/100:2*pi);
    el_x_rot_tmp = el_h_rot + el_a_rot*cos(0:2*pi/100:2*pi);

    el_y_rot = el_x_rot_tmp * sin(-theta) + el_y_rot_tmp * cos(-theta);
    el_x_rot = el_x_rot_tmp * cos(-theta) - el_y_rot_tmp * sin(-theta);

    plot(el_x + CX, el_y + CY, 'g')
    plot(el_x_rot + CX, el_y_rot + CY, 'm')

    legend({'Centroid', 'Ellipse', 'Rotated Ellipse'})

    pause(2)

    %FIGURE 2
    f_order = figure(Color = 'white');
    mask_in_y_ord = msk_in_y(order);
    mask_in_x_ord = msk_in_x(order);
    bw_range = 0: 1/(NIn-1): 1;
    mask_order = zeros(size(msk_in));
    for i = 1:NIn
        mask_order(mask_in_y_ord(i), mask_in_x_ord(i)) = bw_range(i);
    end
    imshow(mask_order)
    hold on
    scatter(CX, CY, '*', 'r')
    title(f_order.Children(1), 'Voxel Order')
    el_y = el_k + el_b*sin(0:2*pi/100:2*pi);
    el_x = el_h + el_a*cos(0:2*pi/100:2*pi);

    el_y_rot_tmp = el_k_rot + el_b_rot*sin(0:2*pi/100:2*pi);
    el_x_rot_tmp = el_h_rot + el_a_rot*cos(0:2*pi/100:2*pi);

    el_y_rot = el_x_rot_tmp * sin(-theta) + el_y_rot_tmp * cos(-theta);
    el_x_rot = el_x_rot_tmp * cos(-theta) - el_y_rot_tmp * sin(-theta);

    plot(el_x + CX, el_y + CY, 'g')
    plot(el_x_rot + CX, el_y_rot + CY, 'm')
    legend({'Centroid', 'Ellipse', 'Rotated Ellipse'})
 end
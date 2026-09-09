%-------------------------------------------------------------------------------%
% GRAPPA RECONSTRUCTION FOR SMS PC-MRI DATA
% By Benjamin Keedwell (2026)
% Split slice-GRAPPA functions by Mark Chiew
% MapVBVD by Philipp Ehses with modifications by Yang Ji
%-------------------------------------------------------------------------------%
% Split slice-GRAPPA reconstruction of SMS data for no in-plane parallel 
% acceleration. Intended for reconstruction of PC-MRI data with a single 
% reference cardiac phase, used for the estimation of PWV. Assumes 75%
% readout partial Fourier.
%-------------------------------------------------------------------------------%
% INPUTS
% meas_direc - Path to meas.dat files (string)
% dataset - Participant ID (string)
% dat_SMS - SMS PC-MRI meas.dat (string)
% dat_calib - ACS meas.dat (string)
% rep - SMS PC-MRI repetition (scalar)
% fig - Generate figures (bool)
% save_data - Save IMG file (bool)
% save_label - Label for IMG file (bool)
% save_direc - Directory to save IMG file (string)
%-------------------------------------------------------------------------------%
% SAVED
% IMG - Complex slice-separated phase difference image
%-------------------------------------------------------------------------------%

function PWV1_SMSRecon(meas_direc, dataset, dat_SMS, dat_calib, rep, fig, save_data, save_label, save_direc)

% OPTIONS/PARAMETERS
kernel_size = [5,5];  % GRAPPA Kernel size
LCph = 35;            % Number of cardiac phases to reconstruct excluding reference (set to 0 to reconstruct all)
ExCph = 27;           % Example cardiac phase for figures, not including reference phase

% SETUP
% Read SMS Data
disp('Reading SMS data...')
[ksp_SMS, twix_SMS] = ReadTWIX(fullfile(meas_direc, dataset, dat_SMS), 'img'); %ksp_SMS: (NFrq NCha NPhs NCph+1 NSet)
ksp_SMS = ksp_SMS(:,:,:,:,rep); %(NFrq NCha NPhs NCph+1)
[NFrq, ~, NPhs, NCph] = size(ksp_SMS);
NCph = NCph - 1;
NSli = twix_SMS.hdr.MeasYaps.sWipMemBlock.alFree{12};
PhaseShiftFactor = twix_SMS.hdr.MeasYaps.sWipMemBlock.alFree{14};
TR = twix_SMS.hdr.Config.TR / 1000;
TE = twix_SMS.hdr.Meas.alTE(1) / 1000;
ksp_SMS = permute(ksp_SMS,[2,1,3,5,4]); %(NCha NFrq NPhs 1 NCph+1)

% Read calibration data
disp('Reading calibration data...')
[ksp_calib, ~] = ReadTWIX(fullfile(meas_direc, dataset, dat_calib), 'acs'); %ksp_calib: (NFrq NCha NPhs_calib NSli)
ksp_calib = permute(ksp_calib, [2,1,3,4]); %(NCha NFrq NPhs_calib NSli)

% Shift calibration data according to phase shift factor
NPhs_calib = size(ksp_calib, 3);
ksp_calib = padarray(ksp_calib, [0, 0, (NPhs - NPhs_calib)/2], 0); %(NCha NFrq NPhs NSli)
if PhaseShiftFactor == 5
    lookup = [1,4,2,5,3];
elseif PhaseShiftFactor == 4
    lookup = [1,3,2,4,5];
else
    lookup = [1,2,3,4,5];
end
for i = 1:NSli
    i_lu = lookup(i);
    ksp_calib(:,:,:,i) = ksp_calib(:,:,:,i) .* reshape(exp(-1i*2*pi*(i_lu-1)*(-NPhs/2:NPhs/2-1)/PhaseShiftFactor), [1,1,NPhs,1]);
end
ksp_calib = ksp_calib(:,:,((NPhs - NPhs_calib)/2 + 1):(NPhs + NPhs_calib)/2,:);

% Limit number of cardiac phases if requested
if LCph > 0
    ksp_SMS = ksp_SMS(:,:,:,:,1:LCph+1); %(NCha NFrq NPhs 1 LCph+1)
else
    LCph = NCph;
end

% APPLY GRAPPA RECON, FT & COMBINE CHANNELS
img = zeros(NFrq, NPhs, NSli, LCph);

% Calibrate kernels
disp('Calibrating weights...')
w   =   weights_spsg(ksp_calib(:,(NFrq/4 + 1):end,:,:), kernel_size); 
disp('GRAPPA Reconstructing...')
for i = 1:LCph+1
    disp([num2str(i) ' / ' num2str(LCph+1)])
    ksp_SMS_i = ksp_SMS(:,(NFrq/4 + 1):end,:,:,i);

    % Apply GRAPPA Reconstruction
    ksp_sg_i = apply_weights(ksp_SMS_i, w);

    % Fermi filter to smooth asymmetric echo
    ksp_sg_i = padarray(ksp_sg_i, [0,NFrq/4,0,0], 0, 'pre'); %(NCha NFrq NPhs NSli)
    W_frq = NFrq/200;
    fermi_frq = flip((1.0 ./ (1.0 + exp((((1:NFrq))-(NFrq - NFrq/4 - 3*W_frq)) / W_frq))));
    ksp_sg_i = ksp_sg_i .* fermi_frq;

    % Fourier Transform
    img_sg_i = IFFTdim(ksp_sg_i, [2,3]);
    clear('ksp_sg_i')

    % Combine across channels and compute phase difference
    if i == 1
        img_sg_ref = img_sg_i; %(NCha NFrq NPhs NSli)
    else
        img(:,:,:,i-1) = squeeze(sum(img_sg_i .* conj(img_sg_ref), 1)); %(NFrq NPhs NSli LCph)
        img(:,:,:,i-1) = sqrt(abs(img(:,:,:,i-1))) .* exp(1i * angle(img(:,:,:,i-1)));
    end
    clear('img_sg_i')
end

% PROCESS & SAVE IMAGES
% Remove FOV shift
ksp = FFTdim(img, [1,2]);
for i = 1:NSli
    i_lu = lookup(i);
    ksp(:,:,i,:) = ksp(:,:,i,:) .* reshape(exp(1i*2*pi*(i_lu-1)*(-NPhs/2:NPhs/2-1)/PhaseShiftFactor), [1,NPhs,1,1]);
end
img = IFFTdim(ksp, [1,2]);

% Rotate images by 180 degrees and flip slice order
img = rot90(img, 2);
img = flip(img, 3); %(NFrq NPhs NSli LCph) 

% Save reconstructed images
disp('Reconstruction Complete.')
SaveName = fullfile(save_direc, dataset, ['IMG_' dataset(end) '_' save_label '.mat']);
if save_data == 1 && ~isfile(SaveName)
    save(SaveName, 'img', 'TR', 'TE', '-v7.3');
elseif save_data == 1
    error('File with the requested name already exists.')
end

% FIGURES
if fig==1

    % Example channel-combined magnitude images
    figure(color = 'white');
    t_img_mag = tiledlayout('TileSpacing', 'tight');
    for i = 1:NSli
        nexttile
        ImgScale(abs(img(:, :, i, ExCph)), 2, 98)
        title(['Slice ' num2str(i)])
    end
    title(t_img_mag, ['Example channel-combined reconstructed magnitude image for c. phase = ' num2str(ExCph)])

    % Example channel-combined phase difference images
    figure(color = 'white');
    t_img_phs = tiledlayout('TileSpacing', 'tight');
    for i = 1:NSli
        nexttile
        imshow(angle(img(:, :, i, ExCph)), [-pi, pi])
        title(['Slice ' num2str(i)])
    end
    title(t_img_phs, ['Example channel-combined reconstructed phase diff. image for c. phase = ' num2str(ExCph)])

end








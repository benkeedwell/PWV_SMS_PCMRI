%-------------------------------------------------------------------------------%
% Reconstruction and analysis for "Pulse wave velocity measurements in the 
% internal carotid arteries using simultaneous multi-slice phase contrast MRI
% as an assessment of local arterial stiffness"
%
% By Benjamin Keedwell, Aaron T. Hess, Thomas W. Okell, Iulius Dragonu, 
% Yang Ji, and Peter Jezzard
%
% Script by Benjamin Keedwell (2026)
% 
% SECTIONS
% PWV1 - SMS Reconstruction
% PWV2 - Extract Time Series
% PWV3 - Estimate PWV
% PWV4 - Compare Participants
% Generate Figures
%
% Data are provided to run sections PWV3 and PWV4. Due to concerns regarding 
% participant anonymity, data are not provided for sections PWV1 and PWV2.
% However, all the necessary functions are provided to run these sections.
%
% Similarly, only some figures from the 'Generate Figures' section can be
% generated with the provided data. Figures that cannot be generated are
% commented out but all necessary functions are provided.
%-------------------------------------------------------------------------------%


%% PWV1 - SMS Reconstruction 
% Perform SMS reconstruction using split slice-GRAPPA
% Reconstructed images are saved as .mat files (e.g. 'IMG_A_1')
clc; clear; close all; restoredefaultpath
check_set_path

% Participant ID  
DATASET = {'20260306_A', '20260317_B', '20260320_C', ...
           '20260327_D', '20260401_E', '20260409_F', ...
           '20260415_G', '20260416_1_H', '20260416_2_I', ...
           '20260429_J'};                     

% Corresponding SMS PC-MRI data
PC = {'meas_MID00225_FID85843_PC_MB5_256L_x5.dat', 'meas_MID00027_FID86513_PC_MB5_256L_x5.dat', 'meas_MID00029_FID86939_PC_MB5_256L_x5.dat', ...
        'meas_MID00029_FID87500_PC_MB5_256L_x5.dat', 'meas_MID00027_FID87802_PC_MB5_256L_x5.dat', 'meas_MID00151_FID88063_PC_MB5_256L_x5.dat', ...
        'meas_MID00026_FID88245_PC_MB5_256L_x5.dat', 'meas_MID00027_FID88326_PC_MB5_256L_x5.dat', 'meas_MID00047_FID88346_PC_MB5_256L_x5.dat', ...
        'meas_MID00028_FID89030_PC_MB5_256L_x5.dat'};

% Corresponding calibration data
ACS = {'meas_MID00222_FID85840_ACS_MB5.dat', 'meas_MID00023_FID86509_ACS_MB5.dat', 'meas_MID00023_FID86933_ACS_MB5.dat', ...
       'meas_MID00024_FID87495_ACS_MB5.dat', 'meas_MID00023_FID87798_ACS_MB5.dat', 'meas_MID00149_FID88061_ACS_MB5.dat', ...
       'meas_MID00023_FID88242_ACS_MB5.dat', 'meas_MID00023_FID88322_ACS_MB5.dat', 'meas_MID00043_FID88342_ACS_MB5.dat', ...
       'meas_MID00025_FID89027_ACS_MB5.dat'};

% Set repeats to reconstruct
N = 1:5;    

% Save reconstructed data (e.g. 'IMG_A_1.mat')
SAVE_DATA = 0;     

% Path to meas.dat files
MEAS_DIREC = '/Users/benjaminkeedwell/Data_Storage/Physiological_Variation_Study';

% Path to save IMG files
SAVE_DIREC = '/Users/benjaminkeedwell/IMG/PWV_SMS_PCMRI';                          

% Run SMS reconstruction
disp('SMS RECONSTRUCTION')

for i = 1:length(DATASET)
    
    disp(['PARTICIPANT ' DATASET{i}(end)])

    if ~exist(fullfile(SAVE_DIREC, DATASET{i}), 'dir')
        if exist(SAVE_DIREC, 'dir')
            mkdir(fullfile(SAVE_DIREC, DATASET{i}))
        else
            error('Requested save directory does not exist.')
        end
    end
    
    for j = 1:length(N)
        disp(['REPEAT ' num2str(N(j))])
        FIG = (i == length(DATASET)) && (j == length(N));
        PWV1_SMSRecon(MEAS_DIREC, DATASET{i}, PC{i}, ACS{i}, N(j), FIG, SAVE_DATA, num2str(N(j)), SAVE_DIREC)
    end

end
disp('FINISHED')


%% PWV2 - Extract Time Series
% Select vessel ROI and extract time series.
% Time series are saved as .mat files (e.g. 'TS_A_1_LICA')
clc; clear; close all; restoredefaultpath
check_set_path

% Participant ID  
DATASET = {'20260306_A', '20260317_B', '20260320_C', ...
           '20260327_D', '20260401_E', '20260409_F', ...
           '20260415_G', '20260416_1_H', '20260416_2_I', ...
           '20260429_J'}; 

% Target Vessels
VESSEL = {'LICA', 'RICA'};

% Set repeats to extract time series
N = 1:5;

% Save time series (e.g. 'TS_A_1_LICA.mat')
SAVE_DATA = 0;                   

% Path to IMG files
IMG_DIREC = '/Users/benjaminkeedwell/IMG/PWV_SMS_PCMRI';

% Mask setting (0 - Select Individually, 1 - Grow from seed for repeats)
MASK_SETTING = 0;

% Run time series extraction
disp('TIME SERIES EXTRACTION')

for i = 1:length(DATASET)
    
    disp(['PARTICIPANT ' DATASET{i}(end)])
    READ_DIREC = fullfile(IMG_DIREC, DATASET{i});
    SAVE_DIREC = fullfile(pwd, DATASET{i});

    for j = 1:length(VESSEL)

        disp(['VESSEL ' VESSEL{j}])
    
        for k = 1:length(N)
            disp(['REPEAT ' num2str(N(k))])
            IMG = ['IMG_' DATASET{i}(end) '_' num2str(N(k)) '.mat'];
            if MASK_SETTING == 0
                MASK = 0;
            elseif MASK_SETTING == 1
                MASK = 2*(k~=1);
            end
            FIG = 3*(k==length(N));
            PWV2_TimeSeries(IMG, READ_DIREC, VESSEL{j}, MASK, FIG, SAVE_DATA, SAVE_DIREC)
        end

    end

end
disp('FINISHED')


%% PWV3 - Estimate PWV
% Load time series & estimate PWV with various methods
% PWV estimates are saved as .mat files (e.g. 'PWV_A_1_LICA')
clc; clear; close all; restoredefaultpath
check_set_path

% Participant ID  
DATASET = {'20260306_A', '20260317_B', '20260320_C', ...
           '20260327_D', '20260401_E', '20260409_F', ...
           '20260415_G', '20260416_1_H', '20260416_2_I', ...
           '20260429_J'}; 

% Target Vessels
VESSEL = {'LICA', 'RICA'};

% Set repeats to extract time series
N = 1:5;

% Save time series (e.g. 'PWV_A_1_LICA.mat')
SAVE_DATA = 0;

% Pathlength data
PL = 'pl.xlsx';                                        

% Run PWV estimation
disp('PWV ESTIMATION')

for i = 1:length(DATASET)

    disp(['PARTICIPANT ' DATASET{i}(end)])
    DIREC = fullfile(pwd, DATASET{i});

    for j = 1:length(VESSEL)

        disp(['VESSEL ' VESSEL{j}])
        ID = DATASET{i}(end);
        TS = arrayfun(@(a) sprintf(['TS_' ID '_' '%d' '_' VESSEL{j}], a), N, 'UniformOutput', false);
        FIG = (i == 1) && (j == 1);
        PWV3_EstimatePWV(TS, PL, ID, VESSEL{j}, DIREC, SAVE_DATA, FIG)

    end

end
disp('FINISHED')


%% PWV4 - Compare Participants
%Load PWV across participants and compare across methods
clc; clear; close all; restoredefaultpath
check_set_path

% Participant ID 
DATASET = {'20260306_A', '20260317_B', '20260320_C', ...
           '20260327_D', '20260401_E', '20260409_F', ...
           '20260415_G', '20260416_1_H', '20260416_2_I', ...
           '20260429_J'};

% Outlying Participant Excluded
% DATASET = {'20260306_A', '20260317_B', '20260320_C', ...
%            '20260327_D', '20260401_E', '20260409_F', ...
%                          '20260416_1_H', '20260416_2_I', ...
%            '20260429_J'};

% Participant
AGE = 'age.xlsx';

% Pathlength data
PL = 'pl.xlsx';

% Run comparison
PWV4_CompareParticipants(DATASET, AGE, PL)


%% Generate Figures
clc; clear; close all; restoredefaultpath
check_set_path

% Figure 1a - TOF MIP
% GenerateFigures('TOF')

% Figure 1b-c - PC-MRI & Segmentation
% READ_DIREC = '/Users/benjaminkeedwell/IMG/PWV_SMS_PCMRI/20260429_J';
% PWV2_TimeSeries('IMG_J_1.mat', READ_DIREC, 'LICA', 0, 4, 0, '')

% Figure 1d - Waveforms
GenerateFigures('Waveforms')

% Figure 2 - PTT Methods
GenerateFigures('PTT')

% Figure 3 - PTT-Pathlength Plot
TS = arrayfun(@(a) sprintf(['TS_A_' '%d' '_LICA'], a), 1:5, 'UniformOutput', false);
DIREC = fullfile(pwd, '20260306_A');
PWV3_EstimatePWV(TS, 'pl.xlsx', 'A', 'LICA', DIREC, 0, 2)

% Figures 4-5 - Mean PWV, CV PWV & Age
%  - Figures generated by PWV4 above.

% Figures S1 & S2 - Filtering & XCorr Regions
%  - Figures generated by PWV4 above.

% Figure S3 - TTF Methods Examples
GenerateFigures('TTF')


%% FUNCTIONS
% Ensures data are organised correctly when saved.
% If the directory name has been changed, please change here too.
function check_set_path
    run_direc = strsplit(pwd, filesep);
    if strcmp(run_direc{end}, 'PWV_SMS_PCMRI')
        addpath(genpath(pwd))
    else
        error('Current directory name is not as expected')
    end
end
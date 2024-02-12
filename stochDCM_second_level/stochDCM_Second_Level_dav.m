%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/stochDCM/';

%% Experiment information 
subjects = {'001', '002', '003', '004', ...
    '005', '009', '010', '011', '012', ...
    '013', '014', '017', '018', '019', ...
    '020', '022', '023', '024', '025', ...
    '026', '029', '030', '031', '032', ...
    '033', '035', '036', '037', '038', ... 
    '039', '040', '041', '042', '043', ... 
    '045', '046', '047', '048', '049', ...
    '050', '051', '052', '053', '054', ... 
    '055', '056', '057', '058', '060', ... 
    '061', '063', '064', '065', '066', ... 
    '067', '068', '069', '071', '072', ...
    '073', '074', '076', '077', '078', ...
    '079', '080'}; 
numSubs = length(subjects); 
session = {'bromo', 'placebo', 'tolcapone'}; 
ROIs = {'FPI', 'MFG', 'cMFG'};
numROIs = length(ROIs); 

%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
N = length(subjects); 

%output directory 
output_dir = [dcmbase, '/Second_Level/']; 



% ------------------------------------------------------------------------------
% Create a GDCMs
% i.e., a group of reduced DCMs generate at first level using BMR
% ------------------------------------------------------------------------------
% load subject 1 to find number of models
%dcmdir = [datadir,metadata.ParticipantID{1},preprostr,WhichNoise,'spDCM/',WhichSeed,'_',num2str(WhichStri),'_',WhichHemi,extraStr,'/'];
%load([dcmdir,'DCM.mat'],'RCM');
load([dcmbase, 'DCM_models/DCM_sub-003models1-18.mat'], 'RCM'); 

numModels = length(RCM); 
clear RCM 

% Load DCMs for each subject into a cell array
GRCM = cell(0);
GBMC.F = zeros(numSubs,numModels);
GBMC.P = zeros(numSubs,numModels);


fprintf(1, 'Loading DCMs...\n');

%% probably DON'T need this loop dav 15 feb 2019 


    
for i = 1:length(subjects)
    subject=subjects(i);
    sub = char(subject); 
    %sess = session(s); 
    dir = fullfile(['sub-', sub]); 

    load([dcmbase, 'DCM_models/DCM_sub-', sub, 'models1-18.mat'], 'RCM', 'BMC');


    GRCM(i,:) = RCM;
    % GRCM(i,:) = RCM(1,1:3);
    GBMC.F(i,:) = BMC.F;
    GBMC.P(i,:) = BMC.P;
    clear RCM BMC

end

%Temporary code REMOVE later 
%GRCM = GRCM(:, 1:8)

fprintf(1, 'Done with loading DCMs \n');


% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';

% Specify design matrix for N subjects. It should start with a constant column
M.X = ones(N,1);

%M.Xnames = ['Condition', cond];

% Choose field
field = {'B'};

rng('default') % for reproducibility 
clear PEB BMA BMR
close all
%field = {'A'};
PEB = spm_dcm_peb(GRCM,M,field);
[BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,:));
%[BMA,BMR] = spm_dcm_bmc_peb(PEB(1), GRCM(1,:));
%
% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------

GRCMs = GRCM(:,1);
spm_dcm_peb_review(BMA,GRCMs)
% spm_dcm_peb_review(PEB,GRCMs)

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_output.mat'],'GRCMs','GRCM', 'PEB','BMA','BMR')
    
    % ------------------------------------------------------------------------------
    % Do comparison across the three drug sessions 
    % ------------------------------------------------------------------------------
%     cd(dcmbase)
%     X2 = [1 1; 1 -1]
%     bromo_dir = fullfile(['sub-00', int2str(sub), '/session-bromo']);
%     placebo_dir = fullfile(['sub-00', int2str(sub), '/session-placebo']);
%     cd(bromo_dir)
%     name_session1 = ['GCM_sub-00', int2str(sub), '-session-bromo.mat'];
%     cd(dcmbase)
%     cd(placebo_dir)
%     name_session2 = ['GCM_sub-00', int2str(sub), '-session-placebo.mat']; 
%     
%     
%     PEB_subject1 = spm_dcm_peb([name_session1; name_session2],X2)
    

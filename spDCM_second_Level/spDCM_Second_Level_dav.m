%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_fMRIprep/DCM_DA_team/';

%% Experiment information 
subjects = [1:4, 9:12,14,18:20, 22:26, 29:33, 35:43, 45:47, 51:55, 57:58, 60:61, 63:69, 72:73];
numSubs = length(subjects); 
session = {'bromo', 'placebo', 'tolcapone'}; 
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
numROIs = length(ROIs); 

%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
N = length(subjects); 

%output directory 
output_dir = [dcmbase, '/Second_Level/RH_second_spDCM']; 



% ------------------------------------------------------------------------------
% Create a GDCMs
% i.e., a group of reduced DCMs generate at first level using BMR
% ------------------------------------------------------------------------------
% load subject 1 to find number of models
%dcmdir = [datadir,metadata.ParticipantID{1},preprostr,WhichNoise,'spDCM/',WhichSeed,'_',num2str(WhichStri),'_',WhichHemi,extraStr,'/'];
%load([dcmdir,'DCM.mat'],'RCM');
load([dcmbase, 'RH_dcm_models/sub-001/session-bromo/DCM_sub-001-session-bromo.mat'], 'RCM'); 

numModels = length(RCM); 
clear RCM 

% Load DCMs for each subject into a cell array
GRCM = cell(0);
GBMC.F = zeros(numSubs,numModels);
GBMC.P = zeros(numSubs,numModels);


fprintf(1, 'Loading DCMs...\n');

for s = 1:length(session)
    sess = session(s); 
    cond=char(sess);

    
    for i = 1:length(subjects) 
        sub=subjects(i);
        %sess = session(s); 
        dir = fullfile(['sub-00', int2str(sub), '/session-', cond]); 
        
        if sub <= 9 
            load([dcmbase, 'RH_dcm_models/sub-00', int2str(sub), '/session-', cond, '/DCM_sub-00', int2str(sub), '-session-', cond, '.mat'], 'RCM', 'BMC');
        else 
            load([dcmbase, 'RH_dcm_models/sub-0', int2str(sub), '/session-', cond, '/DCM_sub-0', int2str(sub), '-session-', cond, '.mat'], 'RCM', 'BMC');
        end

        GRCM(i,:) = RCM;
        % GRCM(i,:) = RCM(1,1:3);
        GBMC.F(i,:) = BMC.F;
        GBMC.P(i,:) = BMC.P;
        clear RCM BMC

    end
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
    field = {'A'};

    rng('default') % for reproducibility
    clear PEB BMA BMR
    close all
    field = {'A'};
    PEB = spm_dcm_peb(GRCM,M,field);
    [BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,:));

    % ------------------------------------------------------------------------------
    % Review results
    % ------------------------------------------------------------------------------
    pause(2)
    GRCMs = GRCM(:,1);
    %spm_dcm_peb_review(BMA,GRCMs)
    % spm_dcm_peb_review(PEB,GRCMs)

    % ------------------------------------------------------------------------------
    % Save DCM outputs
    % ------------------------------------------------------------------------------
    cd(output_dir); 
    save(['PEB_output', cond, '.mat'],'GRCMs','GRCM', 'PEB','BMA','BMR')
    
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
    
end
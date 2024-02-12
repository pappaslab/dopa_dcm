%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; clc; close all; 

%% Directories
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
dcmbase = [base, 'DCM_models_only_winning_conn']; 
%output directory 
output_dir = [base, '/Second_Level/']; 

SPM_path = '/home/despoC/dvogel/Toolbox/spm12_v7487/';
addpath(genpath(SPM_path)); 

subFndr = dir([base '/sub*']);
numSub = length(subFndr);

%% Experiment information 
% subjects = {'001', '002', '003', '004', ...
%     '005'} 
%numSubs = length(subjects); 
sessions = {'bromo', 'placebo', 'tolcapone'}; 
ROIs = {'FPI', 'MFG', 'cMFG'};
numROIs = length(ROIs); 

%% Info about dataset DA_team
DA_team = 0; %if 1 then dopa team data if 0 then other dataset (e.g. TBP)
%% Info about dataset DA_team
if DA_team
    TR = 2;
    TE = 24; %TE is in ms 
    num_slices = 36; 
else
    TR = 1.6;
    TE = 27; 
    num_slices = 28;
end 


%load the DCMs per subject and for all three drug sessions     
% % for i = 1:length(subjects)
% %     subj=subjects(i);
% %     sub = char(subj); 
% %     subject = ['sub-' sub]; 
% %     %sess = session(s); 
% %     dir = fullfile(['sub-', sub]);
for subIdx = 1:numSub
    subject = subFndr(subIdx).name;    
    for j = 1:length(sessions)
        sess = sessions(j); 
        ses = char(sess); 
        session = ['session-' ses]; 
        cd(dcmbase); 
        
        clear S1
        %load each of the models 
        S1(j) = load(['/DCM_', subject, '-', session, '-models.mat'], 'DCM');
        %put it into a cell array so spm_dcm_peb can read it 
        GRCM{subIdx,j} = S1(j).DCM{1,1}; 
        

        
    end
end 

fprintf(1, 'Done with loading DCMs \n');

% % GRCM_1 = GRCM(:,1);
% % GRCM_2 = GRCM(:,2);
% % GRCM_3 = GRCM(:,3);
% % GRCM_new = [GRCM_1; GRCM_2; GRCM_3]; 


% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';
% for subs = 1:numSub
%     M.X(subs,1) = [1 1 1;1 -1 0;0 -1 1]; 
% end 
% Specify design matrix for N subjects. It should start with a constant column
%M.X{subIdx} = [1 1 1;1 -1 0;0 -1 1]; 
Xb = ones(numSub, 1); 
Xw = [1 1 1;1 -1 0;0 -1 1];
M.X = kron(Xb, Xw); 

%M.X(:,1) = ones(numSub,1)%[1 1 1;1 -1 0;0 -1 1]; 
% % M.X(:,2) = [ones(66, 1); -1*ones(66, 1); zeros(66,1)];
% % M.X(:,3) = [zeros(66, 1); -1*ones(66, 1); ones(66,1)];
%M.X(:,2) = -1*ones(numSub,1);
%M.X(:,3) = zeros(numSub,1);
% M.X(:,4) = ones(numSub,1);
% M.X(:,5) = ones(numSub,1);
% M.X(:,6) = ones(numSub,1);
% M.X(:,4) = ones(numSub,1); 
% M.X(:,5) = [ones(numSub,1)*-1];
% M.X(:,6) = zeros(numSub,1); 
%M.X(:,2) = [ones(numSub,1)*-1]; 
%M.X(:,3) = zeros(numSub,1); 
%M.X(:,1) = ones(numSub,1); 

%M.Xnew = repmat(M.X, numSub,1)
%M.X = [0 1 0; 0 0 0; 0 0 0]; 
%M.Xnames = {'group', 'brom_vs_plac', 'tolc_vs_plac'}; 
M.Xnames = {'Group', 'bromo_vs_placebo', 'tolc_vs_plac'}; 

% Choose field
field = {'A'};
rng('default') % for reproducibility 
cd(dcmbase); 

% RUN the first PEB (so second level PEB to compare the sessions for
% each subject 
PEB = spm_dcm_peb(GRCM,M,field);
%[BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,:));
[BMA BMR] = spm_dcm_peb_bmc(PEB(1)); 
% % PEBs{subIdx,1} = PEB %put it into a cell array 
% % 
% % %end
% % clear M 
% % 
% % % Specify PEB model settings (see batch editor for help on each setting)
% % M = struct();
% % M.alpha = 1;
% % M.beta  = 16;
% % M.hE    = 0;
% % M.hC    = 1/16;
% % M.Q     = 'single';
% % M.X = ones(numSub,1); 
% % M.Xnames = {'group'}; 
% % 
% % % Choose field
% % field = {'A'};
% % 
% % rng('default') % for reproducibility 
% % 
% % %Run the PEB of PEBs
% % PEB_group = spm_dcm_peb(PEBs, M, field);
% % 
% % BMA = spm_dcm_peb_bmc(PEB_group); 

 

% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------
%spm_dcm_peb_review(BMA, GRCM_new);
spm_dcm_peb_review(BMA, GRCM);

% ------------------------------------------------------------------------------
% Cross Validation procedure 
% ------------------------------------------------------------------------------
M.X = ones(numSub,1); 
spm_dcm_loo(GRCM,M,field); 

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_of_PEBS_output.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    


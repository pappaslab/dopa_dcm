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
% ROIs = {'FPI', 'MFG', 'cMFG'};
% numROIs = length(ROIs); 

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

% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';
% Specify design matrix
Xb = ones(numSub, 1); 
Xw = [1 1 1;1 -1 0;0 -1 1];
M.X = kron(Xb, Xw); 
M.Xnames = {'Group', 'bromo_vs_placebo', 'tolc_vs_plac'}; 

% Choose field
field = {'A'};
rng('default') % for reproducibility 
cd(dcmbase); 

% ------------------------------------------------------------------------------
% RUN PEB
% ------------------------------------------------------------------------------
PEB = spm_dcm_peb(GRCM,M,field);

% ------------------------------------------------------------------------------
% Run BMA 
% ------------------------------------------------------------------------------
%[BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,:));
[BMA BMR] = spm_dcm_peb_bmc(PEB(1)); 


% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------
%spm_dcm_peb_review(BMA, GRCM_new);
spm_dcm_peb_review(BMA, GRCM);

% ------------------------------------------------------------------------------
% Cross Validation procedure 
% ------------------------------------------------------------------------------
%M.X = ones(numSub,1); 
spm_dcm_loo(GRCM,M,field); 

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_of_PEBS_output.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    


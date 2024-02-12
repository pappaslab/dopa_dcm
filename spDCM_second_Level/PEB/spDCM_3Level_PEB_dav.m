%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang spring 2019 
%%-----------------------------------------------------------------------

clear all; 
clc; 
close all; 

%% Directories
base_root = '/home/despo/DCM_hierachy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';

base='/home/despo/DCM_hierachy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/DCM_models_only_winning_conn/';
%base = '/home/despo/DCM_hierachy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/afni_spDCM_task/';
dcmbase = [base];% 'DCM_models']; 
%dcmbase = [base, 'DCM_models_only_winning_conn']; 
%output directory 

base2='/home/despo/ioannis/';
output_dir = [base2, '/DCM_Second_Level_Nov2021/']; 
mkdir(output_dir);
SPM_path = '/home/despo/DCM_hierachy/Toolbox/spm12_v7487/';
addpath(genpath(SPM_path)); 

subFndr = dir([base_root '/sub*']);
numSub = length(subFndr);

%% Experiment information 
% % subjects = {'001', '002', '004', ...
% %     '013', '014', '037', '038', ...
% %     '039', '040', '041', '043', ...
% %     '045', '046', '048', '049', ...
% %     '053', '055', '057', '060', ...
% %     '064', '065', '068', '069', ...
% %     '073', '076', '078', '080'} 
% % numSub = length(subjects); 
sessions = {'bromo', 'placebo', 'tolcapone'};
%sessions = {'bromo', 'placebo'}; 


%% Info about dataset DA_team
DA_team = 1; %if 1 then dopa team data if 0 then other dataset (e.g. TBP)
%% Info about dataset DA_team
if DA_team
    TR = 2;
    TE = 24; %TE is in ms 
    num_slices = 36; 
else
    TR = 1.6;
    TE = 27
    num_slices = 28;
end 


for subIdx = 1:numSub
    subject = subFndr(subIdx).name; 
%     subj = subjects{subIdx};
%     subject = ['sub-' subj];
    for j = 1:length(sessions)
        sess = sessions(j); 
        ses = char(sess); 
        session = ['session-' ses]; 
        cd(dcmbase); 
        
        clear S1
        %load each of the models 
        S1(j) = load(['DCM_', subject, '-', session, '-models.mat'], 'DCM');
        %put it into a cell array so spm_dcm_peb can read it 
        GCM{j,1} = S1(j).DCM{1,1}; 

        
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
    %M.X = [1 1 1;1 -1 0;1 0 -1]'; 
    %M.X = [1 1 1;1 -1 0;0 -1 1; 1 0 -1]; 
    M.X =[1 1 0;1 -1 -1;1 0 1];
    %M.X = [1 1 1;1 -1 0;1 0 -1];
    %M.X = [1 1 1]; 
    M.Xnames = {'group', 'brom_vs_plac', 'tolc_vs_plac'}; 
    %M.Xnames = {'group'}; 

    % Choose field
    field = {'A'};
    rng('default') % for reproducibility 
    cd(dcmbase); 
    
    % RUN the first PEB (so second level PEB to compare the sessions for
    % each subject 
    PEB = spm_dcm_peb(GCM,M,field);
    PEBs{subIdx,1} = PEB; %put it into a cell array 

end
clear M 

% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';
M.X = ones(numSub,1); 
M.Xnames = {'group'}; 

% Choose field
field = {'A'};

rng('default') % for reproducibility 

% ------------------------------------------------------------------------------
% Run the PEB of PEBs
% ------------------------------------------------------------------------------
PEB_group = spm_dcm_peb(PEBs, M, field);


% ------------------------------------------------------------------------------
% Run Bayesian Model Comparison 
% ------------------------------------------------------------------------------
BMA = spm_dcm_peb_bmc(PEB_group); 


% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCM);

% ------------------------------------------------------------------------------
% Cross Validation procedure 
% ------------------------------------------------------------------------------
% M.X(:,1) = ones(numSub,1); 
% M.X(:,2) = ones(numSub,1); 
% M.X(34:66,2) = -1; 
% M.X(:,3) = ones(numSub,1);  
% M.X(23:44,3) = -1; 
% M.X(45:66,3) = 0; 
% 
% % Xb = ones(numSub, 1); 
% % Xw = [1 1 1;1 -1 0;0 -1 1];
% % M.X = kron(Xb, Xw); 
% spm_dcm_loo(PEBs,M,field); 

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_of_PEBS_Nove2021.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    


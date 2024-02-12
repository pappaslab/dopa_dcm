%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang spring 2019 
%%-----------------------------------------------------------------------

clear all; 
clc; 
close all; 

%% Directories
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DSST_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%dcmbase = [base, 'DCM_models_only_winning_conn']; 
dcmbase = [base, 'DCM_models']; 
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
%sessions = {'bromo', 'placebo'}; 
ROIs = {'FPI', 'MFG', 'cMFG'};
numROIs = length(ROIs); 

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
    M.X = [1 1 1;1 -1 0;0 -1 1]; 
    %M.X = [1 1; 1 -1]; 
    M.Xnames = {'group', 'brom_vs_plac', 'tolc_vs_plac'}; 
    %M.Xnames = {'group', 'brom_vs_plac'}; 

    % Choose field
    field = {'A'};
    rng('default') % for reproducibility 
    cd(dcmbase); 
    
    % RUN the first PEB (so second level PEB to compare the sessions for
    % each subject 
    PEB = spm_dcm_peb(GCM,M,field);
    PEBs{subIdx,1} = PEB %put it into a cell array 

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

%Run the PEB of PEBs
PEB_group = spm_dcm_peb(PEBs, M, field);

BMA = spm_dcm_peb_bmc(PEB_group); 

 

% ------------------------------------------------------------------------------
% Review results
% ------------------------------------------------------------------------------
spm_dcm_peb_review(BMA, GCM);

% ------------------------------------------------------------------------------
% Cross Validation procedure 
% ------------------------------------------------------------------------------
M.X = ones(numSub,1); 

% Xb = ones(numSub, 1); 
% Xw = [1 1 1;1 -1 0;0 -1 1];
% M.X = kron(Xb, Xw); 
spm_dcm_loo(PEBs,M,field); 

% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_of_PEBS_output_only_winning_conn.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    


%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% PEB of PEBS approach 
% D. A. Vogelsang 
%%-----------------------------------------------------------------------

clear all; 
clc; 
close all; 

%% Directories
%base = '/home/despo/DCM_hierachy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despo/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs_right_HEMISPERE';
base = '/home/despoC/dvogel/DCM_hierarchy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs_Schaeffer_ROIs'
%subject_base = '/home/despo/DCM_hierachy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despo/dvogel/NotBackedUp/spDCM_PEB_of_PEBs/';

dcmbase = [base, '/DCM_models']; 

%base2='/home/despo/ioannis/';
%output directory 
output_dir = [base, '/Second_Level_Nov2021/']; 
mkdir(output_dir)


% add SPM 
SPM_path = '/home/despo/DCM_hierachy/Toolbox/spm12_v7487/';
addpath(genpath(SPM_path)); 

%% Subject information 
subFndr = dir([base '/sub*']);
numSub = length(subFndr);

%% Experiment information 
sessions = {'bromo', 'placebo', 'tolcapone'};
%sessions = {'001', '002', '003'};


% loop to create a PEB for all sessions and subjects 
for subIdx = 1:numSub
    subject = subFndr(subIdx).name; 

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
    M.X = [1 1 0;1 0 0;1 0 1]; 
    %M.X = [1 0 0;1 1 0;1 0 1]; 
    %M.X = [1; 1; 1];
    M.Xnames = {'placebo', 'brom_vs_plac', 'tolc_vs_plac'}; 
    %M.Xnames = {'group', 'plac_vs_brom', 'tolc_vs_brom'}; 
    %M.Xnames = {'group'}
    

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
% Save DCM outputs
% ------------------------------------------------------------------------------
%cd(output_dir); 
%save(['PEBofPEBS_results.mat'],'GCM', 'PEBs','BMA', 'PEB_group'); 
    


%This is a script that concatenates sessions 
%It calls the function from Derek Nee collapseSessions
% Derek's function takes in an SPM design
% collapses multiple sessions together
% scans is a vector that specifies the scan lengths for the collapsed model
% For example, if SPM.nscan = [200 200 200 200]
% scans = [400 400] to collapse sessions 1 and 2 together and 3 and 4
% together
% scans = [800] would collapse all sessions into a single session
%
% Recommended use: first create and estimate a design matrix with each run 
% as a session as you would normally
% Then supply that SPM.mat along with scans
%
% Scaling and masking will be copied over to new design
% Run regressors will be created for each run of the original design
% Temporal non-sphericity correction and high-pass filter will be applied 
% for each run of the original design
% Edge artifact regressor will be created to deal with edge effects from
% high-pass filter
% For best results, runs should start and end with un-modeled fixation

clear all clc 

%base info 
%base = '/home/despoC/HierarchyThetaBeta/fMRI_Experiment/Session2_Baseline/fMRIprep/data/fmriprep/fmriprep'; 
%preproc_base = '/home/despoC/HierarchyThetaBeta/fMRI_Experiment/Session2_Baseline/fMRIprep/data/preproc';
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/stochDCM';

%design_matrix_base = '/home/despoC/HierarchyThetaBeta/fMRI_Experiment/Session2_Baseline/DesignMatrices'; 

collapseSession_path = '/home/despoC/dvogel/fMRI_script_Despolab/Concatenate_sessions/Denee_script/'; 
addpath(genpath(collapseSession_path)); 

%Experimental info 
subjects = {'001'}; 




%% RUN concatenate 
for i = 1:length(subjects) 

    sub=subjects{i};
    subbie = sub(2:3); 
    
    subs = ['sub-' sub];
    

    SPM = [dcmbase, '/', subs, '/SPM.mat']; 

    load(SPM); 
    scans = [540]; % 3 times 180 epi scans 
    
    SPM = collapseSessions(SPM,scans); 
    
    
    
    disp('---Finished concatenating sessions... now save the SPM---');
    
    save([first_level_base, '/MNI_space/', subs, '/SPM.mat'], 'SPM'); 
    
end 

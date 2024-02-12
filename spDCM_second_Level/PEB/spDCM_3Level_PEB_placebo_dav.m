%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; 
clc; 
close all; 

%% Directories
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DSST_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
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
%sessions = {'bromo', 'placebo', 'tolcapone'};
% sessions = {'placebo'}; 
% ROIs = {'FPI', 'MFG', 'cMFG'};
% numROIs = length(ROIs); 

%% Info about dataset DA_team
DA_team = 1; %if 1 then dopa team data if 0 then other dataset (e.g. TBP)
%% Info about dataset DA_team
% % if DA_team
% %     TR = 2;
% %     TE = 24; %TE is in ms 
% %     num_slices = 36; 
% % else
% %     TR = 1.6;
% %     TE = 27
% %     num_slices = 28;
% % end 


%load the DCMs per subject and for all three drug sessions     
% % for i = 1:length(subjects)
% %     subj=subjects(i);
% %     sub = char(subj); 
% %     subject = ['sub-' sub]; 
% %     %sess = session(s); 
% %     dir = fullfile(['sub-', sub]);
for subIdx = 1:numSub
    subject = subFndr(subIdx).name;    
%     for j = 1:length(sessions)
%         sess = sessions(j); 
%         ses = char(sess); 
%         session = ['session-' ses]; 
%         cd(dcmbase); 
%         
% %         clear S1
% %         %load each of the models 
% %         S1(j) = load(['/DCM_', subject, '-', session, '-models.mat'], 'RCM');
% %         %put it into a cell array so spm_dcm_peb can read it 
% %         GCM{j,1} = S1(j).DCM{1,1}; 
        cd(dcmbase); 
        load(['/DCM_', subject, '-session-placebo-models.mat'], 'DCM');
        %GRCM(subIdx,:) = RCM;
        GRCM{subIdx,:} = DCM{1,1};

        
    %end

% %     fprintf(1, 'Done with loading DCMs \n');
% % 
% % 
% %     % Specify PEB model settings (see batch editor for help on each setting)
% %     M = struct();
% %     M.alpha = 1;
% %     M.beta  = 16;
% %     M.hE    = 0;
% %     M.hC    = 1/16;
% %     M.Q     = 'single';
% %     % Specify design matrix for N subjects. It should start with a constant column
% %     M.X = [1 1 1;1 -1 0;0 -1 1]; 
% %     %M.X = [0 1 0; 0 0 0; 0 0 0]; 
% %     M.Xnames = {'group', 'brom_vs_plac', 'tolc_vs_plac'}; 
% %     %M.Xnames = {'placebo', 'plac', 'plac2'}; 
% % 
% %     % Choose field
% %     field = {'A'};
% %     rng('default') % for reproducibility 
% %     cd(dcmbase); 
% %     
% %     % RUN the first PEB (so second level PEB to compare the sessions for
% %     % each subject 
% %     PEB = spm_dcm_peb(GCM,M,field);
% %     PEBs{subIdx,1} = PEB %put it into a cell array 

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
M.X = ones(numSub,1);

%M.Xnames = ['Condition', cond];

% Choose field
field = {'A'};

rng('default') % for reproducibility 
clear PEB BMA BMR
close all
%field = {'A'};
PEB = spm_dcm_peb(GRCM,M,field);
%[BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,1));
[BMA,BMR] = spm_dcm_peb_bmc(PEB);
spm_dcm_peb_review(BMA, GRCM);


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
%spm_dcm_peb_review(BMA, GCM);


% ------------------------------------------------------------------------------
% Save DCM outputs
% ------------------------------------------------------------------------------
cd(output_dir); 
save(['PEB_placebo_output_cost18.mat'], 'PEB','BMA', 'GRCM'); 
    


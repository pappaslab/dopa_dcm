% This script computes the hierarchical strengths of the frontal ROIs 

%% D. A. Vogelsang October 2018
clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';

%% Experiment information 
%DA_team
subFndr = dir([dcmbase '/sub*']);
numSub = length(subFndr);

%subjects = [1]; 
%TBP subs
%subjects = [2,4:6,8,13:14,16,18:20, 23:24]
% numSubs = length(subjects); 
%session = {'bromo', 'placebo', 'tolcapone'};
session = {'placebo'}
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
numROIs = length(ROIs); 


for s = 1:length(session)
    sess = session(s); 
    cond=char(sess);
     for subj = 1:numSub%length(subjects)
         subject = subFndr(subj).name; 
            %sub=subjects(i); 
 
         %load([dcmbase, '/', subject, '/session-', cond, '/DCM_13-Jun-2019.mat'], 'DCM');
         load([dcmbase, 'DCM_models_4order/DCM_', subject, '-session-', cond, '-models.mat'], 'DCM');

                
                



                

                DCM_matrix = DCM{1,1}.Ep.A;
                %DCM_matrix = max(0,DCM_matrix);

                %FPI
                FPI_out(subj) = (DCM_matrix(1,1) + DCM_matrix(2,1) + DCM_matrix(3,1) + DCM_matrix(5,1)); 
                FPI_in(subj) = (DCM_matrix(1,1) + DCM_matrix(1,2) + DCM_matrix(1,3) + DCM_matrix(1,5)); 
                %FPI_out(i) = (DCM_matrix(1,1) + DCM_matrix(2,1) + DCM_matrix(3,1) + DCM_matrix(4,1) + DCM_matrix(5,1) + DCM_matrix(6,1)); 
                %FPI_in(i) = (DCM_matrix(1,1) + DCM_matrix(1,2) + DCM_matrix(1,3) + DCM_matrix(1,4) + DCM_matrix(1,5) + DCM_matrix(1,6)); 
                FPI_hier = FPI_out' - FPI_in';
                FPI_hier_mean = mean(FPI_hier); 

                %MFG
                MFG_out(subj) = (DCM_matrix(1,2) + DCM_matrix(2,2) + DCM_matrix(3,2) + DCM_matrix(5,2)); 
                MFG_in(subj) = (DCM_matrix(2,1) + DCM_matrix(2,2) + DCM_matrix(2,3)  + DCM_matrix(2,5)); 
                %MFG_out(sub) = (DCM_matrix(1,2) + DCM_matrix(2,2) + DCM_matrix(3,2) + DCM_matrix(4,2) + DCM_matrix(5,2) + DCM_matrix(6,2)); 
                %MFG_in(sub) = (DCM_matrix(2,1) + DCM_matrix(2,2) + DCM_matrix(2,3) + DCM_matrix(2,4) + DCM_matrix(2,5) + DCM_matrix(2,6)); 
                MFG_hier = MFG_out' - MFG_in';
                MFG_hier_mean = mean(MFG_hier); 

                %cMFG
                cMFG_out(subj) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3)); 
                cMFG_in(subj) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4)); 
                %cMFG_out(sub) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3) + DCM_matrix(5,3) + DCM_matrix(6,3)); 
                %cMFG_in(sub) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4) + DCM_matrix(3,5) + DCM_matrix(3,6)); 
                cMFG_hier = cMFG_out' - cMFG_in';
                cMFG_hier_mean = mean(cMFG_hier); 

                %SFS
                SFS_out(subj) = (DCM_matrix(3,4) + DCM_matrix(4,4)); 
                SFS_in(subj) = (DCM_matrix(4,3) + DCM_matrix(4,4)); 
                %SFS_out(sub) = (DCM_matrix(1,4) + DCM_matrix(2,4) + DCM_matrix(3,4) + DCM_matrix(4,4) + DCM_matrix(5,4) + DCM_matrix(6,4)); 
                %SFS_in(sub) = (DCM_matrix(4,1) + DCM_matrix(4,2) + DCM_matrix(4,3) + DCM_matrix(4,4) + DCM_matrix(4,5) + DCM_matrix(4,6)); 
                SFS_hier = SFS_out' - SFS_in';
                SFS_hier_mean = mean(SFS_hier); 

                %IFS
                IFS_out(subj) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(5,5) + DCM_matrix(6,5)); 
                IFS_in(subj) = (DCM_matrix(5,1) + DCM_matrix(5,2) + DCM_matrix(5,5) + DCM_matrix(5,6)); 
                %IFS_out(sub) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(3,5) + DCM_matrix(4,5) + DCM_matrix(5,5) + DCM_matrix(6,5)); 
                %IFS_in(sub) = (DCM_matrix(5,1) + DCM_matrix(5,2) + DCM_matrix(5,3) + DCM_matrix(5,4) + DCM_matrix(5,5) + DCM_matrix(5,6)); 
                IFS_hier = IFS_out' - IFS_in';
                IFS_hier_mean = mean(IFS_hier); 

                %IFJ
                IFJ_out(subj) = (DCM_matrix(5,6) + DCM_matrix(6,6)); 
                IFJ_in(subj) = (DCM_matrix(6,5) + DCM_matrix(6,6)); 
                %IFJ_out(sub) = (DCM_matrix(1,6) + DCM_matrix(2,6) + DCM_matrix(3,6) + DCM_matrix(4,6) + DCM_matrix(5,6) + DCM_matrix(6,6)); 
                %IFJ_in(sub) = (DCM_matrix(6,1) + DCM_matrix(6,2) + DCM_matrix(6,3) + DCM_matrix(6,4) + DCM_matrix(6,5) + DCM_matrix(6,6)); 
                IFJ_hier = IFJ_out' - IFJ_in';
                IFJ_hier_mean = mean(IFJ_hier); 
                
                
      
     end






end





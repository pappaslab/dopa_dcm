%%-----------------------------------------------------------------------
% SECOND level analysis for  spDCM 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

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
session = {'bromo', 'placebo', 'tolcapone'};
%session = {'placebo'}
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
numROIs = length(ROIs); 


for s = 1:length(session)
    sess = session(s); 
    cond=char(sess);
     for i = 1:numSub%length(subjects)
         subject = subFndr(i).name; 
            %sub=subjects(i); 
 
         load([dcmbase, '/', subject, '/session-', cond, '/DCM_13-Jun-2019.mat'], 'DCM');

        clear DCM_matrix
        DCM_matrix = DCM.Ep.A;
        %DCM_matrix = max(0,DCM_matrix);
        
        %TEST FPI to other regions and vice versa 
        FPI_to_MFG(i,s) =  DCM_matrix(2,1);
        FPI_to_cMFG(i,s) = DCM_matrix(3,1);
        %FPI_to_SFS(sub) = DCM_matrix(4,1);
        FPI_to_IFS(i,s) = DCM_matrix(5,1);
        %FPI_to_IFJ(sub) = DCM_matrix(6,1); 
        
        %t-test 
        [h,p_fpi_mfg] = ttest(FPI_to_MFG);
        [h,p_fpi_cmfg] = ttest(FPI_to_cMFG);
        %[h,p_fpi_sfs] = ttest(FPI_to_SFS);
        [h,p_fpi_ifs] = ttest(FPI_to_IFS);
        %[h,p_fpi_ifj] = ttest(FPI_to_IFJ);
        mean_FPI_to_MFG = mean(FPI_to_MFG);
        mean_FPI_to_cMFG = mean(FPI_to_cMFG); 
        mean_FPI_to_IFS = mean(FPI_to_IFS); 
     
        
        %FROM MFG to others 
        MFG_to_FPI (i,s) =  DCM_matrix(1,2);
        MFG_to_cMFG(i,s) = DCM_matrix(3,2);
        %MFG_to_SFS(i,s) = DCM_matrix(4,2);
        MFG_to_IFS(i,s) = DCM_matrix(5,2);
        %MFG_to_IFJ(sub) = DCM_matrix(6,2); 

        [h,p_mfg_fpi] = ttest(MFG_to_FPI);
        [h,p_mfg_cmfg] = ttest(MFG_to_cMFG);
        %[h,p_mfg_sfs] = ttest(MFG_to_SFS);
        [h,p_mfg_ifs] = ttest(MFG_to_IFS);
        %[h,p_mfg_ifj] = ttest(MFG_to_IFJ);
        mean_MFG_to_FPI = mean(MFG_to_FPI); 
        mean_MFG_to_cMFG = mean(MFG_to_cMFG); 
        %mean_MFG_to_SFS = mean(MFG_to_SFS); 
        mean_MFG_to_IFS = mean(MFG_to_IFS); 

        %FROM cMFG to others 
        cMFG_to_FPI(i,s) = DCM_matrix(1,3);
        cMFG_to_MFG(i,s) = DCM_matrix(2,3);
        cMFG_to_SFS(i,s) = DCM_matrix(4,3);
        %cMFG_to_IFS(i,s) = DCM_matrix(5,3);
        %cMFG_to_IFJ(sub) = DCM_matrix(6,3); 

        [h,p_cmfg_fpi] = ttest(cMFG_to_FPI);
        [h,p_cmfg_mfg] = ttest(cMFG_to_MFG);
        [h,p_cmfg_sfs] = ttest(cMFG_to_SFS);
        %[h,p_cmfg_ifs] = ttest(cMFG_to_IFS);
        %[h,p_cmfg_ifj] = ttest(cMFG_to_IFJ);
        mean_cMFG_to_FPI = mean(cMFG_to_FPI); 
        mean_cMFG_to_MFG = mean(cMFG_to_MFG); 
        mean_cMFG_to_SFS = mean(cMFG_to_SFS); 
        %mean_cMFG_to_IFS = mean(cMFG_to_IFS); 

        %FROM SFS to others 
        %SFS_to_FPI(i,s) = DCM_matrix(1,4);
        %SFS_to_MFG(i,s) = DCM_matrix(2,4);
        SFS_to_cMFG(i,s) = DCM_matrix(3,4);
        %SFS_to_IFS(sub) = DCM_matrix(5,4);
        %SFS_to_IFJ(i,s) = DCM_matrix(6,4); 

        %[h,p_sfs_fpi] = ttest(SFS_to_FPI);
        [h,p_sfs_mfg] = ttest(SFS_to_MFG);
        [h,p_sfs_cmfg] = ttest(SFS_to_cMFG);
        %[h,p_sfs_ifs] = ttest(SFS_to_IFS);
        [h,p_sfs_ifj] = ttest(SFS_to_IFJ);
        mean_SFS_to_MFG = mean(SFS_to_MFG); 
        mean_SFS_to_cMFG = mean(SFS_to_cMFG); 
        mean_SFS_to_IFJ = mean(SFS_to_IFJ); 

        %FROM IFS to others 
        IFS_to_FPI(i,s) = DCM_matrix(1,5);
        IFS_to_MFG(i,s) = DCM_matrix(2,5);
        IFS_to_cMFG(i,s) = DCM_matrix(3,5);
        %IFS_to_SFS(sub) = DCM_matrix(4,5);
        IFS_to_IFJ(i,s) = DCM_matrix(6,5); 

        [h,p_ifs_fpi] = ttest(IFS_to_FPI);
        [h,p_ifs_mfg] = ttest(IFS_to_MFG);
        [h,p_ifs_cmfg] = ttest(IFS_to_cMFG);
        %[h,p_ifs_sfs] = ttest(IFS_to_SFS);
        [h,p_ifs_ifj] = ttest(IFS_to_IFJ);
        mean_IFS_to_FPI = mean(IFS_to_FPI); 
        mean_IFS_to_MFG = mean(IFS_to_MFG); 
        mean_IFS_to_cMFG = mean(IFS_to_cMFG); 
        mean_IFS_to_IFJ = mean(IFS_to_IFJ); 
        

        %FROM IFS to others 
        %IFJ_to_FPI(sub) = DCM_matrix(1,6);
        %IFJ_to_MFG(sub) = DCM_matrix(2,6);
        %IFJ_to_cMFG(sub) = DCM_matrix(3,6);
        IFJ_to_SFS(i,s) = DCM_matrix(4,6);
        IFJ_to_IFS(i,s) = DCM_matrix(5,6); 

        %[h,p_ifj_fpi] = ttest(IFJ_to_FPI);
        %[h,p_ifj_mfg] = ttest(IFJ_to_MFG);
        %[h,p_ifj_cmfg] = ttest(IFJ_to_cMFG);
        [~,p_ifj_sfs] = ttest(IFJ_to_SFS);
        [h,p_ifj_ifs] = ttest(IFJ_to_IFS);
        mean_IFJ_to_SFS = mean(IFJ_to_SFS);
        mean_IFJ_to_IFS = mean(IFJ_to_IFS); 
        

     end
end



 
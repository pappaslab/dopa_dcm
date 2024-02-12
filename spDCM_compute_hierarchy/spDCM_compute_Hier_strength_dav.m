%Script to compute hierarchical strength
clear all 
clc

%base 
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/DCM_DA_team_rsfMRI_fmriprep_afni/';


GRCM_base = [base 'Second_Level/'];

sessions = {'bromo', 'placebo', 'tolcapone'}; 

% % subFndr = dir([base '/sub*']);
% % numSub = length(subFndr);

load([GRCM_base, 'PEB_placebo_output.mat'], 'GRCM'); 
numSub = length(GRCM(:,1)); 

for subIdx = 1:numSub
%     subject = subFndr(subIdx).name;
%     
%     DCM_file = [base, subject, '/session-placebo/GCM_', subject, '-session-placebo.mat']; 
%     if exist(DCM_file,'file')==2
%         
%     
%         load([base, subject, '/session-placebo/GCM_', subject, '-session-placebo.mat'], 'DCM'); 
%         GRCM(subIdx) = DCM{1,1}; 
%         RCM = DCM; 
    

        %DCM_matrix = RCM{1,1}.Ep.A; 
        DCM_matrix = GRCM{subIdx,1}.Ep.A;

        %TEST FPI to other regions and vice versa 
        FPI_to_MFG(subIdx) =  DCM_matrix(2,1);
        FPI_to_cMFG(subIdx) = DCM_matrix(3,1);
        FPI_to_SFS(subIdx) = DCM_matrix(4,1);
        FPI_to_IFS(subIdx) = DCM_matrix(5,1);
        FPI_to_IFJ(subIdx) = DCM_matrix(6,1); 


        %FROM MFG to others 
        MFG_to_FPI (subIdx) =  DCM_matrix(1,2);
        MFG_to_cMFG(subIdx) = DCM_matrix(3,2);
        MFG_to_SFS(subIdx) = DCM_matrix(4,2);
        MFG_to_IFS(subIdx) = DCM_matrix(5,2);
        MFG_to_IFJ(subIdx) = DCM_matrix(6,2); 


        %FROM cMFG to others 
        cMFG_to_FPI(subIdx) = DCM_matrix(1,3);
        cMFG_to_MFG(subIdx) = DCM_matrix(2,3);
        cMFG_to_SFS(subIdx) = DCM_matrix(4,3);
        cMFG_to_IFS(subIdx) = DCM_matrix(5,3);
        cMFG_to_IFJ(subIdx) = DCM_matrix(6,3); 


        %FROM SFS to others 
        SFS_to_FPI(subIdx) = DCM_matrix(1,4);
        SFS_to_MFG(subIdx) = DCM_matrix(2,4);
        SFS_to_cMFG(subIdx) = DCM_matrix(3,4);
        SFS_to_IFS(subIdx) = DCM_matrix(5,4);
        SFS_to_IFJ(subIdx) = DCM_matrix(6,4); 


        %FROM IFS to others 
        IFS_to_FPI(subIdx) = DCM_matrix(1,5);
        IFS_to_MFG(subIdx) = DCM_matrix(2,5);
        IFS_to_cMFG(subIdx) = DCM_matrix(3,5);
        IFS_to_SFS(subIdx) = DCM_matrix(4,5);
        IFS_to_IFJ(subIdx) = DCM_matrix(6,5); 


        %FROM IFS to others 
        IFJ_to_FPI(subIdx) = DCM_matrix(1,6);
        IFJ_to_MFG(subIdx) = DCM_matrix(2,6);
        IFJ_to_cMFG(subIdx) = DCM_matrix(3,6);
        IFJ_to_SFS(subIdx) = DCM_matrix(4,6);
        IFJ_to_IFS(subIdx) = DCM_matrix(5,6); 
% %     else 
% %         fprintf('Skipping this subject %d and %s\n\n',subject);
% %     end 


     %end 
end

[h,p_fpi_mfg] = ttest(FPI_to_MFG);
[h,p_fpi_cmfg] = ttest(FPI_to_cMFG);
[h,p_fpi_sfs] = ttest(FPI_to_SFS);
[h,p_fpi_ifs] = ttest(FPI_to_IFS);
[h,p_fpi_ifj] = ttest(FPI_to_IFJ);
mean_FPI_to_MFG = mean(FPI_to_MFG);

[h,p_mfg_fpi] = ttest(MFG_to_FPI);
[h,p_mfg_cmfg] = ttest(MFG_to_cMFG);
[h,p_mfg_sfs] = ttest(MFG_to_SFS);
[h,p_mfg_ifs] = ttest(MFG_to_IFS);
[h,p_mfg_ifj] = ttest(MFG_to_IFJ);


[h,p_cmfg_fpi] = ttest(cMFG_to_FPI);
[h,p_cmfg_mfg] = ttest(cMFG_to_MFG);
[h,p_cmfg_sfs] = ttest(cMFG_to_SFS);
[h,p_cmfg_ifs] = ttest(cMFG_to_IFS);
[h,p_cmfg_ifj] = ttest(cMFG_to_IFJ);

[h,p_sfs_fpi] = ttest(SFS_to_FPI);
[h,p_sfs_mfg] = ttest(SFS_to_MFG);
[h,p_sfs_cmfg] = ttest(SFS_to_cMFG);
[h,p_sfs_ifs] = ttest(SFS_to_IFS);
[h,p_sfs_ifj] = ttest(SFS_to_IFJ);

[h,p_ifs_fpi] = ttest(IFS_to_FPI);
[h,p_ifs_mfg] = ttest(IFS_to_MFG);
[h,p_ifs_cmfg] = ttest(IFS_to_cMFG);
[h,p_ifs_sfs] = ttest(IFS_to_SFS);
[h,p_ifs_ifj] = ttest(IFS_to_IFJ);

[h,p_ifj_fpi] = ttest(IFJ_to_FPI);
[h,p_ifj_mfg] = ttest(IFJ_to_MFG);
[h,p_ifj_cmfg] = ttest(IFJ_to_cMFG);
[~,p_ifj_sfs] = ttest(IFJ_to_SFS);
[h,p_ifj_ifs] = ttest(IFJ_to_IFS);

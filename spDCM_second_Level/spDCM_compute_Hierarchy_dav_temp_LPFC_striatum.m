% This script computes the hierarchical strengths of the frontal ROIs 

%% D. A. Vogelsang October 2018
clear all; clc; close all; 

%% Directories
%core_base
core_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/';  
%Datasets = {'DA_team', 'TBP', 'DSST', 'Dist_Exp', 'Idibrom'}; 
Datasets = {'DA_team', 'TBP'}; 
num_datasets = length(Datasets);

%DA_team
dateam_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/DCM_striatum_fmriprep_afni';
%TBP
tbp_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/DCM_striatum_fmriprep_afni';
%DSST
dsst_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DSST_rsfMRI/DCM_DSST_rsfMRI_fmriprep_afni';
%Dist_Exp
dist_exp_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/Dist_Exp_rsfMRI/DCM_Dist_Exp_rsfMRI_fmriprep_afni';
%Idibrom
idibrom_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/Idibrom_rsfMRI/DCM_Idibrom_rsfMRI_fmriprep_afni';

%% Experiment information 
%DA_team subs
%subjects = [1:5,9, 11:14,17:20, 22:26, 29:43, 45:58, 60:61, 63:69, 71:74, 76:80];
%subjects = [1:5,9, 11:14, 17:20, 22:26, 29:33, 35:43, 45:58, 60:61, 63:69, 71:74, 76:80]; 

%TBP subs
%subjects = [1:4, 6, 10, 12, 15:16, 18:24]
%DSST

%Dist_Exp
%subjects = [1:20]; 
%numSubs = length(subjects); 
session = {'bromo', 'placebo', 'tolcapone'}; 
%session = {'bromo', 'placebo'}; 
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'Caudate_Head', 'Caudate_body'};
numROIs = length(ROIs); 

%datasets = [1:6];


for d = 1:num_datasets
    
    currDataset = Datasets{d};
    
    %file = ([base, currDataset, '/FunImg/DCM_models']); 
    %cd(file)
    for s = 1:length(session)
        sess = session(s); 
        cond=char(sess);
            

            if d == 1%currDataset=='DA_team'; 
                 
                
            
                subjects_dateam = [1:5,9, 11, 13:14, 17:20, 22:26, 29:33, 35:43, 45:58, 60:61, 63:69, 71:74, 76:80]; %fmriprep_afni
                %subjects_dateam = [1, 3:9, 11, 13:14, 17:28, 30:39, 41:58, 60:61, 63:71, 73:74, 76:77]; %afni
                subjects = subjects_dateam; 

     

                
            elseif d == 2 %currDataset=='TBP'; 
                
                subjects_tbp = [1:2, 4:6, 8, 10, 12:24];
                subjects = subjects_tbp; 
            elseif d == 3 %currDataset=='DSST'; 
                
                subjects_dsst = [1:10, 12, 14:16, 18:24];
                subjects = subjects_dsst; 
            elseif d == 4 % currDataset=='Dist_Exp'; 
                
                subjects_dist_exp = [1:5, 9:12, 14:18, 20:21, 23]; 
                subjects = subjects_dist_exp; 
            elseif d == 5 %currDataset=='Idibrom'; 
                
                subjects_idibrom = [1:20]; 
                subjects = subjects_idibrom; 
            end 
                
                
            for i = 1:length(subjects) 
                sub=subjects(i);
                
                %load the dataset 
                if sub <= 9 
                   %DCM_model = strcat([core_base, currDataset, '_rsfMRI/DCM_', currDataset, '_rsfMRI_fmriprep_afni/sub-00', int2str(sub),'/session-', cond, '/DCM_sub-00', int2str(sub), '-session-', cond, '.mat']);
                   DCM_model = strcat([core_base, currDataset, '_rsfMRI/DCM_striatum_fmriprep_afni/sub-00', int2str(sub),'/session-', cond, '/DCM_sub-00', int2str(sub), '-session-', cond, '.mat']);
                else 
                   %DCM_model = strcat([core_base, currDataset, '_rsfMRI/DCM_', currDataset, '_rsfMRI_fmriprep_afni/sub-0', int2str(sub),'/session-', cond, '/DCM_sub-0', int2str(sub), '-session-', cond, '.mat']);
                   DCM_model = strcat([core_base, currDataset, '_rsfMRI/DCM_striatum_fmriprep_afni/sub-0', int2str(sub),'/session-', cond, '/DCM_sub-0', int2str(sub), '-session-', cond, '.mat']);
                end
                

                %sub = i;
                 
                if d == 1;
                    
                    subj = i; 

                elseif d == 2
                    if s == 1
                        subj = subj + 1;
                    elseif s == 2
                        subj = length(subjects_dateam) + i; 
                    elseif s == 3
                        subj = length(subjects_dateam) + i; 
                    end 
                elseif d == 3
                    if s == 1
                        subj = subj + 1;
                    elseif s == 2
                        subj = length(subjects_dateam) + length(subjects_tbp) + i; 
                    elseif s == 3
                        subj = length(subjects_dateam) + length(subjects_tbp) + i; 
                    end 
                elseif d == 4
                    if s == 1
                        subj = subj + 1;
                    elseif s == 2
                        subj = length(subjects_dateam) + length(subjects_tbp) + length(subjects_dsst) + i; 
                    elseif s == 3
                        subj = length(subjects_dateam) + length(subjects_tbp) + length(subjects_dsst) + i; 
                    end  
                elseif d == 5
                    if s == 1
                        subj = subj + 1;
                    elseif s == 2
                        subj = length(subjects_dateam) + length(subjects_tbp) + length(subjects_dsst) + length(subjects_dist_exp) + i; 
                    elseif s == 3
                        subj = length(subjects_dateam) + length(subjects_tbp) + length(subjects_dsst) + length(subjects_dist_exp) + i; 
                    end 



                end    


                load(DCM_model);

                DCM_matrix = DCM{1}.Ep.A;
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
                cMFG_out(subj) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3) + DCM_matrix(5,3)); 
                cMFG_in(subj) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4) + DCM_matrix(3,5)); 
                %cMFG_out(sub) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3) + DCM_matrix(5,3) + DCM_matrix(6,3)); 
                %cMFG_in(sub) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4) + DCM_matrix(3,5) + DCM_matrix(3,6)); 
                cMFG_hier = cMFG_out' - cMFG_in';
                cMFG_hier_mean = mean(cMFG_hier); 

                %SFS
                SFS_out(subj) = (DCM_matrix(3,4) + DCM_matrix(4,4) + DCM_matrix(6,4)); 
                SFS_in(subj) = (DCM_matrix(4,3) + DCM_matrix(4,4) + DCM_matrix(4,6)); 
                %SFS_out(sub) = (DCM_matrix(1,4) + DCM_matrix(2,4) + DCM_matrix(3,4) + DCM_matrix(4,4) + DCM_matrix(5,4) + DCM_matrix(6,4)); 
                %SFS_in(sub) = (DCM_matrix(4,1) + DCM_matrix(4,2) + DCM_matrix(4,3) + DCM_matrix(4,4) + DCM_matrix(4,5) + DCM_matrix(4,6)); 
                SFS_hier = SFS_out' - SFS_in';
                SFS_hier_mean = mean(SFS_hier); 

                %IFS
                IFS_out(subj) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(3,5) + DCM_matrix(6,5)); 
                IFS_in(subj) = (DCM_matrix(5,1) + DCM_matrix(5,2) + DCM_matrix(5,3) + DCM_matrix(5,6)); 
                %IFS_out(sub) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(3,5) + DCM_matrix(4,5) + DCM_matrix(5,5) + DCM_matrix(6,5)); 
                %IFS_in(sub) = (DCM_matrix(5,1) + DCM_matrix(5,2) + DCM_matrix(5,3) + DCM_matrix(5,4) + DCM_matrix(5,5) + DCM_matrix(5,6)); 
                IFS_hier = IFS_out' - IFS_in';
                IFS_hier_mean = mean(IFS_hier); 

                %IFJ
                IFJ_out(subj) = (DCM_matrix(4,6) + DCM_matrix(5,6) + DCM_matrix(6,6)); 
                IFJ_in(subj) = (DCM_matrix(6,4) + DCM_matrix(6,5) + DCM_matrix(6,6)); 
                %IFJ_out(sub) = (DCM_matrix(1,6) + DCM_matrix(2,6) + DCM_matrix(3,6) + DCM_matrix(4,6) + DCM_matrix(5,6) + DCM_matrix(6,6)); 
                %IFJ_in(sub) = (DCM_matrix(6,1) + DCM_matrix(6,2) + DCM_matrix(6,3) + DCM_matrix(6,4) + DCM_matrix(6,5) + DCM_matrix(6,6)); 
                IFJ_hier = IFJ_out' - IFJ_in';
                IFJ_hier_mean = mean(IFJ_hier); 
                
                %Compute strenghts of the connections 
                %TEST FPI to other regions and vice versa 


                if s == 1 
                    FPI_to_MFG_bromo(subj) =  DCM_matrix(2,1);
                    FPI_to_cMFG_bromo(subj) = DCM_matrix(3,1);
                    FPI_to_SFS_bromo(subj) = DCM_matrix(4,1);
                    FPI_to_IFS_bromo(subj) = DCM_matrix(5,1);
                    FPI_to_IFJ_bromo(subj) = DCM_matrix(6,1); 
                    FPI_to_CH_bromo(subj) = DCM_matrix(7,1); %CH = caudate head 
                    FPI_to_CB_bromo(subj) = DCM_matrix(8,1); %CB = caudate body 
                    [h,p_fpi_mfg_bromo] = ttest(FPI_to_MFG_bromo);
                    [h,p_fpi_cmfg_bromo] = ttest(FPI_to_cMFG_bromo);
                    [h,p_fpi_sfs_bromo] = ttest(FPI_to_SFS_bromo);
                    [h,p_fpi_ifs_bromo] = ttest(FPI_to_IFS_bromo);
                    [h,p_fpi_ifj_bromo] = ttest(FPI_to_IFJ_bromo);
                    [h,p_fpi_CH_bromo] = ttest(FPI_to_CH_bromo);
                    [h,p_fpi_CB_bromo] = ttest(FPI_to_CB_bromo);
                    %mean_FPI_to_MFG_bromo = mean(FPI_to_MFG);
                elseif s == 2
                    FPI_to_MFG_placebo(subj) =  DCM_matrix(2,1);
                    FPI_to_cMFG_placebo(subj) = DCM_matrix(3,1);
                    FPI_to_SFS_placebo(subj) = DCM_matrix(4,1);
                    FPI_to_IFS_placebo(subj) = DCM_matrix(5,1);
                    FPI_to_IFJ_placebo(subj) = DCM_matrix(6,1); 
                    FPI_to_CH_placebo(subj) = DCM_matrix(7,1); %CH = caudate head 
                    FPI_to_CB_placebo(subj) = DCM_matrix(8,1); %CB = caudate body 
                    [h,p_fpi_mfg_placebo] = ttest(FPI_to_MFG_placebo);
                    [h,p_fpi_cmfg_placebo] = ttest(FPI_to_cMFG_placebo);
                    [h,p_fpi_sfs_placebo] = ttest(FPI_to_SFS_placebo);
                    [h,p_fpi_ifs_placebo] = ttest(FPI_to_IFS_placebo);
                    [h,p_fpi_ifj_placebo] = ttest(FPI_to_IFJ_placebo);
                    [h,p_fpi_CH_placebo] = ttest(FPI_to_CH_placebo);
                    [h,p_fpi_CB_placebo] = ttest(FPI_to_CB_placebo);
                elseif s == 3
                    FPI_to_MFG_tolcapone(subj) =  DCM_matrix(2,1);
                    FPI_to_cMFG_tolcapone(subj) = DCM_matrix(3,1);
                    FPI_to_SFS_tolcapone(subj) = DCM_matrix(4,1);
                    FPI_to_IFS_tolcapone(subj) = DCM_matrix(5,1);
                    FPI_to_IFJ_tolcapone(subj) = DCM_matrix(6,1); 
                    FPI_to_CH_tolcapone(subj) = DCM_matrix(7,1); %CH = caudate head 
                    FPI_to_CB_tolcapone(subj) = DCM_matrix(8,1); %CB = caudate body 
                    [h,p_fpi_mfg_tolcapone] = ttest(FPI_to_MFG_tolcapone);
                    [h,p_fpi_cmfg_tolcapone] = ttest(FPI_to_cMFG_tolcapone);
                    [h,p_fpi_sfs_tolcapone] = ttest(FPI_to_SFS_tolcapone);
                    [h,p_fpi_ifs_tolcapone] = ttest(FPI_to_IFS_tolcapone);
                    [h,p_fpi_ifj_tolcapone] = ttest(FPI_to_IFJ_tolcapone);
                    [h,p_fpi_CH_tolcapone] = ttest(FPI_to_CH_tolcapone);
                    [h,p_fpi_CB_tolcapone] = ttest(FPI_to_CB_tolcapone);
                    
                end 
                    


                
                if s == 1 
                    MFG_to_FPI_bromo(subj) =  DCM_matrix(1,2);
                    MFG_to_cMFG_bromo(subj) = DCM_matrix(3,2);
                    MFG_to_SFS_bromo(subj) = DCM_matrix(4,2);
                    MFG_to_IFS_bromo(subj) = DCM_matrix(5,2);
                    MFG_to_IFJ_bromo(subj) = DCM_matrix(6,2); 
                    MFG_to_CH_bromo(subj) = DCM_matrix(7,2); %CH = caudate head 
                    MFG_to_CB_bromo(subj) = DCM_matrix(8,2); %CB = caudate body 

                    [h,p_mfg_fpi_bromo] = ttest(MFG_to_FPI_bromo);
                    [h,p_mfg_cmfg_bromo] = ttest(MFG_to_cMFG_bromo);
                    [h,p_mfg_sfs_bromo] = ttest(MFG_to_SFS_bromo);
                    [h,p_mfg_ifs_bromo] = ttest(MFG_to_IFS_bromo);
                    [h,p_mfg_ifj_bromo] = ttest(MFG_to_IFJ_bromo);
                    [h,p_mfg_CH_bromo] = ttest(MFG_to_CH_bromo);
                    [h,p_mfg_CB_bromo] = ttest(MFG_to_CB_bromo);
                elseif s == 2
                    MFG_to_FPI_placebo(subj) =  DCM_matrix(1,2);
                    MFG_to_cMFG_placebo(subj) = DCM_matrix(3,2);
                    MFG_to_SFS_placebo(subj) = DCM_matrix(4,2);
                    MFG_to_IFS_placebo(subj) = DCM_matrix(5,2);
                    MFG_to_IFJ_placebo(subj) = DCM_matrix(6,2); 
                    MFG_to_CH_placebo(subj) = DCM_matrix(7,2); %CH = caudate head 
                    MFG_to_CB_placebo(subj) = DCM_matrix(8,2); %CB = caudate body
                    
                    [h,p_mfg_fpi_placebo] = ttest(MFG_to_FPI_placebo);
                    [h,p_mfg_cmfg_placebo] = ttest(MFG_to_cMFG_placebo);
                    [h,p_mfg_sfs_placebo] = ttest(MFG_to_SFS_placebo);
                    [h,p_mfg_ifs_placebo] = ttest(MFG_to_IFS_placebo);
                    [h,p_mfg_ifj_placebo] = ttest(MFG_to_IFJ_placebo);
                    [h,p_mfg_CH_placebo] = ttest(MFG_to_CH_placebo);
                    [h,p_mfg_CB_placebo] = ttest(MFG_to_CB_placebo);
                    
                elseif s == 3
                    MFG_to_FPI_tolcapone(subj) =  DCM_matrix(1,2);
                    MFG_to_cMFG_tolcapone(subj) = DCM_matrix(3,2);
                    MFG_to_SFS_tolcapone(subj) = DCM_matrix(4,2);
                    MFG_to_IFS_tolcapone(subj) = DCM_matrix(5,2);
                    MFG_to_IFJ_tolcapone(subj) = DCM_matrix(6,2); 
                    MFG_to_CH_tolcapone(subj) = DCM_matrix(7,2); %CH = caudate head 
                    MFG_to_CB_tolcapone(subj) = DCM_matrix(8,2); %CB = caudate body
                    
                    [h,p_mfg_fpi_tolcapone] = ttest(MFG_to_FPI_tolcapone);
                    [h,p_mfg_cmfg_tolcapone] = ttest(MFG_to_cMFG_tolcapone);
                    [h,p_mfg_sfs_tolcapone] = ttest(MFG_to_SFS_tolcapone);
                    [h,p_mfg_ifs_tolcapone] = ttest(MFG_to_IFS_tolcapone);
                    [h,p_mfg_ifj_tolcapone] = ttest(MFG_to_IFJ_tolcapone);
                    [h,p_mfg_CH_tolcapone] = ttest(MFG_to_CH_tolcapone);
                    [h,p_mfg_CB_tolcapone] = ttest(MFG_to_CB_tolcapone);
                    
                end 


                
                if s == 1
                     %FROM cMFG to others 
                    cMFG_to_FPI_bromo(subj) = DCM_matrix(1,3);
                    cMFG_to_MFG_bromo(subj) = DCM_matrix(2,3);
                    cMFG_to_SFS_bromo(subj) = DCM_matrix(4,3);
                    cMFG_to_IFS_bromo(subj) = DCM_matrix(5,3);
                    cMFG_to_IFJ_bromo(subj) = DCM_matrix(6,3); 
                    cMFG_to_CH_bromo(subj) = DCM_matrix(7,3); %CH = caudate head 
                    cMFG_to_CB_bromo(subj) = DCM_matrix(8,3); %CB = caudate body 
                    

                    [h,p_cmfg_fpi_bromo] = ttest(cMFG_to_FPI_bromo);
                    [h,p_cmfg_mfg_bromo] = ttest(cMFG_to_MFG_bromo);
                    [h,p_cmfg_sfs_bromo] = ttest(cMFG_to_SFS_bromo);
                    [h,p_cmfg_ifs_bromo] = ttest(cMFG_to_IFS_bromo);
                    [h,p_cmfg_ifj_bromo] = ttest(cMFG_to_IFJ_bromo);
                    [h,p_cmfg_CH_bromo] = ttest(cMFG_to_CH_bromo);
                    [h,p_cmfg_CB_bromo] = ttest(cMFG_to_CB_bromo);
                    
                elseif s == 2
                    
                    cMFG_to_FPI_placebo(subj) = DCM_matrix(1,3);
                    cMFG_to_MFG_placebo(subj) = DCM_matrix(2,3);
                    cMFG_to_SFS_placebo(subj) = DCM_matrix(4,3);
                    cMFG_to_IFS_placebo(subj) = DCM_matrix(5,3);
                    cMFG_to_IFJ_placebo(subj) = DCM_matrix(6,3); 
                    cMFG_to_CH_placebo(subj) = DCM_matrix(7,3); %CH = caudate head 
                    cMFG_to_CB_placebo(subj) = DCM_matrix(8,3); %CB = caudate body 
                    
                    [h,p_cmfg_fpi_placebo] = ttest(cMFG_to_FPI_placebo);
                    [h,p_cmfg_mfg_placebo] = ttest(cMFG_to_MFG_placebo);
                    [h,p_cmfg_sfs_placebo] = ttest(cMFG_to_SFS_placebo);
                    [h,p_cmfg_ifs_placebo] = ttest(cMFG_to_IFS_placebo);
                    [h,p_cmfg_ifj_placebo] = ttest(cMFG_to_IFJ_placebo);
                    [h,p_cmfg_CH_placebo] = ttest(cMFG_to_CH_placebo);
                    [h,p_cmfg_CB_placebo] = ttest(cMFG_to_CB_placebo);
                    
                elseif s == 3
                    cMFG_to_FPI_tolcapone(subj) = DCM_matrix(1,3);
                    cMFG_to_MFG_tolcapone(subj) = DCM_matrix(2,3);
                    cMFG_to_SFS_tolcapone(subj) = DCM_matrix(4,3);
                    cMFG_to_IFS_tolcapone(subj) = DCM_matrix(5,3);
                    cMFG_to_IFJ_tolcapone(subj) = DCM_matrix(6,3); 
                    cMFG_to_CH_tolcapone(subj) = DCM_matrix(7,3); %CH = caudate head 
                    cMFG_to_CB_tolcapone(subj) = DCM_matrix(8,3); %CB = caudate body 
                    
                    [h,p_cmfg_fpi_tolcapone] = ttest(cMFG_to_FPI_tolcapone);
                    [h,p_cmfg_mfg_tolcapone] = ttest(cMFG_to_MFG_tolcapone);
                    [h,p_cmfg_sfs_tolcapone] = ttest(cMFG_to_SFS_tolcapone);
                    [h,p_cmfg_ifs_tolcapone] = ttest(cMFG_to_IFS_tolcapone);
                    [h,p_cmfg_ifj_tolcapone] = ttest(cMFG_to_IFJ_tolcapone);
                    [h,p_cmfg_CH_tolcapone] = ttest(cMFG_to_CH_tolcapone);
                    [h,p_cmfg_CB_tolcapone] = ttest(cMFG_to_CB_tolcapone);
                end 


                
                if s == 1 
                     %FROM SFS to others 
                    SFS_to_FPI_bromo(subj) = DCM_matrix(1,4);
                    SFS_to_MFG_bromo(subj) = DCM_matrix(2,4);
                    SFS_to_cMFG_bromo(subj) = DCM_matrix(3,4);
                    SFS_to_IFS_bromo(subj) = DCM_matrix(5,4);
                    SFS_to_IFJ_bromo(subj) = DCM_matrix(6,4); 
                    SFS_to_CH_bromo(subj) = DCM_matrix(7,4); %CH = caudate head 
                    SFS_to_CB_bromo(subj) = DCM_matrix(8,4); %CB = caudate body 

                    [h,p_sfs_fpi_bromo] = ttest(SFS_to_FPI_bromo);
                    [h,p_sfs_mfg_bromo] = ttest(SFS_to_MFG_bromo);
                    [h,p_sfs_cmfg_bromo] = ttest(SFS_to_cMFG_bromo);
                    [h,p_sfs_ifs_bromo] = ttest(SFS_to_IFS_bromo);
                    [h,p_sfs_ifj_bromo] = ttest(SFS_to_IFJ_bromo);
                    [h,p_sfs_CH_bromo] = ttest(SFS_to_CH_bromo);
                    [h,p_sfs_CB_bromo] = ttest(SFS_to_CB_bromo);
                    
                elseif s == 2
                    %FROM SFS to others 
                    SFS_to_FPI_placebo(subj) = DCM_matrix(1,4);
                    SFS_to_MFG_placebo(subj) = DCM_matrix(2,4);
                    SFS_to_cMFG_placebo(subj) = DCM_matrix(3,4);
                    SFS_to_IFS_placebo(subj) = DCM_matrix(5,4);
                    SFS_to_IFJ_placebo(subj) = DCM_matrix(6,4);
                    SFS_to_CH_placebo(subj) = DCM_matrix(7,4); %CH = caudate head 
                    SFS_to_CB_placebo(subj) = DCM_matrix(8,4); %CB = caudate body 
                    
                    [h,p_sfs_fpi_placebo] = ttest(SFS_to_FPI_placebo);
                    [h,p_sfs_mfg_placebo] = ttest(SFS_to_MFG_placebo);
                    [h,p_sfs_cmfg_placebo] = ttest(SFS_to_cMFG_placebo);
                    [h,p_sfs_ifs_placebo] = ttest(SFS_to_IFS_placebo);
                    [h,p_sfs_ifj_placebo] = ttest(SFS_to_IFJ_placebo);
                    [h,p_sfs_CH_placebo] = ttest(SFS_to_CH_placebo);
                    [h,p_sfs_CB_placebo] = ttest(SFS_to_CB_placebo);
                    
                elseif s == 3
                    SFS_to_FPI_tolcapone(subj) = DCM_matrix(1,4);
                    SFS_to_MFG_tolcapone(subj) = DCM_matrix(2,4);
                    SFS_to_cMFG_tolcapone(subj) = DCM_matrix(3,4);
                    SFS_to_IFS_tolcapone(subj) = DCM_matrix(5,4);
                    SFS_to_IFJ_tolcapone(subj) = DCM_matrix(6,4);
                    SFS_to_CH_tolcapone(subj) = DCM_matrix(7,4); %CH = caudate head 
                    SFS_to_CB_tolcapone(subj) = DCM_matrix(8,4); %CB = caudate body 
                    
                    [h,p_sfs_fpi_tolcapone] = ttest(SFS_to_FPI_tolcapone);
                    [h,p_sfs_mfg_tolcapone] = ttest(SFS_to_MFG_tolcapone);
                    [h,p_sfs_cmfg_tolcapone] = ttest(SFS_to_cMFG_tolcapone);
                    [h,p_sfs_ifs_tolcapone] = ttest(SFS_to_IFS_tolcapone);
                    [h,p_sfs_ifj_tolcapone] = ttest(SFS_to_IFJ_tolcapone);
                    [h,p_sfs_CH_tolcapone] = ttest(SFS_to_CH_tolcapone);
                    [h,p_sfs_CB_tolcapone] = ttest(SFS_to_CB_tolcapone);
                    
                end 


                
                if s == 1
                    %FROM IFS to others 
                    IFS_to_FPI_bromo(subj) = DCM_matrix(1,5);
                    IFS_to_MFG_bromo(subj) = DCM_matrix(2,5);
                    IFS_to_cMFG_bromo(subj) = DCM_matrix(3,5);
                    IFS_to_SFS_bromo(subj) = DCM_matrix(4,5);
                    IFS_to_IFJ_bromo(subj) = DCM_matrix(6,5); 
                    IFS_to_CH_bromo(subj) = DCM_matrix(7,5); %CH = caudate head 
                    IFS_to_CB_bromo(subj) = DCM_matrix(8,5); %CB = caudate body 

                    [h,p_ifs_fpi_bromo] = ttest(IFS_to_FPI_bromo);
                    [h,p_ifs_mfg_bromo] = ttest(IFS_to_MFG_bromo);
                    [h,p_ifs_cmfg_bromo] = ttest(IFS_to_cMFG_bromo);
                    [h,p_ifs_sfs_bromo] = ttest(IFS_to_SFS_bromo);
                    [h,p_ifs_ifj_bromo] = ttest(IFS_to_IFJ_bromo);
                    [h,p_ifs_CH_bromo] = ttest(IFS_to_CH_bromo);
                    [h,p_ifs_CB_bromo] = ttest(IFS_to_CB_bromo);
                elseif s == 2
                    %FROM IFS to others 
                    IFS_to_FPI_placebo(subj) = DCM_matrix(1,5);
                    IFS_to_MFG_placebo(subj) = DCM_matrix(2,5);
                    IFS_to_cMFG_placebo(subj) = DCM_matrix(3,5);
                    IFS_to_SFS_placebo(subj) = DCM_matrix(4,5);
                    IFS_to_IFJ_placebo(subj) = DCM_matrix(6,5); 
                    IFS_to_CH_placebo(subj) = DCM_matrix(7,5); %CH = caudate head 
                    IFS_to_CB_placebo(subj) = DCM_matrix(8,5); %CB = caudate body 
                    
                    [h,p_ifs_fpi_placebo] = ttest(IFS_to_FPI_placebo);
                    [h,p_ifs_mfg_placebo] = ttest(IFS_to_MFG_placebo);
                    [h,p_ifs_cmfg_placebo] = ttest(IFS_to_cMFG_placebo); 
                    [h,p_ifs_sfs_placebo] = ttest(IFS_to_SFS_placebo);
                    [h,p_ifs_ifj_placebo] = ttest(IFS_to_IFJ_placebo);
                    [h,p_ifs_CH_placebo] = ttest(IFS_to_CH_placebo);
                    [h,p_ifs_CB_placebo] = ttest(IFS_to_CB_placebo);
                elseif s == 3
                    IFS_to_FPI_tolcapone(subj) = DCM_matrix(1,5);
                    IFS_to_MFG_tolcapone(subj) = DCM_matrix(2,5);
                    IFS_to_cMFG_tolcapone(subj) = DCM_matrix(3,5);
                    IFS_to_SFS_tolcapone(subj) = DCM_matrix(4,5);
                    IFS_to_IFJ_tolcapone(subj) = DCM_matrix(6,5); 
                    IFS_to_CH_tolcapone(subj) = DCM_matrix(7,5); %CH = caudate head 
                    IFS_to_CB_tolcapone(subj) = DCM_matrix(8,5); %CB = caudate body 
                    
                    [h,p_ifs_fpi_tolcapone] = ttest(IFS_to_FPI_tolcapone);
                    [h,p_ifs_mfg_tolcapone] = ttest(IFS_to_MFG_tolcapone);
                    [h,p_ifs_cmfg_tolcapone] = ttest(IFS_to_cMFG_tolcapone); 
                    [h,p_ifs_sfs_tolcapone] = ttest(IFS_to_SFS_tolcapone);
                    [h,p_ifs_ifj_tolcapone] = ttest(IFS_to_IFJ_tolcapone);
                    [h,p_ifs_CH_tolcapone] = ttest(IFS_to_CH_tolcapone);
                    [h,p_ifs_CB_tolcapone] = ttest(IFS_to_CB_tolcapone);
                    
                end 


                
                if s == 1
                     %FROM IFS to others 
                    IFJ_to_FPI_bromo(subj) = DCM_matrix(1,6);
                    IFJ_to_MFG_bromo(subj) = DCM_matrix(2,6);
                    IFJ_to_cMFG_bromo(subj) = DCM_matrix(3,6);
                    IFJ_to_SFS_bromo(subj) = DCM_matrix(4,6);
                    IFJ_to_IFS_bromo(subj) = DCM_matrix(5,6); 
                    IFJ_to_CH_bromo(subj) = DCM_matrix(7,6); %CH = caudate head 
                    IFJ_to_CB_bromo(subj) = DCM_matrix(8,6); %CB = caudate body 

                    [h,p_ifj_fpi_bromo] = ttest(IFJ_to_FPI_bromo);
                    [h,p_ifj_mfg_bromo] = ttest(IFJ_to_MFG_bromo);
                    [h,p_ifj_cmfg_bromo] = ttest(IFJ_to_cMFG_bromo);
                    [~,p_ifj_sfs_bromo] = ttest(IFJ_to_SFS_bromo);
                    [h,p_ifj_ifs_bromo] = ttest(IFJ_to_IFS_bromo);
                    [h,p_ifj_CH_bromo] = ttest(IFJ_to_CH_bromo);
                    [h,p_ifj_CB_bromo] = ttest(IFJ_to_CB_bromo);
                    
                elseif s == 2 
                    %FROM IFS to others 
                    IFJ_to_FPI_placebo(subj) = DCM_matrix(1,6);
                    IFJ_to_MFG_placebo(subj) = DCM_matrix(2,6);
                    IFJ_to_cMFG_placebo(subj) = DCM_matrix(3,6);
                    IFJ_to_SFS_placebo(subj) = DCM_matrix(4,6);
                    IFJ_to_IFS_placebo(subj) = DCM_matrix(5,6); 
                    IFJ_to_CH_placebo(subj) = DCM_matrix(7,6); %CH = caudate head 
                    IFJ_to_CB_placebo(subj) = DCM_matrix(8,6); %CB = caudate body 
                    
                    [h,p_ifj_fpi_placebo] = ttest(IFJ_to_FPI_placebo);
                    [h,p_ifj_mfg_placebo] = ttest(IFJ_to_MFG_placebo);
                    [h,p_ifj_cmfg_placebo] = ttest(IFJ_to_cMFG_placebo);
                    [~,p_ifj_sfs_placebo] = ttest(IFJ_to_SFS_placebo);
                    [h,p_ifj_ifs_placebo] = ttest(IFJ_to_IFS_placebo);
                    [h,p_ifj_CH_placebo] = ttest(IFJ_to_CH_placebo);
                    [h,p_ifj_CB_placebo] = ttest(IFJ_to_CB_placebo);
                    
                elseif s ==3 
                    IFJ_to_FPI_tolcapone(subj) = DCM_matrix(1,6);
                    IFJ_to_MFG_tolcapone(subj) = DCM_matrix(2,6);
                    IFJ_to_cMFG_tolcapone(subj) = DCM_matrix(3,6);
                    IFJ_to_SFS_tolcapone(subj) = DCM_matrix(4,6);
                    IFJ_to_IFS_tolcapone(subj) = DCM_matrix(5,6); 
                    IFJ_to_CH_tolcapone(subj) = DCM_matrix(7,6); %CH = caudate head 
                    IFJ_to_CB_tolcapone(subj) = DCM_matrix(8,6); %CB = caudate body 
                    
                    
                    [h,p_ifj_fpi_tolcapone] = ttest(IFJ_to_FPI_tolcapone);
                    [h,p_ifj_mfg_tolcapone] = ttest(IFJ_to_MFG_tolcapone);
                    [h,p_ifj_cmfg_tolcapone] = ttest(IFJ_to_cMFG_tolcapone);
                    [~,p_ifj_sfs_tolcapone] = ttest(IFJ_to_SFS_tolcapone);
                    [h,p_ifj_ifs_tolcapone] = ttest(IFJ_to_IFS_tolcapone);
                    [h,p_ifj_CH_tolcapone] = ttest(IFJ_to_CH_tolcapone);
                    [h,p_ifj_CB_tolcapone] = ttest(IFJ_to_CB_tolcapone);
                    
                end
                
                
                if s == 1
                     %FROM CH to others 
                    CH_to_FPI_bromo(subj) = DCM_matrix(1,7);
                    CH_to_MFG_bromo(subj) = DCM_matrix(2,7);
                    CH_to_cMFG_bromo(subj) = DCM_matrix(3,7);
                    CH_to_SFS_bromo(subj) = DCM_matrix(4,7);
                    CH_to_IFS_bromo(subj) = DCM_matrix(5,7); 
                    CH_to_IFJ_bromo(subj) = DCM_matrix(6,7); %CH = caudate head 
                    CH_to_CB_bromo(subj) = DCM_matrix(8,7); %CB = caudate body 

                    [h,p_CH_fpi_bromo] = ttest(CH_to_FPI_bromo);
                    [h,p_CH_mfg_bromo] = ttest(CH_to_MFG_bromo);
                    [h,p_CH_cmfg_bromo] = ttest(CH_to_cMFG_bromo);
                    [~,p_CH_sfs_bromo] = ttest(CH_to_SFS_bromo);
                    [h,p_CH_ifs_bromo] = ttest(CH_to_IFS_bromo);
                    [h,p_CH_ifj_bromo] = ttest(CH_to_IFJ_bromo);
                    [h,p_CH_CB_bromo] = ttest(CH_to_CB_bromo);
                    
                elseif s == 2 
                    CH_to_FPI_placebo(subj) = DCM_matrix(1,7);
                    CH_to_MFG_placebo(subj) = DCM_matrix(2,7);
                    CH_to_cMFG_placebo(subj) = DCM_matrix(3,7);
                    CH_to_SFS_placebo(subj) = DCM_matrix(4,7);
                    CH_to_IFS_placebo(subj) = DCM_matrix(5,7); 
                    CH_to_IFJ_placebo(subj) = DCM_matrix(6,7); %CH = caudate head 
                    CH_to_CB_placebo(subj) = DCM_matrix(8,7); %CB = caudate body 

                    [h,p_CH_fpi_placebo] = ttest(CH_to_FPI_placebo);
                    [h,p_CH_mfg_placebo] = ttest(CH_to_MFG_placebo);
                    [h,p_CH_cmfg_placebo] = ttest(CH_to_cMFG_placebo);
                    [~,p_CH_sfs_placebo] = ttest(CH_to_SFS_placebo);
                    [h,p_CH_ifs_placebo] = ttest(CH_to_IFS_placebo);
                    [h,p_CH_ifj_placebo] = ttest(CH_to_IFJ_placebo);
                    [h,p_CH_CB_placebo] = ttest(CH_to_CB_placebo);
                    
                elseif s ==3 
                    CH_to_FPI_tolcapone(subj) = DCM_matrix(1,7);
                    CH_to_MFG_tolcapone(subj) = DCM_matrix(2,7);
                    CH_to_cMFG_tolcapone(subj) = DCM_matrix(3,7);
                    CH_to_SFS_tolcapone(subj) = DCM_matrix(4,7);
                    CH_to_IFS_tolcapone(subj) = DCM_matrix(5,7); 
                    CH_to_IFJ_tolcapone(subj) = DCM_matrix(6,7); %CH = caudate head 
                    CH_to_CB_tolcapone(subj) = DCM_matrix(8,7); %CB = caudate body 

                    [h,p_CH_fpi_tolcapone] = ttest(CH_to_FPI_tolcapone);
                    [h,p_CH_mfg_tolcapone] = ttest(CH_to_MFG_tolcapone);
                    [h,p_CH_cmfg_tolcapone] = ttest(CH_to_cMFG_tolcapone);
                    [~,p_CH_sfs_tolcapone] = ttest(CH_to_SFS_tolcapone);
                    [h,p_CH_ifs_tolcapone] = ttest(CH_to_IFS_tolcapone);
                    [h,p_CH_ifj_tolcapone] = ttest(CH_to_IFJ_tolcapone);
                    [h,p_CH_CB_tolcapone] = ttest(CH_to_CB_tolcapone);
                    
                end 
                
                if s == 1
                     %FROM CB to others 
                    CB_to_FPI_bromo(subj) = DCM_matrix(1,8);
                    CB_to_MFG_bromo(subj) = DCM_matrix(2,8);
                    CB_to_cMFG_bromo(subj) = DCM_matrix(3,8);
                    CB_to_SFS_bromo(subj) = DCM_matrix(4,8);
                    CB_to_IFS_bromo(subj) = DCM_matrix(5,8); 
                    CB_to_IFJ_bromo(subj) = DCM_matrix(6,8); %CH = caudate head 
                    CB_to_CH_bromo(subj) = DCM_matrix(7,8); %CB = caudate body 

                    [h,p_CB_fpi_bromo] = ttest(CB_to_FPI_bromo);
                    [h,p_CB_mfg_bromo] = ttest(CB_to_MFG_bromo);
                    [h,p_CB_cmfg_bromo] = ttest(CB_to_cMFG_bromo);
                    [~,p_CB_sfs_bromo] = ttest(CB_to_SFS_bromo);
                    [h,p_CB_ifs_bromo] = ttest(CB_to_IFS_bromo);
                    [h,p_CB_ifj_bromo] = ttest(CB_to_IFJ_bromo);
                    [h,p_CB_CH_bromo] = ttest(CB_to_CH_bromo);
                    
                elseif s == 2 
                    CB_to_FPI_placebo(subj) = DCM_matrix(1,8);
                    CB_to_MFG_placebo(subj) = DCM_matrix(2,8);
                    CB_to_cMFG_placebo(subj) = DCM_matrix(3,8);
                    CB_to_SFS_placebo(subj) = DCM_matrix(4,8);
                    CB_to_IFS_placebo(subj) = DCM_matrix(5,8); 
                    CB_to_IFJ_placebo(subj) = DCM_matrix(6,8); %CH = caudate head 
                    CB_to_CH_placebo(subj) = DCM_matrix(7,8); %CB = caudate body 

                    [h,p_CB_fpi_placebo] = ttest(CB_to_FPI_placebo);
                    [h,p_CB_mfg_placebo] = ttest(CB_to_MFG_placebo);
                    [h,p_CB_cmfg_placebo] = ttest(CB_to_cMFG_placebo);
                    [~,p_CB_sfs_placebo] = ttest(CB_to_SFS_placebo);
                    [h,p_CB_ifs_placebo] = ttest(CB_to_IFS_placebo);
                    [h,p_CB_ifj_placebo] = ttest(CB_to_IFJ_placebo);
                    [h,p_CB_CH_placebo] = ttest(CB_to_CH_placebo);
                    
                elseif s ==3 
                    CB_to_FPI_tolcapone(subj) = DCM_matrix(1,8);
                    CB_to_MFG_tolcapone(subj) = DCM_matrix(2,8);
                    CB_to_cMFG_tolcapone(subj) = DCM_matrix(3,8);
                    CB_to_SFS_tolcapone(subj) = DCM_matrix(4,8);
                    CB_to_IFS_tolcapone(subj) = DCM_matrix(5,8); 
                    CB_to_IFJ_tolcapone(subj) = DCM_matrix(6,8); %CH = caudate head 
                    CB_to_CH_tolcapone(subj) = DCM_matrix(7,8); %CB = caudate body 

                    [h,p_CB_fpi_tolcapone] = ttest(CB_to_FPI_tolcapone);
                    [h,p_CB_mfg_tolcapone] = ttest(CB_to_MFG_tolcapone);
                    [h,p_CB_cmfg_tolcapone] = ttest(CB_to_cMFG_tolcapone);
                    [~,p_CB_sfs_tolcapone] = ttest(CB_to_SFS_tolcapone);
                    [h,p_CB_ifs_tolcapone] = ttest(CB_to_IFS_tolcapone);
                    [h,p_CB_ifj_tolcapone] = ttest(CB_to_IFJ_tolcapone);
                    [h,p_CB_CH_tolcapone] = ttest(CB_to_CH_tolcapone);
                    
                end 
                    



            end

            if s == 1
                %FPI
                FPI_hier_bromo = FPI_hier; 
                FPI_hier_bromo_mean = FPI_hier_mean;
                [H P_fpi_hier_bromo] = ttest(FPI_hier_bromo);
                %MFG
                MFG_hier_bromo = MFG_hier; 
                MFG_hier_bromo_mean = MFG_hier_mean;
                [H P_mfg_hier_bromo] = ttest(MFG_hier_bromo);
                %cMFG
                cMFG_hier_bromo = cMFG_hier; 
                cMFG_hier_bromo_mean = cMFG_hier_mean;
                [H P_cmfg_hier_bromo] = ttest(cMFG_hier_bromo);
                %SFS
                SFS_hier_bromo = SFS_hier; 
                SFS_hier_bromo_mean = SFS_hier_mean;
                [H P_sfs_hier_bromo] = ttest(SFS_hier_bromo);
                %IFS
                IFS_hier_bromo = IFS_hier; 
                IFS_hier_bromo_mean = IFS_hier_mean;
                [H P_ifs_hier_bromo] = ttest(IFS_hier_bromo);
                %IFJ
                IFJ_hier_bromo = IFJ_hier; 
                IFJ_hier_bromo_mean = IFJ_hier_mean;
                [H P_ifj_hier_bromo] = ttest(IFJ_hier_bromo);

            elseif s == 2
                %FPI
                FPI_hier_placebo = FPI_hier; 
                FPI_hier_placebo_mean = FPI_hier_mean;
                [H P_fpi_hier_placebo] = ttest(FPI_hier_placebo);
                %MFG
                MFG_hier_placebo = MFG_hier; 
                MFG_hier_placebo_mean = MFG_hier_mean;
                [H P_mfg_hier_placebo] = ttest(MFG_hier_placebo);
                %cMFG
                cMFG_hier_placebo = cMFG_hier; 
                cMFG_hier_placebo_mean = cMFG_hier_mean;
                [~, P_cmfg_hier_placebo] = ttest(cMFG_hier_placebo);
                %SFS
                SFS_hier_placebo = SFS_hier; 
                SFS_hier_placebo_mean = SFS_hier_mean;
                [H P_sfs_hier_placebo] = ttest(SFS_hier_placebo);
                %IFS
                IFS_hier_placebo = IFS_hier; 
                IFS_hier_placebo_mean = IFS_hier_mean;
                [H P_ifs_hier_placebo] = ttest(IFS_hier_placebo);
                %IFJ
                IFJ_hier_placebo = IFJ_hier; 
                IFJ_hier_placebo_mean = IFJ_hier_mean;
                [H P_ifj_hier_placebo] = ttest(IFJ_hier_placebo);

            elseif s == 3
                %FPI
                FPI_hier_tolcapone = FPI_hier; 
                FPI_hier_tolcapone_mean = FPI_hier_mean;
                [H P_fpi_hier_tolcapone] = ttest(FPI_hier_tolcapone);
                %MFG
                MFG_hier_tolcapone = MFG_hier; 
                MFG_hier_tolcapone_mean = MFG_hier_mean;
                [H P_mfg_hier_tolcapone] = ttest(MFG_hier_tolcapone);
                %cMFG
                cMFG_hier_tolcapone = cMFG_hier; 
                cMFG_hier_tolcapone_mean = cMFG_hier_mean;
                [H P_cmfg_hier_tolcapone] = ttest(cMFG_hier_tolcapone);
                %SFS
                SFS_hier_tolcapone = SFS_hier; 
                SFS_hier_tolcapone_mean = SFS_hier_mean; 
                [H P_sfs_hier_tolcapone] = ttest(SFS_hier_tolcapone);
                %IFS
                IFS_hier_tolcapone = IFS_hier; 
                IFS_hier_tolcapone_mean = IFS_hier_mean;
                [H P_ifs_hier_tolcapone] = ttest(IFS_hier_tolcapone);
                %IFJ
                IFJ_hier_tolcapone = IFJ_hier; 
                IFJ_hier_tolcapone_mean = IFJ_hier_mean;
                [H P_ifj_hier_tolcapone] = ttest(IFJ_hier_tolcapone);
            end
            clear FPI_hier_mean MFG_hier_mean cMFG_hier_mean SFS_hier_mean IFS_hier_mean IFJ_hier_mean FPI_hier MFG_hier cMFG_hier SFS_hier IFS_hier IFJ_hier







    end
end 

%Overall comparisons
[H P_cmfg_plac_brom]= ttest(cMFG_hier_placebo, cMFG_hier_bromo); 


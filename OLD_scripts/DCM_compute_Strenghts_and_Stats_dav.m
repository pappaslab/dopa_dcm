%Script to compute hierarchical strength
clear all 

%I HAVE NOT CHECKED WHETHER THIS NEW SCRIPT WORKS!!!!!!! IT SHOULD BUT I AM
%NOT SURE 
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_DSST/'; 

Datasets = {'DCM_DSST', 'DCM_Dist_Exp', 'DCM_Idibrom', 'DCM_Guanfa_Bromo', 'DCM_TBP', 'DCM_DA_team'}; 
Datasets = {'DCM_DSST_Step1'}; 
num_datasets = length(Datasets);

%DSST
subjects_dsst = [1:23];
%Dist_Exp
subjects_dist_exp = [1:23];
%Idibrom
subjects_idibrom = [1:23]; 
%Guanfa_bromo
subjects_guanfa_bromo = [1:27];
%TBP
subjects_tbp = [1:24];
%DA_team
subjects_dateam = [1:39]; 


for d = 1:num_datasets
    currDataset = Datasets{d};
    
    file = ([base, currDataset, '/FunImg/DCM_models']); 
    cd(file)
    
    if d == 1;
        subjects = subjects_dsst;
    elseif d == 2;
        subjects = subjects_dist_exp; 
    elseif d == 3
        subjects = subjects_idibrom;
    elseif d == 4

        subjects = subjects_guanfa_bromo; 
    elseif d == 5
        subjects = subjects_tbp;
    elseif d == 6
        subjects = subjects_dateam; 
    end
     for i = 1:length(subjects)
            sub=subjects(i); 

  



            DCM_model = strcat(base,currDataset, '/FunImg/DCM_models/GCM_example.mat');
            load(DCM_model);

            DCM_matrix = GCM{sub,1}.Ep.A; 
            %remove all negative connections
            %DCM_matrix = max(0,DCM_matrix);

    % %         %now compute hierarchical strength for ingoing and outgoing per
    % %         %LPFC ROI
    % %         %FPI
    % %         FPI_out(sub) = (DCM_matrix(1,1) + DCM_matrix(2,1) + DCM_matrix(3,1) + DCM_matrix(5,1)); 
    % %         FPI_in(sub) = (DCM_matrix(1,1) + DCM_matrix(1,2) + DCM_matrix(1,3) + DCM_matrix(1,5)); 
    % %         FPI_hier = FPI_out' - FPI_in';
    % %         FPI_hier_mean = mean(FPI_hier); 
    % %         %MFG
    % %         MFG_out(sub) = (DCM_matrix(1,2) + DCM_matrix(2,2) + DCM_matrix(3,2) + DCM_matrix(5,2)); 
    % %         MFG_in(sub) = (DCM_matrix(2,1) + DCM_matrix(2,2) + DCM_matrix(2,3) + DCM_matrix(2,5)); 
    % %         MFG_hier = MFG_out' - MFG_in';
    % %         MFG_hier_mean = mean(MFG_hier); 
    % %         
    % %         %cMFG
    % %         cMFG_out(sub) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3)); 
    % %         cMFG_in(sub) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4)); 
    % %         cMFG_hier = cMFG_out' - cMFG_in';
    % %         cMFG_hier_mean = mean(cMFG_hier); 
    % %         
    % %         %SFS
    % %         SFS_out(sub) = (DCM_matrix(3,4) + DCM_matrix(4,4)); 
    % %         SFS_in(sub) = (DCM_matrix(4,3) + DCM_matrix(4,4)); 
    % %         SFS_hier = SFS_out' - SFS_in';
    % %         SFS_hier_mean = mean(SFS_hier); 
    % %         
    % %         %IFS
    % %         IFS_out(sub) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(6,5)); 
    % %         IFS_in(sub) = (DCM_matrix(5,1) + DCM_matrix(5,2)+ DCM_matrix(5,6)); 
    % %         IFS_hier = IFS_out' - IFS_in';
    % %         IFS_hier_mean = mean(IFS_hier); 
    % %         
    % %         %IFJ
    % %         IFJ_out(sub) = (DCM_matrix(5,6) + DCM_matrix(6,6)); 
    % %         IFJ_in(sub) = (DCM_matrix(6,5) + DCM_matrix(6,6)); 
    % %         IFJ_hier = IFJ_out' - IFJ_in';
    % %         IFJ_hier_mean = mean(IFJ_hier); 
            sub = i; 
            if d == 1;
                sub = sub + 0;
            elseif d == 2
                sub = sub + length(subjects_dsst);
            elseif d == 3
                sub = sub + length(subjects_dsst) + length(subjects_dist_exp);
            elseif d == 4
                sub = sub + length(subjects_dsst) + length(subjects_dist_exp) + length(subjects_idibrom); 
            elseif d == 5
                sub = sub + length(subjects_dsst) + length(subjects_dist_exp) + length(subjects_idibrom) + length(subjects_guanfa_bromo); 
            elseif d == 6
                sub = sub + length(subjects_dsst) + length(subjects_dist_exp) + length(subjects_idibrom) + length(subjects_guanfa_bromo) + length(subjects_tbp); 

            end
            %TEST FPI to other regions and vice versa 
            FPI_to_MFG(sub) =  DCM_matrix(2,1);
            FPI_to_cMFG(sub) = DCM_matrix(3,1);
            FPI_to_SFS(sub) = DCM_matrix(4,1);
            FPI_to_IFS(sub) = DCM_matrix(5,1);
            FPI_to_IFJ(sub) = DCM_matrix(6,1); 

            [h,p_fpi_mfg] = ttest(FPI_to_MFG);
            [h,p_fpi_cmfg] = ttest(FPI_to_cMFG);
            [h,p_fpi_sfs] = ttest(FPI_to_SFS);
            [h,p_fpi_ifs] = ttest(FPI_to_IFS);
            [h,p_fpi_ifj] = ttest(FPI_to_IFJ);
            mean_FPI_to_MFG = mean(FPI_to_MFG);

            %FROM MFG to others 
            MFG_to_FPI (sub) =  DCM_matrix(1,2);
            MFG_to_cMFG(sub) = DCM_matrix(3,2);
            MFG_to_SFS(sub) = DCM_matrix(4,2);
            MFG_to_IFS(sub) = DCM_matrix(5,2);
            MFG_to_IFJ(sub) = DCM_matrix(6,2); 

            [h,p_mfg_fpi] = ttest(MFG_to_FPI);
            [h,p_mfg_cmfg] = ttest(MFG_to_cMFG);
            [h,p_mfg_sfs] = ttest(MFG_to_SFS);
            [h,p_mfg_ifs] = ttest(MFG_to_IFS);
            [h,p_mfg_ifj] = ttest(MFG_to_IFJ);

            %FROM cMFG to others 
            cMFG_to_FPI(sub) = DCM_matrix(1,3);
            cMFG_to_MFG(sub) = DCM_matrix(2,3);
            cMFG_to_SFS(sub) = DCM_matrix(4,3);
            cMFG_to_IFS(sub) = DCM_matrix(5,3);
            cMFG_to_IFJ(sub) = DCM_matrix(6,3); 

            [h,p_cmfg_fpi] = ttest(cMFG_to_FPI);
            [h,p_cmfg_mfg] = ttest(cMFG_to_MFG);
            [h,p_cmfg_sfs] = ttest(cMFG_to_SFS);
            [h,p_cmfg_ifs] = ttest(cMFG_to_IFS);
            [h,p_cmfg_ifj] = ttest(cMFG_to_IFJ);

            %FROM SFS to others 
            SFS_to_FPI(sub) = DCM_matrix(1,4);
            SFS_to_MFG(sub) = DCM_matrix(2,4);
            SFS_to_cMFG(sub) = DCM_matrix(3,4);
            SFS_to_IFS(sub) = DCM_matrix(5,4);
            SFS_to_IFJ(sub) = DCM_matrix(6,4); 

            [h,p_sfs_fpi] = ttest(SFS_to_FPI);
            [h,p_sfs_mfg] = ttest(SFS_to_MFG);
            [h,p_sfs_cmfg] = ttest(SFS_to_cMFG);
            [h,p_sfs_ifs] = ttest(SFS_to_IFS);
            [h,p_sfs_ifj] = ttest(SFS_to_IFJ);

            %FROM IFS to others 
            IFS_to_FPI(sub) = DCM_matrix(1,5);
            IFS_to_MFG(sub) = DCM_matrix(2,5);
            IFS_to_cMFG(sub) = DCM_matrix(3,5);
            IFS_to_SFS(sub) = DCM_matrix(4,5);
            IFS_to_IFJ(sub) = DCM_matrix(6,5); 

            [h,p_ifs_fpi] = ttest(IFS_to_FPI);
            [h,p_ifs_mfg] = ttest(IFS_to_MFG);
            [h,p_ifs_cmfg] = ttest(IFS_to_cMFG);
            [h,p_ifs_sfs] = ttest(IFS_to_SFS);
            [h,p_ifs_ifj] = ttest(IFS_to_IFJ);

            %FROM IFS to others 
            IFJ_to_FPI(sub) = DCM_matrix(1,6);
            IFJ_to_MFG(sub) = DCM_matrix(2,6);
            IFJ_to_cMFG(sub) = DCM_matrix(3,6);
            IFJ_to_SFS(sub) = DCM_matrix(4,6);
            IFJ_to_IFS(sub) = DCM_matrix(5,6); 

            [h,p_ifj_fpi] = ttest(IFJ_to_FPI);
            [h,p_ifj_mfg] = ttest(IFJ_to_MFG);
            [h,p_ifj_cmfg] = ttest(IFJ_to_cMFG);
            [~,p_ifj_sfs] = ttest(IFJ_to_SFS);
            [h,p_ifj_ifs] = ttest(IFJ_to_IFS);

     end 
end

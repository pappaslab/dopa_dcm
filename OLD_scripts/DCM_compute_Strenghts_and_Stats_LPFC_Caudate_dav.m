%Script to compute hierarchical strength LPFC and caudate ROIs
clear all 

%I HAVE NOT CHECKED WHETHER THIS NEW SCRIPT WORKS!!!!!!! IT SHOULD BUT I AM
%NOT SURE 
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/'; 

Datasets = {'DCM_DSST', 'DCM_Dist_Exp', 'DCM_Idibrom', 'DCM_Guanfa_Bromo', 'DCM_TBP', 'DCM_DA_team'}; 
%Datasets = {'DCM_DSST', 'DCM_Dist_Exp','DCM_Idibrom', 'DCM_Guanfa_Bromo', 'DCM_DA_team'}; 
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
subjects_tbp = [1:25];
%DA_team
subjects_dateam = [1:39]; 


for d = 1:num_datasets
    currDataset = Datasets{d};
    
    file = ([base, currDataset, '/FunImg/DCM_models_LPFC_caudate']); 
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

  



            DCM_model = strcat(base,currDataset, '/FunImg/DCM_models_LPFC_caudate/GCM_example.mat');
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
            FPI_to_IFS(sub) = DCM_matrix(5,1);
            FPI_to_MCA(sub) = DCM_matrix(7,1);
            FPI_to_CCA(sub) = DCM_matrix(8,1);
            

            [h,p_fpi_mfg] = ttest(FPI_to_MFG);
            [h,p_fpi_cmfg] = ttest(FPI_to_cMFG);
            [h,p_fpi_ifs] = ttest(FPI_to_IFS);
            [h,p_fpi_mca] = ttest(FPI_to_MCA);
            [h,p_fpi_cca] = ttest(FPI_to_CCA);

            %FROM MFG to others 
            MFG_to_FPI (sub) =  DCM_matrix(1,2);
            MFG_to_cMFG(sub) = DCM_matrix(3,2);
            MFG_to_IFS(sub) = DCM_matrix(5,2);
            MFG_to_MCA(sub) = DCM_matrix(7,2);
            MFG_to_CCA(sub) = DCM_matrix(8,2);

            [h,p_mfg_fpi] = ttest(MFG_to_FPI);
            [h,p_mfg_cmfg] = ttest(MFG_to_cMFG);
            [h,p_mfg_ifs] = ttest(MFG_to_IFS);
            [h,p_mfg_mca] = ttest(MFG_to_MCA);
            [h,p_mfg_cca] = ttest(MFG_to_CCA);

            %FROM cMFG to others 
            cMFG_to_FPI(sub) = DCM_matrix(1,3);
            cMFG_to_MFG(sub) = DCM_matrix(2,3);
            cMFG_to_SFS(sub) = DCM_matrix(4,3);
            cMFG_to_IFS(sub) = DCM_matrix(5,3);
            cMFG_to_MCA(sub) = DCM_matrix(7,3);
            cMFG_to_CCA(sub) = DCM_matrix(8,3);

            [h,p_cmfg_fpi] = ttest(cMFG_to_FPI);
            [h,p_cmfg_mfg] = ttest(cMFG_to_MFG);
            [h,p_cmfg_sfs] = ttest(cMFG_to_SFS);
            [h,p_cmfg_ifs] = ttest(cMFG_to_IFS);
            [h,p_cmfg_mca] = ttest(cMFG_to_MCA);
            [h,p_cmfg_cca] = ttest(cMFG_to_CCA);
            
            %FROM SFS to others 
            SFS_to_cMFG(sub) = DCM_matrix(3,4);
            SFS_to_IFJ(sub) = DCM_matrix(6,4); 
            SFS_to_MCA(sub) = DCM_matrix(7,4); 
            SFS_to_CCA(sub) = DCM_matrix(8,4); 


            [h,p_sfs_cmfg] = ttest(SFS_to_cMFG);
            [h,p_sfs_ifj] = ttest(SFS_to_IFJ);
            [h,p_sfs_mca] = ttest(SFS_to_MCA);
            [h,p_sfs_cca] = ttest(SFS_to_CCA);
            
            %FROM IFS to others 
            IFS_to_FPI(sub) = DCM_matrix(1,5);
            IFS_to_MFG(sub) = DCM_matrix(2,5);
            IFS_to_cMFG(sub) = DCM_matrix(3,5);
            IFS_to_IFJ(sub) = DCM_matrix(6,5); 
            IFS_to_MCA(sub) = DCM_matrix(7,5); 
            IFS_to_CCA(sub) = DCM_matrix(8,5);             

            [h,p_ifs_fpi] = ttest(IFS_to_FPI);
            [h,p_ifs_mfg] = ttest(IFS_to_MFG);
            [h,p_ifs_cmfg] = ttest(IFS_to_cMFG);
            [h,p_ifs_ifj] = ttest(IFS_to_IFJ);
            [h,p_ifs_mca] = ttest(IFS_to_MCA);
            [h,p_ifs_cca] = ttest(IFS_to_CCA);

            %FROM IFJ to others 
            IFJ_to_SFS(sub) = DCM_matrix(4,6);
            IFJ_to_IFS(sub) = DCM_matrix(5,6); 
            IFJ_to_MCA(sub) = DCM_matrix(7,6);
            IFJ_to_CCA(sub) = DCM_matrix(8,6); 


            [~,p_ifj_sfs] = ttest(IFJ_to_SFS);
            [h,p_ifj_ifs] = ttest(IFJ_to_IFS);
            [h,p_ifj_mca] = ttest(IFJ_to_MCA);
            [h,p_ifj_cca] = ttest(IFJ_to_CCA);
            
            %FROM MCA to others 
            MCA_to_FPI(sub) = DCM_matrix(1,7);
            MCA_to_MFG(sub) = DCM_matrix(2,7);
            MCA_to_cMFG(sub) = DCM_matrix(3,7);
            MCA_to_IFS(sub) = DCM_matrix(5,7); 
            MCA_to_SFS(sub) = DCM_matrix(4,7); 
            MCA_to_IFJ(sub) = DCM_matrix(6,7);  
            MCA_to_CCA(sub) = DCM_matrix(8,7);      

            [h,p_mca_fpi] = ttest(MCA_to_FPI);
            [h,p_mca_mfg] = ttest(MCA_to_MFG);
            [h,p_mca_cmfg] = ttest(MCA_to_cMFG);
            [h,p_mca_ifs] = ttest(MCA_to_IFS);
            [h,p_mca_sfs] = ttest(MCA_to_SFS);
            [h,p_mca_ifj] = ttest(MCA_to_IFJ);
            [h,p_mca_cca] = ttest(MCA_to_CCA);
            
            %FROM CCA to others 
            CCA_to_FPI(sub) = DCM_matrix(1,8);
            CCA_to_MFG(sub) = DCM_matrix(2,8);
            CCA_to_cMFG(sub) = DCM_matrix(3,8);
            CCA_to_IFS(sub) = DCM_matrix(5,8); 
            CCA_to_SFS(sub) = DCM_matrix(4,8); 
            CCA_to_IFJ(sub) = DCM_matrix(6,8);  
            CCA_to_MCA(sub) = DCM_matrix(7,8);      

            [h,p_cca_fpi] = ttest(CCA_to_FPI);
            [h,p_cca_mfg] = ttest(CCA_to_MFG);
            [h,p_cca_cmfg] = ttest(CCA_to_cMFG);
            [h,p_cca_ifs] = ttest(CCA_to_IFS);
            [h,p_cca_sfs] = ttest(CCA_to_SFS);
            [h,p_cca_ifj] = ttest(CCA_to_IFJ);
            [h,p_cca_mca] = ttest(CCA_to_MCA);

     end 
end

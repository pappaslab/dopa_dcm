%Script to compute hierarchical strength

clear all 
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/FunImg/';

base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/'; %DCM_TBP/S2_FunImg/DCM_models';


Datasets = {'DCM_Idibrom', 'DCM_TBP', 'DCM_DA_team'}; 
num_datasets = length(Datasets);

%DSST
%subjects = [1:23];
%Idibrom GENETICS
%A1+ 
%subjects = [1,4:6, 9:10, 15:16, 18:20];
%A1-
%subjects = [2:3, 7:8, 11:14, 17, 21:23]; 

%TBP GENETICS
%A1+ 
%subjects = [2:3, 6:7, 10:14, 16, 18, 21, 23:24]; 
%A1-
subjects = [1,4:5, 8:9, 15, 17, 19:20, 22, 25];

%DA_team GENETICS
%A1+
%subjects = [2:3, 6, 8:10, 13, 15, 21:22, 24, 26, 29:34, 36, 39];
%A1-
%subjects = [1,4:5, 7, 11:12, 14, 16:20, 23, 25, 27:28, 35, 37:38]; 

for d = 1:num_datasets
    currDataset = Datasets{d};
    
    file = ([base, currDataset, '/FunImg/DCM_models']); 
    cd(file)
    
    if d == 1;
        %Idibrom GENETICS
        %A1+ 
        %subjects1 = [1,4:6, 9:10, 15:16, 18:20];
        %A1-
        subjects1 = [2:3, 7:8, 11:14, 17, 21:23]; 
        subjects = subjects1;
    elseif d == 2;
        
        %TBP GENETICS
        %A1+ 
        %subjects2 = [2:3, 6:7, 10:14, 16, 18, 21, 23:24]; 
        %A1-
        subjects2 = [1,4:5, 8:9, 15, 17, 19:20, 22, 25];
        subjects = subjects2; 
    elseif d == 3
        %DA_team GENETICS
        %A1+
        %subjects3 = [2:3, 6, 8:10, 13, 15, 21:22, 24, 26, 29:34, 36, 39];
        %A1-
        subjects3 = [1,4:5, 7, 11:12, 14, 16:20, 23, 25, 27:28, 35, 37:38]; 
        subjects = subjects3; 
    end
    
        

    % % %Idibrom
    % % subjects = [1:23]; 
    N = length(subjects); 



     for i = 1:length(subjects)
            sub2=subjects(i); 

    % %         if sub<=9
    % %             DCM_model = strcat(base,'/Sub_00',int2str(sub), 'DCM_model1.mat');
    % %         else
    % %             DCM_model = strcat(base,'/Sub_0',int2str(sub), 'DCM_model1.mat');
    % %         end

            DCM_model = strcat(file,'/GCM_example.mat');
            load(DCM_model);

            DCM_matrix = GCM{sub2,1}.Ep.A; 
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
                sub = sub + length(subjects1);
            elseif d == 3
                sub = sub + length(subjects1) + length(subjects2);
            end
            
            %sub = sub + 12;
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
cd([base, 'DCM_results_ALL']);

%Script to compute hierarchical strength

clear all 
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/FunImg/';

base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/FunImg/DCM_models1';

%DSST
subjects = [1:23];

% % %Idibrom
% % subjects = [1:23]; 
N = length(subjects); 



 for i = 1:length(subjects)
        sub=subjects(i); 
        
% %         if sub<=9
% %             DCM_model = strcat(base,'/Sub_00',int2str(sub), 'DCM_model1.mat');
% %         else
% %             DCM_model = strcat(base,'/Sub_0',int2str(sub), 'DCM_model1.mat');
% %         end

        DCM_model = strcat(base,'/GCM_example.mat');
        load(DCM_model);
        
        DCM_matrix = GCM{sub,1}.Ep.A; 
        %remove all negative connections
        %DCM_matrix = max(0,DCM_matrix);
        
        %now compute hierarchical strength for ingoing and outgoing per
        %LPFC ROI
        %FPI
        FPI_out(sub) = (DCM_matrix(1,1) + DCM_matrix(2,1) + DCM_matrix(3,1) + DCM_matrix(5,1)); 
        FPI_in(sub) = (DCM_matrix(1,1) + DCM_matrix(1,2) + DCM_matrix(1,3) + DCM_matrix(1,5)); 
        FPI_hier = FPI_out' - FPI_in';
        FPI_hier_mean = mean(FPI_hier); 
        %MFG
        MFG_out(sub) = (DCM_matrix(1,2) + DCM_matrix(2,2) + DCM_matrix(3,2) + DCM_matrix(5,2)); 
        MFG_in(sub) = (DCM_matrix(2,1) + DCM_matrix(2,2) + DCM_matrix(2,3) + DCM_matrix(2,5)); 
        MFG_hier = MFG_out' - MFG_in';
        MFG_hier_mean = mean(MFG_hier); 
        
        %cMFG
        cMFG_out(sub) = (DCM_matrix(1,3) + DCM_matrix(2,3) + DCM_matrix(3,3) + DCM_matrix(4,3)); 
        cMFG_in(sub) = (DCM_matrix(3,1) + DCM_matrix(3,2) + DCM_matrix(3,3) + DCM_matrix(3,4)); 
        cMFG_hier = cMFG_out' - cMFG_in';
        cMFG_hier_mean = mean(cMFG_hier); 
        
        %SFS
        SFS_out(sub) = (DCM_matrix(3,4) + DCM_matrix(4,4)); 
        SFS_in(sub) = (DCM_matrix(4,3) + DCM_matrix(4,4)); 
        SFS_hier = SFS_out' - SFS_in';
        SFS_hier_mean = mean(SFS_hier); 
        
        %IFS
        IFS_out(sub) = (DCM_matrix(1,5) + DCM_matrix(2,5) + DCM_matrix(6,5)); 
        IFS_in(sub) = (DCM_matrix(5,1) + DCM_matrix(5,2)+ DCM_matrix(5,6)); 
        IFS_hier = IFS_out' - IFS_in';
        IFS_hier_mean = mean(IFS_hier); 
        
        %IFJ
        IFJ_out(sub) = (DCM_matrix(5,6) + DCM_matrix(6,6)); 
        IFJ_in(sub) = (DCM_matrix(6,5) + DCM_matrix(6,6)); 
        IFJ_hier = IFJ_out' - IFJ_in';
        IFJ_hier_mean = mean(IFJ_hier); 
        
        %TEST FPI to MFG and vice versa 
        FPI_to_MFG(sub) = DCM_matrix(2,1);
        MFG_to_FPI(sub) = DCM_matrix(1,2);
 end 
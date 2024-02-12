% This script computes the hierarchical strengths of the frontal ROIs 

%% D. A. Vogelsang October 2018
clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_fMRIprep/DCM_DA_team/LH_dcm_models';

%% Experiment information 
subjects = [1:4,9:12,14,18:20, 22:26, 29:33, 35:43, 45:47, 51:55, 57:58, 60:61, 63:69, 72:73];
numSubs = length(subjects); 
session = {'bromo', 'placebo', 'tolcapone'}; 
ROIs = {'FPI'; 'MFG'; 'cMFG'; 'SFS'; 'IFS'; 'IFJ'};
numROIs = length(ROIs); 

for s = 1:length(session)
    sess = session(s); 
    cond=char(sess);

    %load one subject 
    DCM_model_test = strcat([dcmbase,'/sub-001/session-', cond, '/DCM_sub-001-session-', cond, '.mat']);
    load(DCM_model_test);

    %DCM_matrix = DCM{1}.Ep.A; 
    %find coordinates of the connections 
    [row, col] = find(DCM{1,1}.a==1); %col is from which region, row is to which region 
    
    for idx = 1:length(row);
        %idx_matrix = leng(idx); 
        col_idx = col(idx);
        row_idx = row(idx);
        
        name_matrix = repmat(ROIs,6); 
        name_matrix= name_matrix(1:6, :); 

        roi1 = char(name_matrix(col_idx)); %from region
        roi2 = char(name_matrix(row_idx)); % to region 

        name_roi = (['from_',roi1, '_to_', roi2]); 
        all_name_rois{idx} = name_roi; 
        clear BMA BMC DCM RCM

        for i = 1:length(subjects) 
            sub=subjects(i);
            
            if sub <= 9 
               DCM_model = strcat([dcmbase,'/sub-00', int2str(sub),'/session-', cond, '/DCM_sub-00', int2str(sub), '-session-', cond, '.mat']);
            else 
               DCM_model = strcat([dcmbase,'/sub-0', int2str(sub),'/session-', cond, '/DCM_sub-0', int2str(sub), '-session-', cond, '.mat']);
            end
            

            load(DCM_model);

            DCM_matrix = DCM{1}.Ep.A;
            %A = DCM{1}.Ep.A; 
            %maxA(i) = max(max(abs(A - diag(diag(A)))));

            Hier_strength.ses{s}.roi{idx}.subjectid{sub} = name_roi;

            Hier_strength.ses{s}.roi{idx}.value{sub} =  DCM_matrix(row_idx,col_idx);
            hier_value(i) = Hier_strength.ses{s}.roi{idx}.value{sub}; 
            
%             if col_idx == 1 %so from FPI to other regions
%                 FPI_hier_out_strength = hier_value; 
%                 FPI_hier_out_strength_all = sum(FPI_hier_out_strength); 
%             end
%             if row_idx == 1 %from other regions to FPI 
%                 FPI_hier_in_strength = hier_value; 
%                 FPI_hier_in_strength_all = sum(FPI_hier_in_strength);  
%             end
%             
%             Hierarchy_FPI(s,i) = FPI_hier_out_strength_all - FPI_hier_in_strength_all;

        end
        %disp(maxA); 
        mean_hier_value = mean(hier_value); 
%         All_Hier_results = struct(...
%             'roi_name', name_roi, ...
%             'session', cond, ...
%             'value', hier_value,...
%             'mean_hierarchical_value', mean_hier_value); 
%         path=([dcmbase, '/Second_Level']); 
%         cd(path)
%         save(['Hierarchy_',name_roi, '_', cond],'-struct','All_Hier_results')
        
        %compute FPI hierarchy 
        if col_idx == 1 %so from FPI to other regions
            FPI_hier_out_strength(row_idx) = mean_hier_value; 
            FPI_hier_out_strength_all = sum(FPI_hier_out_strength); 
        end
        if row_idx == 1
            FPI_hier_in_strength(col_idx) = mean_hier_value; 
            FPI_hier_in_strength_all = sum(FPI_hier_in_strength);  
        end
        
        %compute MFG hierarchy 
        if col_idx == 2 %so from FPI to other regions
            MFG_hier_out_strength(row_idx) = mean_hier_value; 
            MFG_hier_out_strength_all = sum(MFG_hier_out_strength); 
        end
        if row_idx == 2
            MFG_hier_in_strength(col_idx) = mean_hier_value; 
            MFG_hier_in_strength_all = sum(MFG_hier_in_strength);  
        end
        
        %compute cMFG hierarchy 
        if col_idx == 3 %so from FPI to other regions
            cMFG_hier_out_strength(row_idx) = mean_hier_value; 
            cMFG_hier_out_strength_all = sum(cMFG_hier_out_strength); 
        end
        if row_idx == 3
            cMFG_hier_in_strength(col_idx) = mean_hier_value; 
            cMFG_hier_in_strength_all = sum(cMFG_hier_in_strength);  
        end
        
        %compute SFS hierarchy 
        if col_idx == 4 %so from FPI to other regions
            SFS_hier_out_strength(row_idx) = mean_hier_value; 
            SFS_hier_out_strength_all = sum(SFS_hier_out_strength); 
        end
        if row_idx == 4
            SFS_hier_in_strength(col_idx) = mean_hier_value; 
            SFS_hier_in_strength_all = sum(SFS_hier_in_strength);  
        end
        
        %compute IFS hierarchy 
        if col_idx == 5 %so from  to other regions
            IFS_hier_out_strength(row_idx) = mean_hier_value; 
            IFS_hier_out_strength_all = sum(IFS_hier_out_strength); 
        end
        if row_idx == 5
            IFS_hier_in_strength(col_idx) = mean_hier_value; 
            IFS_hier_in_strength_all = sum(IFS_hier_in_strength);  
        end
        
        %compute IFJ hierarchy 
        if col_idx == 6 %so from FPI to other regions
            IFJ_hier_out_strength(row_idx) = mean_hier_value; 
            IFJ_hier_out_strength_all = sum(IFJ_hier_out_strength); 
        end
        if row_idx == 6
            IFJ_hier_in_strength(col_idx) = mean_hier_value; 
            IFJ_hier_in_strength_all = sum(IFJ_hier_in_strength);  
        end
        


    end

    compute hierarchy for all ROIs and this is done per drug session 
    Hierarchy_FPI(s) = FPI_hier_out_strength_all - FPI_hier_in_strength_all; %hierarchy per drug session 
    Hierarchy_MFG(s) = MFG_hier_out_strength_all - MFG_hier_in_strength_all;
    Hierarchy_cMFG(s) = cMFG_hier_out_strength_all - cMFG_hier_in_strength_all;
    Hierarchy_SFS(s) = SFS_hier_out_strength_all - SFS_hier_in_strength_all;
    Hierarchy_IFS(s) = IFS_hier_out_strength_all - IFS_hier_in_strength_all;
    Hierarchy_IFJ(s) = IFJ_hier_out_strength_all - IFJ_hier_in_strength_all;
end

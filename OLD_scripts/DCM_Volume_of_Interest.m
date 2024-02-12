%DCM for the Volume of Interst 

clear all 
% base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/';
% rsfMRI_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/';

base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_DSST/DCM_DSST_Step1/';
rsfMRI_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_DSST/DCM_DSST_Step1/';

subjects = [2:23];

fs = filesep;

%ROIs = {'PCC', 'mPFC', 'LIPC', 'RIPC'}; 
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'MCA', 'CCA'}; 
centers = [-44 48 4; -38 28 44; -26 14 52; -52 20 28; -40 10 20; -24 4 54];
%centers = [-44 48 4; -34 10 60; -26 14 52; -52 20 28; -40 10 20; -24 4 54; -16 -2 18; -11 13 7];


session = {'FunImg', 'S2_FunImg'};

for s = 1:length(session)
    sess = session(s); 

      for i = 1:length(subjects)


            sub=subjects(i);
            sess = char(sess); 

            if sub <= 9
                fmdir = strcat([rsfMRI_dir,sess, ['/Sub_00',int2str(sub)]]);
                %newDCM_dir = strcat(rsfMRI_dir,'/S2_FunImg/Sub_00',int2str(sub));
            else
                fmdir = strcat([rsfMRI_dir,sess, ['/Sub_0',int2str(sub)]]);
                %newDCM_dir = strcat(rsfMRI_dir,'/S2_FunImg/Sub_0',int2str(sub));;
            end

            cd(fmdir);
            clear r roi
            for r = 1:length(ROIs)
                roi = ROIs(r); 
                %for s = 1:length(session)
                matlabbatch{1}.spm.util.voi.spmmat = {'SPM.mat'};
                matlabbatch{1}.spm.util.voi.adjust = NaN;
                matlabbatch{1}.spm.util.voi.session = 1;
                matlabbatch{1}.spm.util.voi.name = char(roi);
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = centers(r,:);
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 6;
                matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'mask.nii,1'};
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;% 0.5
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
                %run the job
                spm_jobman('run',matlabbatch);
                clear matlabbatch
                %end 
            end 




      end
end


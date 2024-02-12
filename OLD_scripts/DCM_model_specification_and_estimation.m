%-----------------------------------------------------------------------
% Model specification and Estimation for the DCM rsfMRI analysis 
%-----------------------------------------------------------------------

clear all 
%for idibrom
% base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/Idibrom_rsfMRI/';
% rsfMRI_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/';

%for DSST 
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/DSST_rsfMRI_Step1/';
DCM_dir_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_DSST/DCM_DSST_Step1/';

subjects = [2:9];

%session = {'FunImg', 'S2_FunImg'};
session = {'S2_FunImg'}; 

fs = filesep;


for s = 1:length(session)
    sess = session(s); 
    for i = 1:length(subjects) 

        sub=subjects(i); 

    % %     if sub <= 9
    % %         fmdir = strcat(base,'FunImg/RglobalCFWS/Sub_00',int2str(sub));
    % %     %MPRAGEdir = strcat(rsfMRI_dir, fs, 'T1Img/Sub_00', int2str(sub));
    % %     else
    % %         fmdir = strcat(base,'FunImg/RglobalCFWS/Sub_0',int2str(sub),'/');
    % %     %MPRAGEdir = strcat(base, fs, 'T1Img/Sub_0', int2str(sub),'/');
    % %     end    
    % %     f1 = spm_select('FPList', fmdir,'^*\.nii');
    % %     Bromo_files = cellstr(f1);  
        sess = char(sess); 

        if sub <= 9
            fmdir2 = strcat([base,sess, ['/RWS/global/CF/Sub_00',int2str(sub)]]);
            %MPRAGEdir = strcat(DCM_dir, fs, 'T1Img/Sub_00', int2str(sub));
        else
            fmdir2 = strcat([base,sess,['/RWS/global/CF/Sub_0',int2str(sub)]]);
            %MPRAGEdir = strcat(base, fs, 'T1Img/Sub_0', int2str(sub),'/');
        end
        f2 = spm_select('FPList', fmdir2,'^*\.nii');
        Functional_files = cellstr(f2);
        %%%%%%%fmdir = strcat(base,fs,'FunImg/Sub_00',int2str(sub));

        cd(fmdir2)


        if sub <= 9
            DCM_dir = strcat([DCM_dir_base,sess, ['/Sub_00',int2str(sub)]]);
        %MPRAGEdir = strcat(rsfMRI_dir, fs, 'T1Img/Sub_00', int2str(sub));
        else
            DCM_dir = strcat([DCM_dir_base,sess,['/Sub_0',int2str(sub)]]);
        %MPRAGEdir = strcat(base, fs, 'T1Img/Sub_0', int2str(sub),'/');
        end
        cd(DCM_dir)

        %session = {'bromo', 'placebo'};
        %for k = 1:length(session) 

            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(DCM_dir);
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = Functional_files;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = Placebo_files;
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {''};
    %         matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';


            %Model estimation 
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            %run the job 
            spm_jobman('run',matlabbatch);
            clear matlabbatch DCM_dir
        %end 
    end 
end





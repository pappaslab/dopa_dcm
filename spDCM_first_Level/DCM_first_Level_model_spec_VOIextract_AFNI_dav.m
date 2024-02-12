%%-----------------------------------------------------------------------
% This script is to run a first level analysis for the spDCM
% It will first do: Model specification and Estimation for the DCM rsfMRI analysis 
% Then it will extract VOI time series that can then be used to specify DCM
% models 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; clc; close all; 

%% Directories
% % dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/preproc_fmriprep_spmfilt';
% % PEB_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs';

dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/preproc_afni_task';
PEB_base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/afni_spDCM_task';

%% Experiment information 
subjects = {'021', '022', '023', '024', ...
    '025', '026', '027', '028', '030', ...
    '031', '032', '033'};
%session = { 'placebo', 'tolcapone'};
session = {'bromo'}

%% Info about VOI
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
%Left hemisphere
centers = [-44 48 4; -38 28 44; -26 14 52; -52 20 28; -40 10 20; -24 4 54];


%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
%%

for i = 1:length(subjects) 
    subj=subjects(i);
    sub = char(subj); 
    subject = ['sub-' sub]; 
    
    for s = 1:length(session)
        sess = session(s); 
        cond=char(sess); 


        outputdir=([PEB_base, '/', subject, '/session-', cond, '/']); 

        
%         nifti_file = fullfile(fmriprep_preproc_data_dir, ['/sub-00', int2str(sub), '/ses-', cond]); 
%         
%         char(gunzip([nifti_file, '/*.gz']));
        
        
        SPMfile = [outputdir 'SPM.mat']; 
        %check if directory and SPM.mat file already exists; if it does then skip the model specification and VOI extraction, if not then make directory and do model specification and VOI extraction 
        if exist(SPMfile,'file')==2
            fprintf('Skipping first level analysis for subject %d and session %s\n\n',sub, cond);
        else
      
            fprintf(1,'Initialising outputdir for  subject %d and session %s and now lets run first level analysis \n', sub, cond);
            mkdir(outputdir)


            %load the fmriprep preproc data per subject and session
            preproc_dir = strcat([dcmbase,'/', subject, '/ses-', cond]);

            %load the data 
            %preproc_data = spm_select('expand', fullfile([preproc_dir, '/smooth_', subject, '_ses-', cond, '_task-rest_run-001_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc_nuisreg_4D.nii']));
            %preproc_data = spm_select('expand', fullfile([preproc_dir, '/*.nii']));
            preproc_data_files = cellstr(spm_select('FPList', preproc_dir, '/*.nii'));
                
            % put preproc_data into cellstr so that it can be read by
            % spm.stats 
            %preproc_data_files = cellstr(preproc_data);

            %set-up matlab batch 
            cd(outputdir); 
            matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(outputdir);
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = preproc_data_files;
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = -Inf;
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

% %             %run the job 
            spm_jobman('run',matlabbatch);
            clear matlabbatch 

            %extract VOI time series 
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
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0;% 0.5
                matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
                %run the job
                spm_jobman('run',matlabbatch);
                clear matlabbatch
              
            end 
         end 

        
    end 
end

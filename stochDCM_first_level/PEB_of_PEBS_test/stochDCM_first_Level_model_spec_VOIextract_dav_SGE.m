%%-----------------------------------------------------------------------
% This script is to run a first level analysis for the spDCM
% It will first do: Model specification and Estimation for the DCM rsfMRI analysis 
% Then it will extract VOI time series that can then be used to specify DCM
% models 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------
function stochDCM_first_Level_model_spec_VOIextract_dav_SGE(subject)
%clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/stochDCM';
fmriprep_preproc_data_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/preproc_fmriprep_afni';

%add SPM12 path 
SPM_path = '/home/despoC/dvogel/Toolbox/spm12/'; 
addpath(genpath(SPM_path)); 


%% Experiment information 
%subjects = {'009'}; 
SUBJECTS = {subject}; 
numSub = length(SUBJECTS);

session = {'bromo', 'placebo', 'tolcapone'};


%% Info about VOI
ROIs = {'FPI', 'MFG', 'cMFG'};
%Coordinates
centers = [-44 48 4; -38 28 44; -26 14 52];


%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
%%

for i = 1:numSub%length(subjects) 
    sub = SUBJECTS{i}; 

    %sub=subjects(i);
    %sub = char(sub);
    for s = 1:length(session)
        sess = session(s); 
        cond=char(sess); 

        %specify output directory 
        outputdir=([dcmbase, '/sub-', sub, '/session-', cond, '/']); 

        
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
            preproc_dir = strcat([fmriprep_preproc_data_dir,'/sub-', sub, '/ses-', cond]);

            
            %check first whether a .nii file exists otherwise you need to
            %gunzip it first 
            %preproc_file = ([preproc_dir, '/Fourier_Filt_func_data.nii']);
            preproc_file = ([preproc_dir, '/fmriprep_afni-sub-', sub, '-', cond, '.nii']);
            if exist(preproc_file) == 2
                %preproc_data = spm_select('expand', fullfile(preproc_dir,'/Fourier_Filt_func_data.nii'));
                preproc_data = spm_select('expand', fullfile([preproc_dir, '/fmriprep_afni-sub-', sub, '-', cond, '.nii']));
            else 
                %preproc_file_gzip = ([preproc_dir, '/Fourier_Filt_func_data.nii.gz']);
                preproc_file_gzip = ([preproc_dir, '/fmriprep_afni-sub-00', sub, '-', cond, '.nii.gz']);
                gunzip(preproc_file_gzip);
                %preproc_data = spm_select('expand', fullfile(preproc_dir,'/Fourier_Filt_func_data.nii'));
                preproc_data = spm_select('expand', fullfile([preproc_dir, '/fmriprep_afni-sub-', sub, '-', cond, '.nii']));
            end
                
            % put preproc_data into cellstr so that it can be read by
            % spm.stats 
            preproc_data_files = cellstr(preproc_data);

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

            %run the job 
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

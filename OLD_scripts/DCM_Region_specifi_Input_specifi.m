% SPM12 DCM batch script

% Settings

% -------------------------------------------------------------------------
rsfMRI_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/';

% Template DCMs go here. These could, for example, be created with the GUI.

dcm_files = {'DCM_DCM_m1.mat';
             'DCM_DCM_m2.mat';
             'DCM_DCM_m3.mat'};

% Subject names go here
subjects = {'Sub_001'};

% Scanning session
session = 1;

% Columns numbers in the design matrix to use as task inputs
% E.g. the example below is the 1st and 3rd columns with no parametric
% input. See spm_dcm_U.m.
inputs = {1 0 1};

% -------------------------------------------------------------------------
for s = 1:length(subjects)
    subject = subjects{s};   

    % We'll assume the SPM.mat for this subject is in GLM/<subject>
    % Change this appropriately
    spm_file  = fullfile(rsfMRI_dir,subject,'SPM.mat');


    % We'll assume all VOIs for this subject are in GLM/<subject
    % Change this appropriately
    voi_files = spm_select('FPList',fullfile(rsfMRI_dir,subject),'VOI_.*mat$');

   
    % Create a subdirectory for this subject
    mkdir(subject);

   cd(rsfMRI_dir)
    % Give this subject a duplicate of each template DCM
    subject_dcms = cell(size(dcm_files));   
    for f = 1:length(dcm_files)
        copyfile(subject,dcm_files{f})
        subject_dcms{s} = fullfile(subject,dcm_files{f});
    end  


    % Update this subject's DCMs to use this subject's ROIs
    clear matlabbatch;
    matlabbatch{1}.spm.dcm.fmri.regions.dcmmat = subject_dcms;
    matlabbatch{1}.spm.dcm.fmri.regions.voimat = voi_files;

   
    % Update this subject's DCMs to use this subject's experimental timing
    matlabbatch{2}.spm.dcm.fmri.inputs.dcmmat  = subject_dcms;
    matlabbatch{2}.spm.dcm.fmri.inputs.spmmat  = cellstr(spm_file);
    matlabbatch{2}.spm.dcm.fmri.inputs.session = session;
    matlabbatch{2}.spm.dcm.fmri.inputs.val     = inputs;


    % Run the above
    spm_jobman('run',matlabbatch);   

end
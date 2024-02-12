% This batch script analyses the resting state fMRI dataset
% available from the SPM website using DCM:
%   http://www.fil.ion.ucl.ac.uk/spm/data/spDCM/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf#Chap:DCM_rsfmri
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Adeel Razi & Guillaume Flandin
% $Id$

% Directory containing the data
%--------------------------------------------------------------------------
data_path = fileparts(mfilename('fullpath'));
if isempty(data_path), data_path = pwd; end
fprintf('%-40s:', 'Downloading Resting State dataset...');
urlwrite('http://www.fil.ion.ucl.ac.uk/spm/data/spDCM/spDCM.zip','spDCM.zip');
unzip(fullfile(data_path,'spDCM.zip'));
data_path = fullfile(data_path,'spDCM');
fprintf(' %30s\n', '...done');

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',true);

% Locate files / directories
%--------------------------------------------------------------------------
f      = spm_select('FPList', fullfile(data_path,'func'), '^sw.*\.img$');

RT = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL GLM FOR EXTRACTING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmdir = fullfile(data_path,'glm');
if ~exist(glmdir,'file'), mkdir(glmdir); end

clear matlabbatch;

% SPM specification
matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);

% SPM estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));

% ROI extraction
matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust  = NaN;
matlabbatch{3}.spm.util.voi.session = 1;
matlabbatch{3}.spm.util.voi.name    = 'CSF';
matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre     = [ 0 -40 -5];
matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius     = 6;
matlabbatch{3}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{3}.spm.util.voi.roi{2}.mask.image        = cellstr(fullfile(glmdir,'mask.nii'));
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{4} = matlabbatch{3};
matlabbatch{4}.spm.util.voi.name = 'WM';
matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [0 -24 -33]; 

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND GLM INCLUDING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glmdir = fullfile(data_path,'glm_corrected');
if ~exist(glmdir,'file'), mkdir(glmdir); end

clear matlabbatch;

% SPM specification
matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
    fullfile(data_path,'glm','rp_rest0000.txt'),...
    fullfile(data_path,'glm','VOI_CSF_1.mat'),...
    fullfile(data_path,'glm','VOI_WM_1.mat'),...
    }';

% SPM estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));

% ROI extraction
matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
matlabbatch{3}.spm.util.voi.adjust  = NaN;
matlabbatch{3}.spm.util.voi.session = 1;
matlabbatch{3}.spm.util.voi.name    = 'PCC';
matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre = [0 -52 26];
matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius = 8;
matlabbatch{3}.spm.util.voi.roi{2}.mask.image    = cellstr(fullfile(glmdir,'mask.nii'));
matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

matlabbatch{4} = matlabbatch{3};
matlabbatch{4}.spm.util.voi.name = 'mPFC';
matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [3 54 -2];

matlabbatch{5} = matlabbatch{3};
matlabbatch{5}.spm.util.voi.name = 'LIPC';
matlabbatch{5}.spm.util.voi.roi{1}.sphere.centre = [-50 -63 32];

matlabbatch{6} = matlabbatch{3};
matlabbatch{6}.spm.util.voi.name = 'RIPC';
matlabbatch{6}.spm.util.voi.roi{1}.sphere.centre = [48 -69 35];

spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY & ESTIMATE DCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear DCM;

% ROIs
cd(glmdir)
load('VOI_PCC_1.mat');
DCM.xY(1) = xY;
load('VOI_mPFC_1.mat');
DCM.xY(2) = xY; 
load('VOI_LIPC_1.mat');
DCM.xY(3) = xY; 
load('VOI_RIPC_1.mat');
DCM.xY(4) = xY; 

% Metadata
v = length(DCM.xY(1).u); % number of time points
n = length(DCM.xY);      % number of regions

DCM.v = v;
DCM.n = n;

% Timeseries
DCM.Y.dt  = RT;
DCM.Y.X0  = DCM.xY(1).X0;
DCM.Y.Q   = spm_Ce(ones(1,n)*v);
for i = 1:DCM.n
    DCM.Y.y(:,i)  = DCM.xY(i).u;
    DCM.Y.name{i} = DCM.xY(i).name;
end

% Task inputs
DCM.U.u    = zeros(v,1);
DCM.U.name = {'null'};         

% Connectivity
DCM.a  = ones(n,n);
DCM.b  = zeros(n,n,0);
DCM.c  = zeros(n,0);
DCM.d  = zeros(n,n,0);

% Timing
DCM.TE     = 0.04;
DCM.delays = repmat(RT,DCM.n,1);

% Options
DCM.options.nonlinear  = 0;
DCM.options.two_state  = 0;
DCM.options.stochastic = 0;
DCM.options.analysis   = 'CSD';
DCM.options.induced   = 1;
% DCM.options.nmax    = 8; 


str = sprintf('DCM_DMN');
DCM.name = str;
save(fullfile(glmdir,str),'DCM');

DCM = spm_dcm_fmri_csd(fullfile(glmdir,str));

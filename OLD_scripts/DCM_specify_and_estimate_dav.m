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
% data_path = fileparts(mfilename('fullpath'));
% if isempty(data_path), data_path = pwd; end
% fprintf('%-40s:', 'Downloading Resting State dataset...');
% urlwrite('http://www.fil.ion.ucl.ac.uk/spm/data/spDCM/spDCM.zip','spDCM.zip');
% unzip(fullfile(data_path,'spDCM.zip'));
% data_path = fullfile(data_path,'spDCM');
% fprintf(' %30s\n', '...done');

clear all 
base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/';

% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',true);

% Locate files / directories
%--------------------------------------------------------------------------
%f      = spm_select('FPList', fullfile(data_path,'func'), '^sw.*\.img$');

RT = 2;

subjects = [1:4];
N = length(subjects); 

fs = filesep;

ROIs = {'PCC', 'mPFC', 'LIPC', 'RIPC'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL GLM FOR EXTRACTING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% glmdir = fullfile(data_path,'glm');
% if ~exist(glmdir,'file'), mkdir(glmdir); end
% 
% clear matlabbatch;
% 
% % SPM specification
% matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
% matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
% matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
% matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);

% % SPM estimation
% matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));
% 
% % ROI extraction
% matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
% matlabbatch{3}.spm.util.voi.adjust  = NaN;
% matlabbatch{3}.spm.util.voi.session = 1;
% matlabbatch{3}.spm.util.voi.name    = 'CSF';
% matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre     = [ 0 -40 -5];
% matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius     = 6;
% matlabbatch{3}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
% matlabbatch{3}.spm.util.voi.roi{2}.mask.image        = cellstr(fullfile(glmdir,'mask.nii'));
% matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';
% 
% matlabbatch{4} = matlabbatch{3};
% matlabbatch{4}.spm.util.voi.name = 'WM';
% matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [0 -24 -33]; 
% 
% spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND GLM INCLUDING WM / CSF REGRESSORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% glmdir = fullfile(data_path,'glm_corrected');
% if ~exist(glmdir,'file'), mkdir(glmdir); end

% clear matlabbatch;
% 
% % SPM specification
% matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
% matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
% matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
% matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
% matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {
%     fullfile(data_path,'glm','rp_rest0000.txt'),...
%     fullfile(data_path,'glm','VOI_CSF_1.mat'),...
%     fullfile(data_path,'glm','VOI_WM_1.mat'),...
%     }';
% 
% % SPM estimation
% matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glmdir,'SPM.mat'));
% 
% % ROI extraction
% matlabbatch{3}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
% matlabbatch{3}.spm.util.voi.adjust  = NaN;
% matlabbatch{3}.spm.util.voi.session = 1;
% matlabbatch{3}.spm.util.voi.name    = 'PCC';
% matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre = [0 -52 26];
% matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius = 8;
% matlabbatch{3}.spm.util.voi.roi{2}.mask.image    = cellstr(fullfile(glmdir,'mask.nii'));
% matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';
% 
% matlabbatch{4} = matlabbatch{3};
% matlabbatch{4}.spm.util.voi.name = 'mPFC';
% matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [3 54 -2];
% 
% matlabbatch{5} = matlabbatch{3};
% matlabbatch{5}.spm.util.voi.name = 'LIPC';
% matlabbatch{5}.spm.util.voi.roi{1}.sphere.centre = [-50 -63 32];
% 
% matlabbatch{6} = matlabbatch{3};
% matlabbatch{6}.spm.util.voi.name = 'RIPC';
% matlabbatch{6}.spm.util.voi.roi{1}.sphere.centre = [48 -69 35];
% 
% spm_jobman('run',matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFY & ESTIMATE DCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i = 1:length(subjects)
        
        
        sub=subjects(i); 
        clear DCM;



        % ROIs
        cd(base)
        load(['Sub_00', int2str(sub),'/VOI_PCC_1.mat']);
        DCM.xY(1) = xY;
        load(['Sub_00',int2str(sub), '/VOI_mPFC_1.mat']);
        DCM.xY(2) = xY; 
        load(['Sub_00',int2str(sub),'/VOI_LIPC_1.mat']);
        DCM.xY(3) = xY; 
        load(['Sub_00', int2str(sub),'/VOI_RIPC_1.mat']);
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
        %DCM.a  = ones(n,n);
        %DCM.a = [1 1 1 1; 1 1 0 0; 0 0 1 1; 1 1 1 1]; 
        DCM.a = [1 1 1 1; 1 0 0 1; 0 1 1 0; 1 0 1 0];
        %DCM.a = [1 1 1 0; 1 1 0 0; 1 0 0 0; 0 0 0 0];
        DCM.b  = zeros(n,n,0);
        DCM.c  = zeros(n,0);
        DCM.d  = zeros(n,n,0);

        % Timing
        DCM.TE     = 0.028;
        DCM.delays = repmat(RT,DCM.n,1);

        % Options
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 0;
        DCM.options.analysis   = 'CSD';
        DCM.options.induced   = 1;
        % DCM.options.nmax    = 8; 


        str = sprintf(['Sub_00',int2str(sub),'DCM_model3']);
        DCM.name = str;
        save(fullfile(base,str),'DCM');

        %DCM = spm_dcm_fmri_csd(fullfile(base,str));% estimates parameters
        %of a DCM using cross spectral fMRI densities 
 end 
 
% Collate DCMs into a GCM file
% GCM = {'Sub_001DCM_DMN.mat','Sub_001DCM_DMN2.mat', 'Sub_001DCM_DMN3.mat';
%        'Sub_002DCM_DMN.mat','Sub_002DCM_DMN2.mat', 'Sub_002DCM_DMN3.mat'};
GCM = {'Sub_001DCM_model1.mat','Sub_001DCM_model2.mat';
       'Sub_002DCM_model1.mat','Sub_002DCM_model2.mat'};   

% Fully estimate model 1
GCM(:,1) = spm_dcm_fit(GCM(:,1)); 

% Use Bayesian Model Reduction to rapidly estimated models 2-N
GCM = spm_dcm_bmr(GCM);

% Alternatively, replace the above lines with this to alternate between estimating
% DCMs and estimating group effects. This is slower, but can draw subjects out 
% of local optima towards the group mean.
GCM = spm_dcm_peb_fit(GCM);

save('GCM_example.mat','GCM');

%% Esimate a second leve PEB model 
% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';

% Specify design matrix for N subjects. It should start with a constant column
M.X = ones(N,1);

% Choose field
field = {'A'};

% Estimate model
PEB = spm_dcm_peb(GCM,M,field);

save('PEB_example.mat','PEB');
%% Compare the full PEB model to nested PEB models to test specific hypotheses
% Compare nested PEB models. Decide which connections to switch off based on the 
% structure of each DCM for subject 1.
%BMA = spm_dcm_peb_bmc(PEB(1), GCM(1,:)); 
BMA = spm_dcm_peb_bmc(PEB(1));
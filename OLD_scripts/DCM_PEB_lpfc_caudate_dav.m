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
%base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_Idibrom/FunImg/';

base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_SPM_DCM/DCM_TBP/FunImg/';

DCM_models = [base, 'DCM_models_LPFC_caudate']; 
% Initialise SPM
%--------------------------------------------------------------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');
%spm_get_defaults('cmdline',true);

% Locate files / directories
%--------------------------------------------------------------------------
%f      = spm_select('FPList', fullfile(data_path,'func'), '^sw.*\.img$');

RT = 1.6;

subjects = [1:25];
N = length(subjects); 
 

fs = filesep;

%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS'}; 
%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'}; 
%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'mid_caudate'};
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'MCA', 'CCA'}; 

%DCM MODELS

%DCMdav.models{1} = [1 1 1 0; 1 1 1 0; 1 1 1 1; 0 0 1 1];
%DCMdav.models{2} = [1 0 0 0; 1 1 0 0; 1 1 1 0; 0 0 0 1];
%DCMdav.models{1} = [1 1 1 0 1 0; 1 1 1 0 1 0; 1 1 1 1 1 0; 0 0 1 1 0 1; 1 1 1 0 1 1; 0 0 0 1 1 1];
%DCMdav.models{2} = [1 0 0 0 0 0; 1 1 0 0 1 0; 1 1 1 1 0 0; 0 0 0 1 0 0; 1 1 0 0 1 1; 0 0 0 0 0 1];
% % DCMdav.models{3} = [1 0 0 0 0 0; 1 1 0 0 1 0; 1 1 1 1 0 0; 0 0 0 1 0 0; 0 1 0 0 1 0; 0 0 0 0 0 1];
% % DCMdav.models{4} = [1 0 0 0 0 0; 1 1 0 0 1 0; 1 1 1 1 0 0; 0 0 0 1 0 0; 0 1 0 0 0 1; 0 0 0 0 0 0];
% % 
%DCMdav.models{1} = [1 1 1 0 1 0 0 0; 1 1 1 0 1 0 0 0; 1 1 1 1 1 0 0 0; 0 0 1 1 0 1 0 0; 1 1 1 0 1 1 0 0; 0 0 0 1 1 1 0 0; 0 1 1 0 0 0 1 1; 0 1 1 0 0 0 1 1];
%DCMdav.models{2} = [1 0 0 0 0 0 0 0; 1 1 0 0 1 0 0 0; 1 1 1 1 0 0 0 0; 0 0 0 1 0 0 0 0; 1 1 0 0 1 1 0 0; 0 0 0 0 0 1 0 0; 0 1 1 0 0 0 1 0; 0 1 1 0 0 0 0 1];%

DCMdav.models{1} = [1 1 1 0 1 0 1 1; 1 1 1 0 1 0 1 1; 1 1 1 1 1 0 1 1; 0 0 1 1 0 1 1 1; 1 1 1 0 1 1 1 1; 0 0 0 1 1 1 1 1; 1 1 1 1 1 1 1 1; 1 1 1 1 1 1 1 1];
DCMdav.models{2} = [1 0 0 0 0 0 0 0; 1 1 0 0 1 0 0 0; 1 1 1 1 0 0 0 0; 0 0 0 1 0 0 0 0; 1 1 0 0 1 1 0 0; 0 0 0 0 0 1 0 0; 0 1 1 0 0 0 1 0; 0 1 1 0 0 0 0 1];
% % DCMdav.models{3} = [1 0 0 0 0 0 0 0; 1 1 0 0 1 0 0 0; 1 1 1 1 0 0 0 0; 0 0 0 1 0 0 0 0; 0 1 0 0 1 0 0 0; 0 0 0 0 0 1 0 0; 0 0 1 0 0 0 0 0; 0 0 1 0 0 0 0 0];
% % DCMdav.models{4} = [1 0 0 0 0 0 0 0; 1 1 0 0 1 0 0 0; 1 1 1 1 0 0 0 0; 0 0 0 1 0 0 0 0; 0 1 0 0 0 1 0 0; 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 1 0 0 0 0 0 0];

% % DCMdav.models{1} = [1 1 1 0 1 0 0; 1 1 1 0 1 0 0; 1 1 1 1 1 0 0; 0 0 1 1 0 1 0; 1 1 1 0 1 1 0; 0 0 0 1 1 1 0; 1 1 1 1 1 1 1];
% % DCMdav.models{2} = [1 0 0 0 0 0 0; 1 1 0 0 1 0 0; 1 1 1 1 0 0 0; 0 0 0 1 0 0 0; 1 1 0 0 1 1 0; 0 0 0 0 0 1 0; 1 1 1 1 1 1 1];
% % DCMdav.models{3} = [1 0 0 0 0 0 0; 1 1 0 0 1 0 0; 1 1 1 1 0 0 0; 0 0 0 1 0 0 0; 0 1 0 0 1 0 0; 0 0 0 0 0 1 0; 1 0 1 0 1 1 0];
% % DCMdav.models{4} = [1 0 0 0 0 0 0; 1 1 0 0 1 0 0; 1 1 1 1 0 0 0; 0 0 0 1 0 0 0; 0 1 0 0 0 1 0; 0 0 0 0 0 0 0; 0 1 1 1 1 0 0];



 for i = 1:length(subjects)
     sub=subjects(i);
        for m = 1:length(DCMdav.models);
             
            clear DCM;
            clear model;
            model = DCMdav.models{m};
            DCM.a = model; 
        
       
            % ROIs
            if sub <= 9
                cd(base)
                load(['Sub_00', int2str(sub),'/VOI_FPI_1.mat']);
                DCM.xY(1) = xY;
                load(['Sub_00',int2str(sub), '/VOI_MFG_1.mat']);
                DCM.xY(2) = xY; 
                load(['Sub_00',int2str(sub),'/VOI_cMFG_1.mat']);
                DCM.xY(3) = xY; 
                load(['Sub_00', int2str(sub),'/VOI_SFS_1.mat']);
                DCM.xY(4) = xY; 
                load(['Sub_00', int2str(sub),'/VOI_IFS_1.mat']);
                DCM.xY(5) = xY;
                load(['Sub_00', int2str(sub),'/VOI_IFJ_1.mat']);
                DCM.xY(6) = xY;
                load(['Sub_00', int2str(sub),'/VOI_MCA_1.mat']);
                DCM.xY(7) = xY;
                load(['Sub_00', int2str(sub),'/VOI_CCA_1.mat']);
                DCM.xY(8) = xY;
            else
                cd(base)
                load(['Sub_0', int2str(sub),'/VOI_FPI_1.mat']);
                DCM.xY(1) = xY;
                load(['Sub_0',int2str(sub), '/VOI_MFG_1.mat']);
                DCM.xY(2) = xY; 
                load(['Sub_0',int2str(sub),'/VOI_cMFG_1.mat']);
                DCM.xY(3) = xY; 
                load(['Sub_0', int2str(sub),'/VOI_SFS_1.mat']);
                DCM.xY(4) = xY; 
                load(['Sub_0', int2str(sub),'/VOI_IFS_1.mat']);
                DCM.xY(5) = xY;
                load(['Sub_0', int2str(sub),'/VOI_IFJ_1.mat']);
                DCM.xY(6) = xY; 
                load(['Sub_0', int2str(sub),'/VOI_MCA_1.mat']);
                DCM.xY(7) = xY;
                load(['Sub_0', int2str(sub),'/VOI_CCA_1.mat']);
                DCM.xY(8) = xY;
            end


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

        
        %create DCM models 

        
        
            % Connectivity full model = ones(n,n) 
            %DCM.a  = [1 1 1 0 1 0; 1 1 1 0 1 0; 1 1 1 1 1 0; 0 0 1 1 0 1; 1 1 1 0 1 1; 0 0 0 1 1 1];%model 1
            %DCM.a = [1 1 0 0 1 0; 1 1 1 0 1 0; 1 1 1 1 1 0; 0 0 0 1 0 0; 1 1 0 0 1 1; 0 0 0 0 1 1];%model 5
            %DCM.a = [1 1 1 0 1 0; 1 1 1 0 1 0; 1 1 1 1 1 0; 0 0 0 1 0 0; 1 1 0 0 1 1; 0 0 0 0 1 1]; %model6 
            %DCM.a = [1 1 1 0 1 0; 1 1 1 0 1 0; 1 1 1 1 1 0; 0 0 1 1 0 0; 1 1 0 0 1 1; 0 0 0 0 1 1]; %model7 
            DCM.b  = zeros(n,n,0);
            DCM.c  = zeros(n,0);
            DCM.d  = zeros(n,n,0);

            % Timing. . 
            DCM.TE     = 0.027;
            DCM.delays = repmat(RT,DCM.n,1);

            % Options
            DCM.options.nonlinear  = 0;
            DCM.options.two_state  = 0;
            DCM.options.stochastic = 0;
            DCM.options.analysis   = 'CSD';
            DCM.options.induced   = 1;
            % DCM.options.nmax    = 8; 

            cd(DCM_models);
            if sub <= 9
                str = sprintf(['Sub_00',int2str(sub),'DCM_model',int2str(m)]);
            else 
                str = sprintf(['Sub_0',int2str(sub),'DCM_model',int2str(m)]);
            end 
            DCM.name = str;
            
            %DCM = spm_dcm_fmri_csd(fullfile(DCM_models,str));% estimates parameters
            %of a DCM using cross spectral fMRI densities 
            save(fullfile(DCM_models,str),'DCM');
            %DCM = spm_dcm_fmri_csd(fullfile(DCM_models,str))
            %clear model
        end 

        %DCM = spm_dcm_fmri_csd(fullfile(DCM_models,str));% estimates parameters
        %of a DCM using cross spectral fMRI densities 
        %save(fullfile(DCM_models,str),'DCM')
 end 
 
% Collate DCMs into a GCM file
GCM = {'Sub_001DCM_model1.mat','Sub_001DCM_model2.mat'; %'Sub_001DCM_model3.mat', 'Sub_001DCM_model4.mat';
       'Sub_002DCM_model1.mat','Sub_002DCM_model2.mat'; %'Sub_002DCM_model3.mat', 'Sub_002DCM_model4.mat';
       'Sub_003DCM_model1.mat','Sub_003DCM_model2.mat'; %'Sub_003DCM_model3.mat', 'Sub_003DCM_model4.mat';
       'Sub_004DCM_model1.mat','Sub_004DCM_model2.mat'; %'Sub_004DCM_model3.mat', 'Sub_004DCM_model4.mat';
       'Sub_005DCM_model1.mat','Sub_005DCM_model2.mat'; %'Sub_005DCM_model3.mat', 'Sub_005DCM_model4.mat'}
       'Sub_006DCM_model1.mat','Sub_006DCM_model2.mat'; %'Sub_006DCM_model3.mat', 'Sub_006DCM_model4.mat';
       'Sub_007DCM_model1.mat','Sub_007DCM_model2.mat'; %'Sub_007DCM_model3.mat', 'Sub_007DCM_model4.mat';
       'Sub_008DCM_model1.mat','Sub_008DCM_model2.mat'; %'Sub_008DCM_model3.mat', 'Sub_008DCM_model4.mat';
       'Sub_009DCM_model1.mat','Sub_009DCM_model2.mat'; %'Sub_009DCM_model3.mat', 'Sub_009DCM_model4.mat';
       'Sub_010DCM_model1.mat','Sub_010DCM_model2.mat'; %'Sub_010DCM_model3.mat', 'Sub_010DCM_model4.mat';
       'Sub_011DCM_model1.mat','Sub_011DCM_model2.mat'; %'Sub_011DCM_model3.mat', 'Sub_011DCM_model4.mat';
       'Sub_012DCM_model1.mat','Sub_012DCM_model2.mat'; %'Sub_012DCM_model3.mat', 'Sub_012DCM_model4.mat';
       'Sub_013DCM_model1.mat','Sub_013DCM_model2.mat'; %'Sub_013DCM_model3.mat', 'Sub_013DCM_model4.mat';
       'Sub_014DCM_model1.mat','Sub_014DCM_model2.mat'; %'Sub_014DCM_model3.mat', 'Sub_014DCM_model4.mat';
       'Sub_015DCM_model1.mat','Sub_015DCM_model2.mat'; %'Sub_015DCM_model3.mat', 'Sub_015DCM_model4.mat';
       'Sub_016DCM_model1.mat','Sub_016DCM_model2.mat'; %'Sub_016DCM_model3.mat', 'Sub_016DCM_model4.mat';
       'Sub_017DCM_model1.mat','Sub_017DCM_model2.mat'; %'Sub_017DCM_model3.mat', 'Sub_017DCM_model4.mat';
       'Sub_018DCM_model1.mat','Sub_018DCM_model2.mat'; %'Sub_018DCM_model3.mat', 'Sub_018DCM_model4.mat';
       'Sub_019DCM_model1.mat','Sub_019DCM_model2.mat'; %'Sub_019DCM_model3.mat', 'Sub_019DCM_model4.mat';
       'Sub_020DCM_model1.mat','Sub_020DCM_model2.mat'; %'Sub_020DCM_model3.mat', 'Sub_020DCM_model4.mat';
       'Sub_021DCM_model1.mat','Sub_021DCM_model2.mat'; %'Sub_021DCM_model3.mat', 'Sub_021DCM_model4.mat';
       'Sub_022DCM_model1.mat','Sub_022DCM_model2.mat'; %'Sub_022DCM_model3.mat', 'Sub_022DCM_model4.mat';
       'Sub_023DCM_model1.mat','Sub_023DCM_model2.mat'; %'Sub_023DCM_model3.mat', 'Sub_023DCM_model4.mat'}; 
       'Sub_024DCM_model1.mat','Sub_024DCM_model2.mat'; %'Sub_024DCM_model3.mat', 'Sub_024DCM_model4.mat';
       'Sub_025DCM_model1.mat','Sub_025DCM_model2.mat'} %'Sub_025DCM_model3.mat', 'Sub_025DCM_model4.mat';
%        'Sub_026DCM_model1.mat','Sub_026DCM_model2.mat'; %'Sub_026DCM_model3.mat', 'Sub_026DCM_model4.mat';
%        'Sub_027DCM_model1.mat','Sub_027DCM_model2.mat'} %'Sub_027DCM_model3.mat', 'Sub_027DCM_model4.mat'; 
%        'Sub_028DCM_model1.mat','Sub_028DCM_model2.mat'; %'Sub_028DCM_model3.mat', 'Sub_028DCM_model4.mat';
%        'Sub_029DCM_model1.mat','Sub_029DCM_model2.mat'; %'Sub_029DCM_model3.mat', 'Sub_029DCM_model4.mat';
%        'Sub_030DCM_model1.mat','Sub_030DCM_model2.mat'; %'Sub_030DCM_model3.mat', 'Sub_030DCM_model4.mat';
%        'Sub_031DCM_model1.mat','Sub_031DCM_model2.mat'; %'Sub_031DCM_model3.mat', 'Sub_031DCM_model4.mat';
%        'Sub_032DCM_model1.mat','Sub_032DCM_model2.mat'; %'Sub_032DCM_model3.mat', 'Sub_032DCM_model4.mat';
%        'Sub_033DCM_model1.mat','Sub_033DCM_model2.mat'; %'Sub_033DCM_model3.mat', 'Sub_033DCM_model4.mat';
%        'Sub_034DCM_model1.mat','Sub_034DCM_model2.mat'; %'Sub_034DCM_model3.mat', 'Sub_034DCM_model4.mat';
%        'Sub_035DCM_model1.mat','Sub_035DCM_model2.mat'; %'Sub_035DCM_model3.mat', 'Sub_035DCM_model4.mat';
%        'Sub_036DCM_model1.mat','Sub_036DCM_model2.mat'; %'Sub_036DCM_model3.mat', 'Sub_036DCM_model4.mat';
%        'Sub_037DCM_model1.mat','Sub_037DCM_model2.mat'; %'Sub_037DCM_model3.mat', 'Sub_037DCM_model4.mat';
%        'Sub_038DCM_model1.mat','Sub_038DCM_model2.mat'; %'Sub_038DCM_model3.mat', 'Sub_038DCM_model4.mat';
%        'Sub_039DCM_model1.mat','Sub_039DCM_model2.mat'} %'Sub_039DCM_model3.mat', 'Sub_039DCM_model4.mat'};
% GCM = {'Sub_001DCM_model1.mat','Sub_001DCM_model4.mat';
%        'Sub_002DCM_model1.mat','Sub_002DCM_model4.mat'};   


% Fully estimate model 1
GCM(:,1) = spm_dcm_fit(GCM(:,1)); 

% Use Bayesian Model Reduction to rapidly estimated models 2-N
GCM = spm_dcm_bmr(GCM);
%[RCM,BMC,BMA] = spm_dcm_bmr(DCM); 

% Alternatively, replace the above lines with this to alternate between estimating
% DCMs and estimating group effects. This is slower, but can draw subjects out 
% of local optima towards the group mean.
%GCM = spm_dcm_peb_fit(GCM);

save('GCM_example.mat','GCM');

%% Esimate a second level PEB model 
% Specify PEB model settings (see batch editor for help on each setting)
M = struct();
M.alpha = 1; %our beliefs about the connections taking on a non-zero value are the same at the group level as at the single subject level
M.beta  = 16; % we expect the between-subject variability for each parameter to be 1/16 of our uncertainty about each parameter ath the single-subject level. 
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
BMA = spm_dcm_peb_bmc(PEB(1), GCM(1,:));
save('BMA_example.mat','BMA');%added by dav in Jan 2018
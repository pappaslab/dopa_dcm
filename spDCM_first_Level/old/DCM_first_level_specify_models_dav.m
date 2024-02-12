%%-----------------------------------------------------------------------
% This script is to run a first level analysis for the spDCM
% to specify all models 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------

clear all; clc; close all; 

%% Directories
dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/DCM_striatum_fmriprep_afni';
%fmriprep_preproc_data_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/data/fmriprep/fmriprep';

%% Experiment information 
subjects = [13];
session = {'bromo', 'placebo', 'tolcapone'}; 
%session = {'bromo', 'placebo'};
%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'Caudate_Head', 'Caudate_body'};

numROIs = length(ROIs); 
%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
N = length(subjects); 
%%
for i = 1:length(subjects) 

    sub=subjects(i);
    for s = 1:length(session)
        sess = session(s); 
        cond=char(sess);
        
        if sub <= 9 
            sub_base = ([dcmbase, '/sub-00', int2str(sub), '/session-', cond, '/']); 
        else
            sub_base = ([dcmbase, '/sub-0', int2str(sub), '/session-', cond, '/']); 
        end 
        
        SPMmat = [sub_base 'DCM_sub-0', int2str(sub), '-session-', cond, '.mat']; 
        if exist(SPMmat,'file')==2
            fprintf('Skipping DCM model specification for subject %d and session %s\n\n',sub, cond);
        else
        
            %load the SPM.mat file 
            load(fullfile(sub_base, 'SPM.mat')); 

            %for each ROI load the VOIs 
            for r = 1:length(ROIs)
                    roi = ROIs(r); 
                    roi = char(roi); 
                    load(fullfile(sub_base, ['VOI_', roi, '_1.mat']), 'xY');
                    dcm.xY(r) = xY; 
            end
            clear xY 

            dcm.n = length(dcm.xY); % number of regions
            dcm.v = length(dcm.xY(1).u); % number of time points

            % Experimental inputs (which there are none)
            dcm.U.u = zeros(dcm.v,1);
            dcm.U.name = {'null'};

            % Time series
            dcm.Y.dt  = SPM.xY.RT;
            dcm.Y.X0  = dcm.xY(1).X0;
            dcm.Y.Q = spm_Ce(ones(1,dcm.n)*dcm.v);
            for i = 1:dcm.n
                dcm.Y.y(:,i) = dcm.xY(i).u;
                dcm.Y.name{i} = dcm.xY(i).name;
            end

            % DCM parameters and options
            dcm.delays = repmat(SPM.xY.RT/2,dcm.n,1);
            dcm.TE = TE;

            dcm.options.nonlinear = 0;
            dcm.options.two_state = 0;
            dcm.options.stochastic = 0;
            dcm.options.analysis   = 'CSD';
            dcm.options.centre = 1;
            dcm.options.induced = 1;
            dcm.options.maxnodes = numROIs; 

            dcm.b = zeros(dcm.n,dcm.n);
            dcm.c = zeros(dcm.n,1);
            dcm.d = double.empty(dcm.n,dcm.n,0); 

            % Build parent model
            main_model = zeros(dcm.n); % initialise
            %main_model = ones(dcm.n); 
            main_model(eye(dcm.n) == 1) = 1; % self cons
            main_model(2,1) = 1; main_model(3,2) = 1; main_model(2,3) = 1; main_model(3,4) = 1; main_model(3,5) = 1;
            main_model(3,1) = 1; main_model(4,2) = 1; main_model(4,3) = 1; main_model(6,4) = 1; main_model(6,5) = 1;
            main_model(5,1) = 1; main_model(5,2) = 1; main_model(5,3) = 1; main_model(1,5) = 1; main_model(4,6) = 1;
            main_model(1,2) = 1; main_model(1,3) = 1; main_model(2,4) = 1; main_model(2,5) = 1; main_model(5,6) = 1;



            %build nested models 
            main_model(eye(dcm.n) == 1) = 2; %
            % flatten
            x = dcm.n * (dcm.n - 1) + dcm.n;

            main_model = reshape(main_model,x,1);
            idx = find(main_model == 1); % index of connections to turn off
            % models
            models = [];

            for i = 1:length(idx)
                models(:,i) = main_model;
                models(idx(i),i) = 0;
            end
            models = [main_model,models];

            % reorganise into NxN matrices, remove 2 values from self cons and store in cell
            TheModels = cell(1,size(models,2));
            for j = 1:size(models,2)
                TheModels{j} = reshape(models(:,j),dcm.n,dcm.n) - eye(dcm.n);
            end

            % ------------------------------------------------------------------------------
            % Specify DCM 
            % ------------------------------------------------------------------------------
            DCM = cell(1,length(TheModels));
            for imod = 1:length(TheModels)
                DCM{imod} = dcm;
                DCM{imod}.a = TheModels{imod};
            end      

            cd(sub_base)
            % ------------------------------------------------------------------------------
            % Invert the parent model and estimate the parent model 
            % ------------------------------------------------------------------------------
            fprintf(1, '\nInverting parent model... \n');
            DCM{1} = spm_dcm_fit(DCM{1});
            DCM{1} = DCM{1}{1};

            % ------------------------------------------------------------------------------
            % Invert the nested model using Bayesian model reduction
            % Use Bayesian Model Reduction to rapidly estimated models 2-N
            % ------------------------------------------------------------------------------
            fprintf(1, '\nInverting nested models using BMR... \n');
            [RCM,BMC,BMA] = spm_dcm_bmr(DCM);

            % ------------------------------------------------------------------------------
            % Save first-level outputs
            % ------------------------------------------------------------------------------
            if sub <= 9 
                save(['DCM_sub-00', int2str(sub), '-session-', cond, '.mat'],'DCM','RCM','BMC','BMA'); 
                save(['GCM_sub-00', int2str(sub), '-session-', cond, '.mat'],'DCM'); 
            else
                save(['DCM_sub-0', int2str(sub), '-session-', cond, '.mat'],'DCM','RCM','BMC','BMA'); 
                save(['GCM_sub-0', int2str(sub), '-session-', cond, '.mat'],'DCM'); 
            end 

            % save('DCM.mat','DCM')

            fprintf(1, '\t\t Finished first level DCM models \n');

            clear DCM dcm BMA BMC RCM 
        end 
       
    end
end

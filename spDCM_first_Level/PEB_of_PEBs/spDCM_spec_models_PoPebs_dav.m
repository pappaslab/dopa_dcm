%%-----------------------------------------------------------------------
% This script is to run a first level analysis for the spDCM
% to specify all models 
% D. A. Vogelsang summer 2018 
%%-----------------------------------------------------------------------


clear all; clc; close all; 

%% Directories
%dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs';
%dcmbase = '/home/despo/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs_right_HEMISPERE';
%dcmbase = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/TBP_rsfMRI/spDCM/spDCM_PEB_of_PEBs';
dcmbase = '/home/despoC/dvogel/DCM_hierarchy/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs_Schaeffer_ROIs';
%fmriprep_preproc_data_dir = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/data/fmriprep/fmriprep';

%% Experiment information 
subjects = { '002', '003', '004', '005',...
    '009', '011', '012', '013', '014',...
    '017', '018', '019', '020', '022',...
    '023', '024', '025', '026', '027',...
    '029', '030', '031', '032', '033',...
    '035', '036', '037', '038', '039',...
    '040', '041', '042', '043', '045',...
    '046', '047', '048', '049', '050',...
    '051', '052', '053', '054', '055',...
    '056', '057', '058', '060', '061',...
    '063', '064', '065', '066', '067',...
    '068', '069', '071', '072', '073',...
    '074', '076', '077', '078', '079',...
    '080'}; 

%subjects = {'001'}; 
numSub = length(subjects); 
% 
% % subFndr = dir([dcmbase '/sub*']);
% % numSub = length(subFndr);

 
sessions = {'bromo', 'placebo', 'tolcapone'}; 


%add SPM12 path 
SPM_path = '/home/despoC/dvogel/Toolbox/spm12_v7219/'; 
%SPM_path = '/home/despoC/dvogel/Toolbox/spm12/';
addpath(genpath(SPM_path)); 

%ROIs
ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ'};
%ROIs = {'FPI', 'MFG', 'cMFG', 'SFS', 'IFS', 'IFJ', 'Caudate_Head', 'Caudate_body'};

numROIs = length(ROIs); 
%% Info about dataset DA_team
TR = 2;
TE = 24; %TE is in ms 
num_slices = 36; 
% TR = 1.6
% TE = 27
% num_slices = 28
%N = length(subjects); 

%% What type of DCM?
stochastic = 0; %if stochastic = 1 then stochDCM else it is spectral DCM 

%% DCM models 
%the a matrices
DCMa.models{1} = [1 1 1 0 1 0;1 1 1 0 1 0;1 1 1 1 0 0;0 0 1 1 0 0;1 1 0 0 1 1;0 0 0 0 1 1]; %cost 18 model
%DCMa.models{1} = [1 0 1 0 0 0;0 1 1 0 0 0;0 0 1 0 0 0;0 0 1 1 0 0;0 0 0 0 1 0;0 0 0 0 1 1];%winning connections model
% DCMa.models{2} = [1 0 0; 1 1 0; 0 1 1];
% DCMa.models{3} = [1 1 0; 1 1 1; 0 1 1];
% DCMa.models{4} = [1 0 0; 1 1 0; 1 1 1];
% DCMa.models{5} = [1 1 0; 0 1 1; 0 0 1];
%[1 1; 1 1]; 
%[1 0; 1 1];
%[1 1; 1 0];


witin_models = {'bromo', 'tolca', 'null'}; 

for subIdx = 1:numSub
    %subject = subFndr(subIdx).name;

    sub=subjects(subIdx);
    subj = char(sub); 
    subject = ['sub-', subj]; 
    
    for j = 1:length(sessions)
        sess = sessions(j); 
        ses = char(sess); 
        session = ['session-' ses]; 

        

        sub_base = ([dcmbase, '/', subject, '/', session, '/']); 

        DCMmat = [sub_base 'DCM_July-2022.mat']; 
        if exist(DCMmat,'file')==2
            fprintf('Skipping DCM model specification for subject %d and %s\n\n',subject, session);
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
            
            % Experimental inputs (which there are none because it is rsfMRI)
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
            if stochastic 
                dcm.options.stochastic = 1;
                %dcm.options.analysis = '';
                dcm.options.nograph = 1; 
                dcm.options.endogenous = 1; 
            else 
                dcm.options.stochastic = 0;
                dcm.options.analysis = 'CSD';
                dcm.options.induced = 1;
                dcm.options.order = 4; %default should be 8 
            end 
            dcm.options.centre = 1;
            dcm.options.maxnodes = numROIs; 

            dcm.b = zeros(dcm.n,dcm.n);

            dcm.c = zeros(dcm.n,1);
            dcm.d = double.empty(dcm.n,dcm.n,0); 
            %DCM = dcm; 

            % ------------------------------------------------------------------------------
            % Specify DCM 
            % ------------------------------------------------------------------------------
            DCM = cell(1,length(DCMa.models));        
            for imod = 1:length(DCMa.models)

                DCM{imod} = dcm;
                DCM{imod}.a = DCMa.models{imod};
                
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
            output_base = [dcmbase, '/DCM_models/']; 
            cd(output_base); 
            %save(['DCM_', subject, '-', session, '-models.mat'],'DCM','RCM','BMC','BMA'); 
            save(['DCM_', subject, '-', session, '-models.mat'],'DCM', 'RCM', 'BMC', 'BMA'); 


            fprintf(1, '\t\t Finished first level DCM models \n');

            clear DCM dcm BMA BMC RCM 
        end 
    end 
       
    %end
end

function SPM = collapseSessions(SPM,scans)
% Takes in an SPM design
% collapses multiple sessions together
% scans is a vector that specifies the scan lengths for the collapsed model
% For example, if SPM.nscan = [200 200 200 200]
% scans = [400 400] to collapse sessions 1 and 2 together and 3 and 4
% together
% scans = [800] would collapse all sessions into a single session
%
% Recommended use: first create and estimate a design matrix with each run 
% as a session as you would normally
% Then supply that SPM.mat along with scans
%
% Scaling and masking will be copied over to new design
% Run regressors will be created for each run of the original design
% Temporal non-sphericity correction and high-pass filter will be applied 
% for each run of the original design
% Edge artifact regressor will be created to deal with edge effects from
% high-pass filter
% For best results, runs should start and end with un-modeled fixation



NEE.xY.RT = SPM.xY.RT;
NEE.xY.P = SPM.xY.P;
NEE.xY.VY = SPM.xY.VY;
NEE.xBF = SPM.xBF;
NEE.nscan = scans;

%% determine collapsing scheme

numSess = length(scans);
oldNumSess = length(SPM.nscan);
oldNscan = SPM.nscan;

sessi = 1;
scani = 1;
collapse = cell(numSess,1);

while sessi <= oldNumSess && scani <=numSess
    csum = cumsum(SPM.nscan(sessi:end));
    idx = find(csum==scans(scani))+sessi-1;
    collapse{scani} = sessi:idx;
    scani = scani+1;
    sessi = idx+1;
end

%need to add in some checks here to make sure everything adds up

%% rearrange SPM.Sess
for sessi = 1:numSess
    
    %first collect all of the conditions
    conds = {};
    for ci = collapse{sessi}
        conds = union(conds,[SPM.Sess(ci).U.name]);
    end
    
    for condi = 1:length(conds)
        NEE.Sess(sessi).U(condi).name = conds(condi);
        NEE.Sess(sessi).U(condi).ons = [];
        NEE.Sess(sessi).U(condi).dur = [];
        NEE.Sess(sessi).U(condi).P.name = 'none'; %temporary, will need to deal with parametric modualtors later
        for ci = collapse{sessi}
            idx = find(strcmp(conds{condi},[SPM.Sess(ci).U.name]));
            if ~isempty(idx)
                if ~isscalar(idx)
                    disp('error: found multiple conditions of the same name');
                end
                onsetOffset = sum(SPM.nscan(collapse{sessi}(1)+1:ci))*(strcmp(SPM.xBF.UNITS,'secs')*SPM.xY.RT);

                NEE.Sess(sessi).U(condi).ons = [NEE.Sess(sessi).U(condi).ons; (SPM.Sess(ci).U(idx).ons+onsetOffset)];
                NEE.Sess(sessi).U(condi).dur = [NEE.Sess(sessi).U(condi).dur; SPM.Sess(ci).U(idx).dur];
            end
        end
    end
    
    %now collect other regressors
    conds = {};
    for ci = collapse{sessi}
        conds = union(conds,[SPM.Sess(ci).C.name]);
    end
    
    if isempty(conds)
        %add regressor to capture filter-related edge effects
        NEE.Sess(sessi).C.name = {'EdgeReg'};
        C = zeros(NEE.nscan(sessi),1);
        tmp1 = cumsum([0 oldNscan(collapse{sessi})]) + 1;
        tmp2 = cumsum(oldNscan(collapse{sessi}));
        t = sort([tmp1(1:end-1) tmp2]);
        C(t) = 1;
        NEE.Sess(sessi).C.C = C;
    else
        NEE.Sess(sessi).C.C = zeros(NEE.nscan(sessi),length(conds));
        for condi = 1:length(conds)
            NEE.Sess(sessi).C.name{condi} = conds{condi};
            csum = cumsum([0 oldNscan(collapse{sessi})]);
            for si = 1:length(collapse{sessi})
                ci = collapse{sessi}(si);
                idx = find(strcmp(conds{condi},[SPM.Sess(ci).C.name]));
                scanidx = (1:oldNscan(ci)) + csum(si);
                if ~isempty(idx)
                    if ~isscalar(idx)
                        disp('error: found multiple conditions of the same name');
                    end
                    NEE.Sess(sessi).C.C(scanidx,condi) = SPM.Sess(ci).C.C(:,idx);
                end
            end
        end
        %add regressor to capture filter-related edge effects
        NEE.Sess(sessi).C.name{end+1} = 'EdgeReg';
        C = zeros(NEE.nscan(sessi),1);
        tmp1 = cumsum([0 oldNscan(collapse{sessi})]) + 1;
        tmp2 = cumsum(oldNscan(collapse{sessi}));
        t = sort([tmp1(1:end-1) tmp2]);
        C(t) = 1;
        NEE.Sess(sessi).C.C = [NEE.Sess(sessi).C.C C];
    end
end

NEE = spm_fMRI_design(NEE,0);

%the rest of the code is adapted from spm_fmri_concatenate

%% add block regressors for run means
Xb = [];
Bn = {};
for sessi = 1:oldNumSess
    Xb = blkdiag(Xb,ones(oldNscan(sessi),1));
    Bn{sessi} = sprintf('Sn(%i) constant',sessi);
end
blockIdx = NEE.xX.iB(1); %index for the first block regressor

NEE.xX.X = [NEE.xX.X(:,1:blockIdx-1) Xb];
NEE.xX.iB = blockIdx:(blockIdx+size(Xb,2)-1);
NEE.xX.name = {NEE.xX.name{1:blockIdx-1} Bn{:}};

%% High-pass filter
s = cumsum([0 oldNscan]);

for sessi = 1:oldNumSess
    K(sessi) = struct('HParam', 128,...
                      'row', s(sessi) + (1:oldNscan(sessi)),...
                      'RT', NEE.xY.RT);
end
NEE.xX.K = spm_filter(K);

%% Temporal non-sphericity
NEE.xVi.form = 'AR(0.2)';
NEE.xVi.Vi = spm_Ce(oldNscan,0.2);


%% Scaling and masking
NEE.xGX = SPM.xGX;
NEE.xM = SPM.xM;
NEE.xVol = SPM.xVol; 
NEE.xX = SPM.xX; 
NEE.Vbeta = SPM.Vbeta; 
NEE.swd = SPM.swd;
NEE.VM = SPM.VM; 
SPM = NEE;
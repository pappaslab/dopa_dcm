%Correlation script 

clear all; 
clc; 
close all; 

base = '/home/despoC/dvogel/rsfMRI_DA_Hierarchy/rsfMRI_DA_Analysis_DAVID/fMRIprep/DA_team_rsfMRI/spDCM/spDCM_PEB_of_PEBs/';
dcmbase = [base, 'DCM_models']; 
%output directory 
output_dir = [base, 'Second_Level/']; 
savedir = [base, 'Correl_values_PEB_Ep/']

cd(output_dir)

load([output_dir 'PEB_of_PEBS_July2019.mat']);

subject = length(PEBs); 

for subIdx = 1:subject
    tolc_vs_plac_cMFG_FPI_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(49);
    tolc_vs_plac_cMFG_MFG_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(50);
    tolc_vs_plac_cMFG_SFS_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(52);
    tolc_vs_plac_IFS_MFG_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(56);
    brom_vs_plac_cMFG_FPI_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(29);
    brom_vs_plac_cMFG_MFG_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(30);
    brom_vs_plac_cMFG_SFS_Ep{subIdx,1} = PEBs{subIdx,1}.Ep(32);
    
end 


save([savedir, 'PEB_PET_correl_output']); 

#!/bin/sh

SUB_ID=${SGE_TASK} 
echo "subject for this spDCM first level VOI extract is $SUB_ID"
cd /home/despoC/dvogel/fMRI_script_Despolab/DCM_scripts/spDCM_first_Level
matlab-2015a-spm12 -nosplash -r "DCM_first_Level_model_spec_VOIextract_AFNI_dav_SGE('$SUB_ID')"

#!/bin/sh

SUB_ID=${SGE_TASK} 
echo "subject is for this stochDCM first level is $SUB_ID"
cd /home/despoC/dvogel/fMRI_script_Despolab/DCM_scripts/stochDCM_first_level/model_spec_VOIextract
matlab-2015a-spm12 -nosplash -r "stochDCM_first_Level_model_spec_dav_SGE('$SUB_ID')"

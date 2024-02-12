#!/bin/sh

SUB_ID=${SGE_TASK} 
echo "subject is for this stochDCM first level is $SUB_ID"
cd /home/despoC/dvogel/fMRI_script_Despolab/DCM_scripts/spDCM_first_Level/PEB_of_PEBs
matlab-2015a-spm12 -nosplash -r "spDCM_spec_models_PoPebs_dav_SGE('$SUB_ID')"


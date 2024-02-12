#!/bin/sh

SUB_ID=${SGE_TASK} 
echo "subject is for this stochDCM first level is $SUB_ID"
cd /home/despoC/dvogel/fMRI_script_Despolab/DCM_scripts/stochDCM_first_level/test_U_u_matrix
matlab-2015a-spm12 -nosplash -r "stochDCM_first_level_specify_models_dav_Test_U_u_hardcode_SGE('$SUB_ID')"


#!/bin/sh

SUB_ID=$(echo ${SGE_TASK_ID} | awk '{printf "%2.2d", $1}')
echo "subject is for this stochDCM first level is $SUB_ID"
cd /home/despoC/dvogel/fMRI_script_Despolab/DCM_scripts/stochDCM_first_level/test_U_u_matrix
matlab-2015a-spm12 -nosplash -r "stochDCM_SGE('$SUB_ID')"


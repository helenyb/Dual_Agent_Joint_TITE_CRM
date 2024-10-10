# Dual_Agent_Joint_TITE_CRM
README for Dual Agent Joint TITE-CRM

# Data Generation:
`data_generation_TTE_v2.R`

This file contains the functions to generate the time-to-event data in the simulation studies.

# Prior Calibration:

## Prior Calibration Functions:
`JOINT_TITE_BLRM_priorcal.R`

`JOINT_TITE_POCRM_priorcal.R`

These two functions each source data_generation_TTE_v2.R. They contain the necessary functions to conduct simulated trials for each of the two methods using only toxicity outcomes, in order to calibrate the priors.

## Prior Calibration Execution:
`calibration_step1A_JOINT_TITE_POCRM.R`

`calibration_step1B_JOINT_TITE_POCRM.R`

These two scripts run the first step of the calibration of the toxicity prior hyperparameter calibration for the Joint TITE-POCRM, calibrating tox only.

`calibration_step2_JOINT_TITE_POCRM.R`

This script runs one iteration of the second step of the calibration of the prior hyperparameter calibration for the Joint TITE-POCRM.

`calibration_step1A _JOINT_TITE_BLRM.R`

`calibration_step1B _JOINT_TITE_BLCRM.R`

These two scripts run the first step of the calibration of the toxicity prior hyperparameter calibration for the Joint TITE-BLRM, calibrating tox only.

`calibration_step2_JOINT_TITE_BLRM.R`

This script runs one iteration of the second step of the calibration of the prior hyperparameter calibration for the Joint TITE-BLRM.

## Prior Calibration Results Extraction:

`calibration_step1A_JOINT_TITE_POCRM_results.R`

`calibration_step1B_JOINT_TITE_POCRM_results.R`

These two scripts extract the results from the first step of the calibration of the toxicity prior hyperparameter calibration for the Joint TITE-POCRM.

`calibration_step1A_JOINT_TITE_BLRM_results.R`

`calibration_step1B_JOINT_TITE_BLCRM.R_results`

These two scripts extract the results from the first step of the calibration of the toxicity prior hyperparameter calibration for the Joint TITE-BLRM.

`calibration_step2_results.R`

This script extracts the results from the second step of the calibration of the prior hyperparameter calibration for the Joint TITE-BLRM and Joint TITE-POCRM.

# Simulations:

## Simulation Functions:

`JOINT_TITE_BLRM.R`

`JOINT_TITE_POCRM.R`

These two functions each source data_generation_TTE_v2.R. They contain the necessary functions to conduct full simulated trials for each of the three methods.

## Simulation Execution:

`simulations_JOINT_TITE_BLRM.R`

`simulations_JOINT_TITE_POCRM.R`

These two scripts each source their respective functions, which in turn source data_generation_TTE_v2.R. They execute the simulations.

## Simulation Results Extraction:

`results_extraction_POCRM_BLRM.R`

This script extracts the simulation results for the Joint TITE-BLRM and Joint TITE-POCRM.


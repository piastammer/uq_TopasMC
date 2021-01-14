# uq_TopasMC
Uncertainty Quantification post-processing for output from TopasMC engine for RT dose calculations


## Example 

Call 

postProcess( './example/EventMaps/', './example/WorkspaceBoxPhantom_Pia_SAD1e6.mat', 1 , [25 60 25] , [9 9 9] , 50,0, 0)

to compute the example with 3mm standard deviation for set-up errors, no range uncertainties and 50 shifts to compute the variance estimate.

## How to run the TOPAS simulation 

- export workspace from matRad (ct, cst, machine, pln, stf, dij, resultGUI)

- Add extension DumpEventMapsDoseToMedium to TOPAS (see TOPAS README)

- Exchange file "forwardMC_TOPAS_setup_generic_protons.txt" from MatRad-Topas interface with the file given here (includes new Scorer)

- run ./forwardMC Workspace_name.mat (starts TOPAS simulation for that workspace) 
--> results and EventMaps will be saved in dataPIDxxxx_Workspace_name

- postProcessing can then be started for those EventMaps and workspace: 
postProcess( '/path/to/EventMaps/', './path/to/Workspace/Workspace_name.mat', number_of_beams , cube_dimensions , set_up error variance (in mm^2)  , number realisations,range error standard deviation (in percent), beam energy spread (in percent))

(To run TOPAS simulation (get expected value) for different uncertainties call 
./runTopasUQ.sh sigma_su_x sigma_su_y sigma_su_z sigma_range data_PIDxxxx_Workspace_name.mat beamSetup_matrad_plan_field

sigma_su_x/y/Z set-up error std. dev. 
sigma_range range error std. dev (in percent))

# uq_TopasMC
Uncertainty Quantification post-processing for output from TopasMC engine for RT dose calculations

Call 

postProcess( './example/EventMaps/', './example/WorkspaceBoxPhantom_Pia_SAD1e6.mat', 1 , [25 60 25] , [9 9 9] , 50,0, 0)

to compute the example with 3mm standard deviation for set-up errors, no range uncertainties and 50 shifts to compute the variance estimate.

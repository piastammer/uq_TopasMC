#!/bin/bash

#$1 Std Dev Beam Position: Sigmax
#$2 Std Dev Beam Position: Sigmay
#$3 Std Dev Beam Position: Sigmaz
#$4 Std Dev Range in %
#$5 Path of folder for base parameter files
#$6 Name of parameter file base (for file_1.txt,...,file_n.txt write "file_")

source ../../../setup_env.sh
cd $5
echo $PWD
#Set value of pi
pi=$(echo "h=10;4*a(1)" | bc -l)

#Go through all beams (each in seperate parameter file)
b=1
for parameterFile in $(find . -name "$6?.txt")
do
  echo "Start loop"
  ##Beam Position uncertainties##
  #Read in Couch and Gantry Angles from $parameterFile
  Phi=$(awk '/GantryAngle = /{match($0, /\ =(.*)\ deg/, a); print a[1];}' $parameterFile)
  Theta=$(awk '/CouchAngle = /{match($0, /\ =(.*)\ deg/, a); print a[1];}' $parameterFile)
  
  echo "Phi="
  echo $Phi
  echo "Theta="
  echo $Theta

  #Read in beam FWHM
  #Read whole vector including number of time steps
  FWHM_Values=$(awk '/FocusFWHM\/Values = /{match($0, /\ =(.*)\ mm/, a); print a[1];}' $parameterFile)
  
  #Cut off everything except one entry (last one)
  
  FWHM=$(echo "$FWHM_Values" | awk '{print $NF}')
  
  echo "FWHM="
  echo $FWHM

  #Transform from grad to rad
  Theta=$(echo "scale=6; $Theta*$pi/180" | bc)
  Phi=$(echo "scale=6; $Phi*$pi/180" | bc)
  
  echo "Phi="
  echo $Phi
  echo "Theta="
  echo $Theta

  #Compute error variances for rotated coordinate system
  SigmaX=$(echo "scale=4; ($1*c($Phi)+$2*s($Phi))*c($Theta) + $3*s($Theta)" | bc -l)
  SigmaZ=$(echo "scale=4; ($1*c($Phi)+$2*s($Phi))*s($Theta) - $3*c($Theta)" | bc -l)
  
  echo "SigmaX="
  echo $SigmaX
  echo "SigmaZ="
  echo $SigmaZ

  #Compute convoluted beam and error variance
  Sigma_convX=$(echo "scale=4; sqrt($SigmaX^2 + (0.424661*$FWHM)^2)" | bc -l)
  Sigma_convZ=$(echo "scale=4; sqrt($SigmaZ^2 + (0.424661*$FWHM)^2)" | bc -l)

  echo "Sigma_convX="
  echo $Sigma_convX
  echo "Sigma_convZ="
  echo $Sigma_convZ

  #Write new values into parameter file
  sed -i "s~d:So/PencilBeam/BeamPositionSpreadX = Tf/Beam/FocusFWHM/Value mm \* Sim/FWHM2SIGMA~d:So/PencilBeam/BeamPositionSpreadX = $Sigma_convX mm~;s~d:So/PencilBeam/BeamPositionSpreadY = Tf/Beam/FocusFWHM/Value mm \* Sim/FWHM2SIGMA~d:So/PencilBeam/BeamPositionSpreadY = $Sigma_convZ mm~" $parameterFile
  
  ##Beam Range uncertainties##
  #Set parameters of Bragg-Kleemann equation
  alpha=0.022
  p=1.77
 
  #Read in beam energy spread
  Sigma_energy=$(awk '/BeamEnergySpread = /{match($0, /\ =(.*)\ #/, a); print a[1];}' $parameterFile)
  echo "Sigma_energy="
  echo $Sigma_energy
  #Compute new energy spread including uncertainty
  Sigma_convEnergy=$(echo "scale=4; sqrt($Sigma_energy^2 + ((1/$p)*$4)^2 )" | bc -l)
  echo "Sigma_convEnergy="
  echo $Sigma_convEnergy
  #Write new value into parameter file
  sed -i "s~u:So/PencilBeam/BeamEnergySpread = .*#~u:So/PencilBeam/BeamEnergySpread = $Sigma_convEnergy~" $parameterFile

 cat >> $parameterFile << EOF 

i:Sc/Dose/BeamNumber = $b

EOF

  echo "Submitting job $parameterFile"	
  ../../topas/topas $parameterFile

  ((b++))
done


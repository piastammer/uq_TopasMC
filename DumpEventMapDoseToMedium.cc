// Scorer for DumpEventMapDoseToMedium
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "DumpEventMapDoseToMedium.hh"
#include "TsTrackInformation.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

DumpEventMapDoseToMedium::DumpEventMapDoseToMedium(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
								 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{   
	//this_id = std::this_thread::get_id();
    //std::stringstream fName;
	//fName << "EventMaps_" << this_id << ".txt";
	
	//fOutFile.open (fName.str(), std::ios::out | std::ios::app); 
	//fOutFile.open ("EventMaps.txt", std::ios::out | std::ios::app); 
	beam_number = fPm->GetIntegerParameter(GetFullParmName("BeamNumber"));
	this_id = std::this_thread::get_id();
    std::stringstream fName;
	fName << "EventMaps_ID_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	
	//fOutFileEvtIds.open (fName.str(), std::ios::out | std::ios::binary | std::ios::app); 
	//fOutFile.open ("EventMaps.txt", std::ios::out | std::ios::app);
	

	//fOutFileEvtIds.close();

	//std::cout << "EventIDs written to " << fName.str() << std::endl;
    
	counter = 0;

	SetUnit("Gy");
	
	pM->SetNeedsTrackingAction();
}




G4bool DumpEventMapDoseToMedium::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

    // Retrieve vertex once per event
    if (!fEvtMap->size()) {

		TsTrackInformation* parentInformation = (TsTrackInformation*)(aStep->GetTrack()->GetUserInformation());
		if (parentInformation) {

			const G4Event* currentEvent = G4RunManager::GetRunManager()->GetCurrentEvent();
            primaryKinEnergy = currentEvent->GetPrimaryVertex(0)->GetPrimary(0)->GetKineticEnergy();

			std::vector<G4ThreeVector> vertexPositions = parentInformation->GetParentTrackVertexPositions();
            if (vertexPositions.size()) {
                vertex = vertexPositions[0];
            }
            else {
                vertex = aStep->GetTrack()->GetVertexPosition();
            }

		} else {
            G4String message = "no user track information";

            G4Exception(__FILE__, "ProcessHits()",
                        FatalException,message);
		}
    }

	G4double edep = aStep->GetTotalEnergyDeposit();
	if ( edep > 0. ) {
		G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();

		ResolveSolid(aStep);

		G4double dose = edep / ( density * fSolid->GetCubicVolume() );
		dose *= aStep->GetPreStepPoint()->GetWeight();

		AccumulateHit(aStep, dose);

		return true;
	}
	return false;
}


void DumpEventMapDoseToMedium::UserHookForEndOfEvent()
{
	G4int evtID = GetEventID();

    G4double x_mm = vertex.x() / mm;
    G4double y_mm = vertex.y() / mm;
    G4double z_mm = vertex.z() / mm;
    
	// Iterator to use with fEvtMap, the map of hits in the event
	std::map<G4int, G4double>::iterator itrX;

	// Bin index
	G4int index;

	// Value from one specific bin
	G4double dose_Gy;

	this->initVertex.push_back(x_mm);
	this->initVertex.push_back(y_mm);
	this->initVertex.push_back(z_mm);	
	this->initEnergy.push_back(primaryKinEnergy/MeV );
	//fOutFileEvtIds << evtId;
	

    G4int i=0;
	for (itrX = fEvtMap->begin(); itrX != fEvtMap->end(); itrX++) {
		index = itrX->first;
		dose_Gy = itrX->second / (joule/kg);
		this->voxValues.push_back(dose_Gy);
    	this->voxIndex.push_back(index);    
		i++;
			
	}
	 this->numEntries.push_back(i);
     this->evtStartIdx.push_back(voxIndex.size()-i);

	this_id = std::this_thread::get_id();
    std::stringstream fName1;
	fName1 << "EventMaps_ID_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileStartIdx.open (fName1.str(), std::ios::out | std::ios::binary | std::ios::app); 
	fOutFileStartIdx.write((char*) &evtStartIdx[evtStartIdx.size()-1],sizeof(G4int));
	fOutFileStartIdx.close();

	std::stringstream fName2;
	fName2 << "EventMaps_numEntries_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileNumEntr.open (fName2.str(), std::ios::out | std::ios::binary | std::ios::app); 
	fOutFileNumEntr.write((char*) &numEntries[numEntries.size()-1], sizeof(G4int));
	fOutFileNumEntr.close();

	std::stringstream fName3;
	fName3 << "EventMaps_voxVal_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileVoxVal.open (fName3.str(), std::ios::out | std::ios::binary | std::ios::app); 
	fOutFileVoxVal.write((char*) &voxValues[voxValues.size()-i],i *sizeof(G4double));
	fOutFileVoxVal.close();

	std::stringstream fName4;
	fName4 << "EventMaps_voxIdx_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileVoxIdx.open (fName4.str(), std::ios::out | std::ios::binary | std::ios::app); 
	fOutFileVoxIdx.write((char*) &voxIndex[voxIndex.size()-i],i*sizeof(G4int));
	fOutFileVoxIdx.close();

	std::stringstream fName5;
	fName5 << "EventMaps_initVert_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileInitVert.open (fName5.str(), std::ios::out | std::ios::binary | std::ios::app); 
	// fOutFileInitVert.write((char*) &x_mm, sizeof(G4double));
	// fOutFileInitVert.write((char*) &y_mm, sizeof(G4double));
	// fOutFileInitVert.write((char*) &z_mm, sizeof(G4double));
	fOutFileInitVert.write((char*) (&initVertex.back()-2),3*sizeof(G4double));
	fOutFileInitVert.close();
	
	std::stringstream fName6;
	fName6 << "EventMaps_initEnergy_" << this_id << "_" << std::to_string(beam_number) <<".bin";
	fOutFileInitEnergy.open (fName6.str(), std::ios::out | std::ios::binary | std::ios::app); 
	fOutFileInitEnergy.write((char*) &initEnergy.back(),sizeof(G4double));
	fOutFileInitEnergy.close();
	
	
	
  
 // std::lock_guard<std::mutex> lock(m);

//  fOutFile << "EvtID"
//            << "\tvertex"
//            << "\tbinIx"
//            << "\tdose"
//            << G4endl;
	// Iterate over all the hits in this event
	/*
	for (itrX = fEvtMap->begin(); itrX != fEvtMap->end(); itrX++) {
		index = itrX->first;
		dose_Gy = itrX->second / (joule/kg);

     fOutFile  << evtID << "\t"
            //    << vertex/mm << "\t"
			   << x_mm << "\t"
			   << y_mm << "\t"
			   << z_mm << "\t"
               << index << "\t"
               << dose_Gy << "\t"
               << G4endl;
	}
	*/

  

}

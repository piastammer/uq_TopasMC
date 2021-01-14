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

#ifndef DumpEventMapDoseToMedium_hh
#define DumpEventMapDoseToMedium_hh
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
//#include <mutex>
#include "TsVBinnedScorer.hh"

//std::mutex m;

class DumpEventMapDoseToMedium : public TsVBinnedScorer
{
public:
	DumpEventMapDoseToMedium(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

    //static ~DumpEventMapDoseToMedium();
    ~DumpEventMapDoseToMedium() {		
	//fOutFileEvtIds.close();
    std::cout<<"Constructor";
    }
    std::ofstream fOutFileStartIdx;
    std::ofstream fOutFileVoxVal;
    std::ofstream fOutFileVoxIdx;
    std::ofstream fOutFileNumEntr;
    std::ofstream fOutFileInitVert;
    std::ofstream fOutFileInitEnergy;
    std::thread::id this_id;
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    void close();
    void UserHookForEndOfEvent();
protected:
    G4ThreeVector vertex;
	G4double primaryKinEnergy;
private:
    std::vector<G4double> voxValues;
    std::vector<G4int> voxIndex;
    std::vector<G4double> initVertex;
    std::vector<G4int> evtStartIdx;
    std::vector<G4int> numEntries;
	std::vector<G4double> initEnergy;
    G4int counter;
    G4int beam_number;


};
#endif

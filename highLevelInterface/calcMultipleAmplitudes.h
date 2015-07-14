#include<iostream>
#include<string>
#include<vector>
#include<complex>
#include<map>
#include<utility>

#include"TFile.h"
#include"TTree.h"
#include"TTreePerfStats.h"
#include"TClonesArray.h"

#include"waveDescription.h"
#include"isobarAmplitude.h"
#include"particleDataTable.h"
#include"eventMetadata.h"
#include"ampIntegralMatrix.h"
namespace rpwa {
	namespace hli {
		std::vector<std::complex<double> > evaluateAmplitudes(std::vector<rpwa::isobarAmplitudePtr> &amplitudes, TClonesArray prodKinematics, TClonesArray decayKinematics);
		bool initAmplitudesKinematics(std::vector<rpwa::isobarAmplitudePtr> &amplitudes, std::vector<std::string> prodNames, std::vector<std::string> decayNames);
		bool checkBinnings(std::vector<std::map<std::string,std::pair<double, double> > > binnings);

		bool calcBinnedIntegralsFromEventTree(  const eventMetadata*                            eventMeta, 
		                                        std::vector<isobarAmplitudePtr>                 &amplitudes, 
		                                        std::vector<ampIntegralMatrix*>                  &matrix,
		                                        std::vector<std::map<std::string,std::pair<double,double> > >binning                    = std::vector<std::map<std::string,std::pair<double,double> > >(0),
		                                        const long int                                  maxNmbEvents                            = -1, 
		                                        const long int                                  startEvent                              =  0,
		                                        bool                                            printProgress                           = true,
		                                        const std::string&                              treePerfStatOutFileName                 = "", 
		                                        const long int                                  treeCacheSize                           = 25000000);
	};
};
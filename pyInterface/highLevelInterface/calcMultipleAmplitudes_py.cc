#include "calcMultipleAmplitudes_py.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

#include<iostream>
namespace bp = boost::python;
bool calcBinnedIntegralsFromEventTree_py( rpwa::eventMetadata&       dataFile,
                                          bp::object&                amplitudesPy,
                                          bp::object&                integralsPy,
                                          const bp::object           tBinningPy,
                                          const long int             maxNmbEvents,
                                          const long int             startEvent,
                                          bool                       printProgress, 
                                          const std::string&         treePerfStatOutFileName,
                                          const long int             treeCacheSize ){

	std::vector<rpwa::isobarAmplitudePtr> amplitudesVector;
	if (not rpwa::py::convertBPObjectToVector<rpwa::isobarAmplitudePtr>(amplitudesPy, amplitudesVector)){
		PyErr_SetString(PyExc_TypeError, "could not extract amplitudes");
		bp::throw_error_already_set();
	};

	std::vector<rpwa::ampIntegralMatrix*> integralsVector;
	if (not rpwa::py::convertBPObjectToVector<rpwa::ampIntegralMatrix*>(integralsPy, integralsVector)){
		PyErr_SetString(PyExc_TypeError, "could not extract integrals");
		bp::throw_error_already_set();
	};

	std::vector<bp::dict> binningIntermediate;
	if (not rpwa::py::convertBPObjectToVector<bp::dict>(tBinningPy,binningIntermediate)){
		PyErr_SetString(PyExc_TypeError, "could not extract binning (step1)");
		bp::throw_error_already_set();
	};
	std::vector<std::map<std::string,std::pair<double,double> > > tBinningFinal;
	for (size_t nn =0; nn<binningIntermediate.size();++nn){
		std::map<std::string, std::pair<double, double> > binningMap;
		{
			bp::list keys = binningIntermediate[nn].keys();
			for(int i = 0; i < bp::len(keys); ++i) {
				std::pair<double, double> element;
				if(not rpwa::py::convertBPObjectToPair<double, double>(binningIntermediate[nn][keys[i]], element)) {
					PyErr_SetString(PyExc_TypeError, "Got invalid pair for binningMap when executing rpwa::calcBinnedIntegralsFromEventTree(...)");
					bp::throw_error_already_set();
				}
				bp::extract<std::string> getString(keys[i]);
				if(not getString.check()) {
					PyErr_SetString(PyExc_TypeError, "Got invalid key for binningMap when executing rpwa::calcBinnedIntegralsFromEventTree(...)");
					bp::throw_error_already_set();
				};
				binningMap[getString()] = element;
			};
		};
		tBinningFinal.push_back(binningMap);
	};
	return rpwa::hli::calcBinnedIntegralsFromEventTree(&dataFile, amplitudesVector, integralsVector, tBinningFinal, maxNmbEvents, startEvent, printProgress, treePerfStatOutFileName, treeCacheSize);
};


void rpwa::py::exportCalcMultipleAmplitudes(){
	bp::def(
		"calcBinnedIntegralsFromEventTree"
		, &calcBinnedIntegralsFromEventTree_py
		, (bp::arg("eventMeta"),
		   bp::arg("amplitudes"),
		   bp::arg("integralMatrices"),
		   bp::arg("binnings"),
		   bp::arg("maxNmbEvents") = -1,
		   bp::arg("startEvent") =  0,
		   bp::arg("printProgress") = false,
		   bp::arg("treePerfStatOutFileName") = "",
		   bp::arg("treeCacheSize") = 25000000)
	);

	bp::register_ptr_to_python<rpwa::isobarAmplitudePtr>();
};

import os
import resource

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def readWaveList(waveListFileName, keyFiles):
	pyRootPwa.utils.printInfo("reading amplitude names and thresholds from wave list file "
	          + "'" + waveListFileName + "'.")
	with open(waveListFileName, 'r') as waveListFile:
#	if (not waveListFile) {
#		printErr << "cannot open file '" << waveListFileName << "'. Aborting..." << endl;
#		throw;
#	}
		waveDescThres = []
		lineNmb = 0
		for line in waveListFile:
			if (line[0] == '#'):  # comments start with #
				continue
			line = line.replace('\n', '')
			lineArray = line.split(" ")
			if(len(lineArray) >= 1 and len(lineArray) <= 2):
				waveName = lineArray[0]
				if(len(lineArray) == 1):
					threshold = 0
				else:
					threshold = lineArray[1]
				waveDesc = pyRootPwa.core.waveDescription()
				waveDesc.parseKeyFile(keyFiles[waveName])
				waveDescThres.append( (waveName, waveDesc, float(threshold)) )
			else:
				pyRootPwa.utils.printWarn("cannot parse line '" + line + "' in wave list file "
				          + "'" + waveListFileName + "'.")
#  			if (_debug):
#  				printDebug("reading line " + lineNmb + 1 + ": " + waveName + ", "
#  				           + "threshold = " + threshold + " MeV/c^2")
			lineNmb += 1
	pyRootPwa.utils.printInfo("read " + str(lineNmb) + " lines from wave list file " + "'" + waveListFileName + "'")
	return waveDescThres


# def compareBinningMaps(const map<string, pair<double, double> >& base, const map<string, pair<double, double> >& moreBins)
# {
# 	map<string, pair<double, double> > additionalVars;
# 	typedef map<string, pair<double, double>>::iterator it_type;
# 	for(it_type iterator = moreBins.begin(); iterator != moreBins.end(); iterator++) {
# 		string key = iterator->first;
# 		double lowerBound = iterator->second.first;
# 		double upperBound = iterator->second.second;
# 		it_type currentBase = base.find(key);
# 		if (currentBase == base.end()) {
# 			additionalVars.insert(pair<string, pair<double,double> >(iterator->first, iterator->second));
# 		}
# 		else {
# 			if(lowerBound != currentBase->first or upperBound != currentBase->second) {
# 				printErr << "LALALA" << endl;
# 			}
# 		}
# 	}
# 	return additionalVars;
# }

#def addAmplitudeFromFileNames(likelihood, waveName, ampFileNameList, addBinningMap, eventMetas):
#		if not len(ampFileNameList) == len(eventMetas) and len(eventMetas) > 0:
#			pyRootPwa.utils.printErr("number of amplitude and event files do not match")
#			return False
#		ampMetas = []
#		for iFile, ampFileName in enumerate(ampFileNameList):
#			ampFile = ROOT.TFile.Open(ampFileName, "READ")
#			if not ampFile:
#				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
#				return False
#			ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
#			if not ampMeta:
#				pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
#				return False
#			ampMetas.append(ampMeta)
#		if (not likelihood.addAmplitude(ampMetas, addBinningMap, eventMetas)):
#			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
#			return False
#		return True
def initLikelihood(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, cauchy=False, cauchyWidth=0.5,  accEventsOverride=0, rank=1, verbose=False, addBinningMap=[], evtFileNameList=[], writeFile = "/nfs/mds/user/fkrinner/rootpwa/testWriteLikelihood.root"):
	waveDescThres = readWaveList(waveListFileName, keyFiles)
	massBinCenter = float(addBinningMap['mass'][1] + addBinningMap['mass'][0]) / 2. /1000 # YOU CAN DO BETTER

	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if (not verbose):
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		printErr("could not initialize likelihood. Aborting...")
		return False

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	if len(normIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	normIntMatrix = normIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addNormIntegral(normIntMatrix)):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return False
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	if len(accIntFile.GetListOfKeys()) != 1:
		pyRootPwa.utils.printWarn("'" + normIntegralFileName + "' does not contain exactly one TKey.")
		return False
	accIntMatrix = accIntFile.Get(pyRootPwa.core.ampIntegralMatrix.integralObjectName)
	if (not likelihood.addAccIntegral(accIntMatrix, accEventsOverride)):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return False
	accIntFile.Close()
	eventMetas = []
	openFiles = []
	for evtFileName in evtFileNameList:
		evtFile = ROOT.TFile.Open(evtFileName, "READ")
		if not evtFile:
			pyRootPwa.utils.printErr("could not open amplitude file '" + evtFileName + "'.")
			return False
		evtMeta = pyRootPwa.core.eventMetadata.readEventFile(evtFile)
		if not evtMeta:
			pyRootPwa.utils.printErr("could not get metadata for event file '" + evtFileName + "'.")
			return False
		eventMetas.append(evtMeta)
		openFiles.append(evtFile)


	for wave in waveDescThres:

		waveName = wave[0]
		ampFileNameList = ampFileList[waveName]
		if not len(ampFileNameList) == len(eventMetas) and len(eventMetas) > 0:
			pyRootPwa.utils.printErr("number of amplitude and event files do not match")
			return False
		ampMetas = []
		for  ampFileName in ampFileNameList:
			ampFile = ROOT.TFile.Open(ampFileName, "READ")
			openFiles.append(ampFile)
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
				return False
			ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not ampMeta:
				pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
				return False
			ampMetas.append(ampMeta)

		print "mem before:",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		os.system("free")
		if (not likelihood.addAmplitude(ampMetas, addBinningMap, eventMetas)):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
		print "mem  after:",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		os.system("free")
#		if not writeFile == "":
#			outRoot = pyRootPwa.ROOT.TFile.Open(writeFile,"RECREATE")
#			outRoot.cd()
#			print "write likelihood to",writeFile
#			likelihood.Write("likelihood")
#			outRoot.Close()


	if (not likelihood.finishInit()):
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return False

	for openFile in openFiles:
		openFile.Close()

#	if not writeFile == "":
#		outRoot = pyRootPwa.ROOT.TFile.Open(writeFile,"RECREATE")
#		outRoot.cd()
#		print "write likelihood to",writeFile
#		likelihood.Write("likelihood")
#		outRoot.Close()


	return likelihood




def pwaFit(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, seed=0, cauchy=False, cauchyWidth=0.5, startValFileName="", accEventsOverride=0, checkHessian=False, saveSpace=False, rank=1, verbose=False, addBinningMap=[], evtFileNameList=[]):

	likelihood = initLikelihood(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, cauchy, cauchyWidth, accEventsOverride, rank, verbose, addBinningMap, evtFileNameList)

	lowerBound = addBinningMap[binningMap.keys()[0]][0]
	upperBound = addBinningMap[binningMap.keys()[0]][1]
	fitResult = pyRootPwa.core.pwaFit(likelihood       = likelihood,
	                                  massBinMin       = lowerBound,
	                                  massBinMax       = upperBound,
	                                  seed             = seed,
	                                  startValFileName = startValFileName,
	                                  checkHessian     = checkHessian,
	                                  saveSpace        = saveSpace,
	                                  verbose          = verbose)
	return fitResult


def pwaNloptFit(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, seed=0, cauchy=False, cauchyWidth=0.5, startValFileName="", accEventsOverride=0, checkHessian=False, saveSpace=False, rank=1, verbose=False, addBinningMap=[], evtFileNameList=[]):

	likelihood = initLikelihood(ampFileList, normIntegralFileName, accIntegralFileName, binningMap, waveListFileName, keyFiles, cauchy, cauchyWidth, accEventsOverride, rank, verbose, addBinningMap, evtFileNameList)

	lowerBound = addBinningMap[binningMap.keys()[0]][0]
	upperBound = addBinningMap[binningMap.keys()[0]][1]

	fitResult = pyRootPwa.core.pwaNloptFit(likelihood       = likelihood,
	                                       massBinMin       = lowerBound,
	                                       massBinMax       = upperBound,
	                                       seed             = seed,
	                                       startValFileName = startValFileName,
	                                       checkHessian     = checkHessian,
	                                       saveSpace        = saveSpace,
	                                       verbose          = verbose)
	return fitResult

#!/usr/bin/env python
from copy import copy
import argparse
import sys
sys.path.append("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/build/pyLib")
        
import pyRootPwa
import pyRootPwa.core


def getResults(rootFileName, valTreeName = "pwa"):
	results = []
	print "open",rootFileName	
	inRoot = pyRootPwa.ROOT.TFile.Open(rootFileName)
	tree = inRoot.Get(valTreeName)
	try:
		for result in tree:
				fitResult = result.fitResult_v2
				results.append(copy(fitResult))
		inRoot.Close()
	except RuntimeError:
		print "error in",rootFileName
		return results
	return results

def getBestResult(rootFileName, valTreeName = "pwa"):
        results = getResults(rootFileName, valTreeName)
        if len(results) == 0:
                return None
        bestLike = results[0].logLikelihood()
        bestResult = 0
        for i in range(1,len(results)):
                if results[i].logLikelihood() < bestLike:
                        bestLike = results[i].logLikelihood()
                        bestResult = i
        return results[bestResult]

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
	                                 description="pwa NLopt fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-B", type=int, metavar="addBin", dest="addBinID",default=-1, help="additional bin index") 
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-R", type=str, metavar="path", dest="resultFileName", default="", help="result to calculate the hessian for")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from acceptance integral file)")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
	parser.add_argument("-g", type=str, metavar="integralPath", dest="genIntFilename", default="", help="phase space integral file override")
	parser.add_argument("-a", type=str, metavar="integralPath", dest="accIntFilename", default="", help="acceptance integral file override")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
	parser.add_argument("-s", "--seed", type=int, dest="seed",help="will be ignored")
	args = parser.parse_args()

	pyRootPwa.utils.stdoutisatty = sys.stdout.isatty()
	pyRootPwa.utils.stderrisatty = sys.stderr.isatty()



	printErr  = pyRootPwa.utils.printErr
	printWarn = pyRootPwa.utils.printWarn
	printSucc = pyRootPwa.utils.printSucc
	printInfo = pyRootPwa.utils.printInfo
	printDebug = pyRootPwa.utils.printDebug

	config = pyRootPwa.rootPwaConfig()
	if not config.initialize(args.configFileName):
		printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	if not args.addBinID > -1:
		printErr("no bin given")
		raise Exception
	if args.resultFileName == "":
		printErr("No result file given")
		raise Exception
	oldResult = getBestResult(args.resultFileName)
	binIDs = fileManager.getBinIDfromAddBinID(args.addBinID)

	mMax = 0.
	mMin = float("inf")

	ampFileList = {}	

	for binID in binIDs:
		ampFileListBin = fileManager.getAmplitudeFilePaths(binID, pyRootPwa.core.eventMetadata.REAL)
		if not ampFileListBin:
			printErr("could not retrieve valid amplitude file list. Aborting...")
			sys.exit(1)
		for key in ampFileListBin:
			if not ampFileList.has_key(key):
				ampFileList[key] = []
			ampFileList[key].append(ampFileListBin[key])

		binningMap = fileManager.getBinFromID(binID)
		mMin = min(mMin, binningMap['mass'][0])
		mMax = max(mMax, binningMap['mass'][1])
	
	massBinCenter = (mMin+mMax) / 2.
	addBinningMap = {}
	if fileManager.binned:
		addBinningMap = fileManager.additionalBinning[args.addBinID]
	elif args.addBinID >0:
		printWarn("addBin > 0, but no binning in fileManager")

	psIntegralPath = fileManager.getIntegralFilePathAdditionalID(args.addBinID, pyRootPwa.core.eventMetadata.GENERATED)
	if not args.genIntFilename == "":
		psIntegralPath = args.genIntFilename

	accIntegralPath = fileManager.getIntegralFilePathAdditionalID(args.addBinID, pyRootPwa.core.eventMetadata.ACCEPTED)
	if not args.accIntFilename == "":
		accIntegralPath = args.accIntFilename

	eventFileNames = []
	for binID in binIDs:
		eventFile = fileManager.getDataFile(binID, pyRootPwa.core.eventMetadata.REAL)
		eventFileName = eventFile.dataFileName
		eventFileNames.append(eventFileName)
		printInfo("evtFile: " + eventFileName)

	keyFilesWithAmpIndex = fileManager.getKeyFiles()
	keyFiles = {}
	for waveName in keyFilesWithAmpIndex:
		keyFiles[waveName] = keyFilesWithAmpIndex[waveName][0]


	fitResults = pyRootPwa.addIntegralsHessian(
					oldResult            = oldResult, 
					ampFileList          = ampFileList, 
					normIntegralFileName = psIntegralPath, 
					accIntegralFileName  = accIntegralPath, 
					binningMap           = binningMap, 
					waveListFileName     = args.waveListFileName, 
					keyFiles             = keyFiles, 
					cauchy               = args.cauchyPriors, 
					cauchyWidth          = args.cauchyPriorWidth, 
					accEventsOverride    = args.accEventsOverride, 
					rank                 = args.rank, 
					verbose              = args.verbose, 
					addBinningMap        = addBinningMap, 
					evtFileNameList      = eventFileNames
					                                  )
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "NEW")
	for fitResult in fitResults:
		if (not fitResult):
			printErr("didn't get a valid fit result. skip...")
			continue
		printInfo("writing result to '" + args.outputFileName + "'")
		valTreeName   = "pwa"
		valBranchName = "fitResult_v2"

		if ((not outputFile) or outputFile.IsZombie()):
			printErr("cannot open output file '" + args.outputFileName + "'. Aborting...")
			sys.exit(1)

		tree = outputFile.Get(valTreeName)
		if (not tree):
			printInfo("file '" + args.outputFileName + "' is empty. "
				+ "creating new tree '" + valTreeName + "' for PWA result.")
			tree = pyRootPwa.ROOT.TTree(valTreeName, valTreeName)
			if not fitResult.branch(tree, valBranchName):
				printErr("failed to create new branch '" + valBranchName + "' in file '" + args.outputFileName + "'.")
				sys.exit(1)
		else:
			fitResult.setBranchAddress(tree, valBranchName)
		tree.Fill()
	nmbBytes = tree.Write()
	outputFile.Close()
	if nmbBytes == 0:
		printErr("problems writing integral to TKey 'fitResult' "
		       + "in file '" + args.outputFileName + "'")
		sys.exit(1)
	else:
		printSucc("wrote integral to TKey 'fitResult' "
		        + "in file '" + args.outputFileName + "'")




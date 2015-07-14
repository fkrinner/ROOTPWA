#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
	                                 description="pwa fit executable"
	                                )

	parser.add_argument("outputFileName", type=str, metavar="fileName", help="path to output file")
	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-b", type=int, metavar="#", dest="binID", default=0, help="bin ID of fit (default: 0)")
	parser.add_argument("-B", dest="addBin", action='append', help="additional binning in the form 'binningVariable;lowerBound;upperBound' (e.g. 'mass;1000;1100')."+
	                                                             "You can use the argument multiple times for multiple binning variables")
	parser.add_argument("-s", type=int, metavar="#", dest="seed", default=0, help="random seed (default: 0)")
	parser.add_argument("-C", "--cauchyPriors", help="use half-Cauchy priors (default: false)", action="store_true")
	parser.add_argument("-P", "--cauchyPriorWidth", type=float, metavar ="WIDTH", default=0.5, help="width of half-Cauchy prior (default: 0.5)")
	parser.add_argument("-w", type=str, metavar="path", dest="waveListFileName", default="", help="path to wavelist file (default: none)")
	parser.add_argument("-S", type=str, metavar="path", dest="startValFileName", default="", help="path to start value fit result file (default: none)")
	parser.add_argument("-r", type=int, metavar="#", dest="rank", default=1, help="rank of spin density matrix (default: 1)")
	parser.add_argument("-A", type=int, metavar="#", dest="accEventsOverride", default=0, help="number of input events to normalize acceptance to (default: use number of events from acceptance integral file)")
	parser.add_argument("-H", "--checkHessian", help="check analytical Hessian eigenvalues (default: false)", action="store_true")
	parser.add_argument("-z", "--saveSpace", help="save space by not saving integral and covariance matrices (default: false)", action="store_true")
	parser.add_argument("-g", type=str, metavar="integralPath", dest="genIntFilename", default="", help="phase space integral file override")
	parser.add_argument("-a", type=str, metavar="integralPath", dest="accIntFilename", default="", help="acceptance integral file override")
	parser.add_argument("-v", "--verbose", help="verbose; print debug output (default: false)", action="store_true")
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

	ampFileList = fileManager.getAmplitudeFilePaths(args.binID, pyRootPwa.core.eventMetadata.REAL)
	if not ampFileList:
		printErr("could not retrieve valid amplitude file list. Aborting...")
		sys.exit(1)
	binningMap = fileManager.getBinFromID(args.binID)
	massBinCenter = (binningMap['mass'][1] + binningMap['mass'][0]) / 2.

	psIntegralPath  = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.GENERATED)
	if not args.genIntFilename == "":
		psIntegralPath = args.genIntFilename
	accIntegralPath = fileManager.getIntegralFilePath(args.binID, pyRootPwa.core.eventMetadata.ACCEPTED)
	if not args.accIntFilename == "":
		accIntegralPath = args.accIntFilename

	addBinningMap = pyRootPwa.utils.binningMapFromArgList(args.addBin)
	if not addBinningMap:
		printWarn("received no valid additional binning map argument")
	eventFile = fileManager.getDataFile(args.binID, pyRootPwa.core.eventMetadata.REAL)
	eventFileName = eventFile.dataFileName
	printInfo("evtFile: " + eventFileName)

	keyFilesWithAmpIndex = fileManager.getKeyFiles()
	keyFiles = {}
	for waveName in keyFilesWithAmpIndex:
		keyFiles[waveName] = keyFilesWithAmpIndex[waveName][0]

	fitResult = pyRootPwa.pwaFit(
	                             ampFileList = ampFileList,
	                             normIntegralFileName = psIntegralPath,
	                             accIntegralFileName = accIntegralPath,
	                             binningMap = binningMap,
	                             waveListFileName = args.waveListFileName,
	                             keyFiles = keyFiles,
	                             seed = args.seed,
	                             cauchy = args.cauchyPriors,
	                             cauchyWidth = args.cauchyPriorWidth,
	                             startValFileName = args.startValFileName,
	                             accEventsOverride = args.accEventsOverride,
	                             checkHessian = args.checkHessian,
	                             saveSpace = args.saveSpace,
	                             rank = args.rank,
	                             verbose = args.verbose,
	                             addBinningMap = addBinningMap,
	                             evtFileName = eventFileName
	                             )
	if (not fitResult):
		printErr("didn't get a valid fit result. Aborting...")
		sys.exit(1)
	printInfo("writing result to '" + args.outputFileName + "'")
	valTreeName   = "pwa"
	valBranchName = "fitResult_v2"
	outputFile = pyRootPwa.ROOT.TFile.Open(args.outputFileName, "UPDATE")
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
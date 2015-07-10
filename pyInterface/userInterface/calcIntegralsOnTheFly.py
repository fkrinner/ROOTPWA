#!/usr/bin/env python

import argparse
import sys

import pyRootPwa
import pyRootPwa.core

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="calculate integral matrices")

	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=-1, help="maximum number of events to process (default: all)")
	parser.add_argument("-b", type=int, metavar="massBin", default=-1, dest="massBin", help="mass bin to be calculated (default: all)")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('real', 'generated' or 'accepted', default: all)")
	parser.add_argument("-w", type=str, metavar="path", dest="weightsFileName", default="", help="path to MC weight file for de-weighting (default: none)")
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
		pyRootPwa.utils.printErr("loading config file '" + args.configFileName + "' failed. Aborting...")
		sys.exit(1)
	pyRootPwa.core.particleDataTable.readFile(config.pdgFileName)
	fileManager = pyRootPwa.loadFileManager(config.fileManagerPath)
	if not fileManager:
		pyRootPwa.utils.printErr("loading the file manager failed. Aborting...")
		sys.exit(1)

	binIDList = fileManager.getBinIDList()
	if (not args.massBin == -1):
		binIDList = [args.massBin]

	eventsTypes = []
	if (args.eventsType == "generated"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED ]
	elif (args.eventsType == "accepted"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.ACCEPTED ]
	elif (args.eventsType == "all"):
		eventsTypes = [ pyRootPwa.core.eventMetadata.GENERATED,
		                pyRootPwa.core.eventMetadata.ACCEPTED ]	
	else:
		pyRootPwa.utils.printErr("Invalid events type given ('" + args.eventsType + "'). Aborting...")
		sys.exit(1)

	for binID in binIDList:
		for eventsType in eventsTypes:
			outputFileName = fileManager.getIntegralFilePath(binID, eventsType)
			keyFiles       = fileManager.getKeyFilePaths()
			inputFileNames  = [fileManager.getDataFile(binID, eventsType).dataFileName]

			if not pyRootPwa.core.getIntegralsFromKeyFiles("integral", outputFileName, keyFiles, inputFileNames, args.nEvents):
				printErr("problem calculating the integrals")
			else:
				printSucc("wrote integral matrix to "+outputFileName)

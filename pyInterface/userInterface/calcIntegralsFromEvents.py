#!/usr/bin/env python

import argparse

import sys
sys.path.append("/nfs/hicran/project/compass/analysis/fkrinner/ROOTPWA/build/pyLib")

import pyRootPwa
import pyRootPwa.core

import os
import datetime

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="calculate integral matrices")

	parser.add_argument("-c", type=str, metavar="configFileName", dest="configFileName", default="./rootpwa.config", help="path to config file (default: './rootpwa.config')")
	parser.add_argument("-n", type=int, metavar="#", dest="nEvents", default=-1, help="maximum number of events to process (default: all)")
	parser.add_argument("-b", type=int, metavar="massBin", default=-1, dest="massBin", help="mass bin to be calculated (default: all)")
	parser.add_argument("-B", type=int, metavar="addBinID", default=-1, dest="addBinID", help="additionalBinID")
	parser.add_argument("-e", type=str, metavar="eventsType", default="all", dest="eventsType", help="events type to be calculated ('real', 'generated' or 'accepted', default: all)")
	parser.add_argument("-w", type=str, metavar="path", dest="weightsFileName", default="", help="path to MC weight file for de-weighting (default: none)")
	parser.add_argument("-s", type=int, metavar="startEvent", default=0, dest="startEvent", help="event index to start integration (default: 0)")
	parser.add_argument("-a", type=str, metavar="appendString", default="", dest="appendString", help="string to append to the file name is several are computed simultnaiously, default: '')")
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

	keyFiles = set(fileManager.getKeyFilePaths())
	amplitudes = []
	waveNames  = []
	for key in keyFiles:
		waveDescription = pyRootPwa.core.waveDescription()
		waveDescription.parseKeyFile(key)
		nAmps = waveDescription.nmbAmplitudes()
		for amp in range(nAmps):
			amplitudes.append(waveDescription.constructAmplitude(amp)[1])
			waveNames.append(waveDescription.waveName(amp))
			amplitudes[-1].init()

	try:
		binnings =  fileManager.additionalBinning
	except AttributeError:
		binnings = [] 

	if args.addBinID > -1:
		if len(binnings) == 0:
			raise Exception("addBinID given, but no additional binning found")
		binnings=[binnings[args.addBinID]]

	integrals = []
	if len(binnings) == 0:
		integrals.append(pyRootPwa.core.ampIntegralMatrix())
		integrals[-1].initialize(waveNames)
	else:
		for binn in binnings:
			integrals.append(pyRootPwa.core.ampIntegralMatrix())
			integrals[-1].initialize(waveNames)
	print datetime.datetime.now()
	for binID in binIDList:
		for eventsType in eventsTypes:

			inputFileName  = fileManager.getDataFile(binID, eventsType).dataFileName
			inputFile = pyRootPwa.ROOT.TFile.Open(inputFileName, "READ")
			if not inputFile:
				printWarn("could not open input file '" + inputFileName + "'.")
				sys.exit(1)
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(inputFile, True)
			if not eventMeta:
				printWarn("could not read metadata from input file '" + inputFileName + "'.")
				sys.exit(1)

			if not pyRootPwa.core.calcBinnedIntegralsFromEventTree(eventMeta, amplitudes, integrals, binnings, args.nEvents, args.startEvent):
				printErr("problem calculating the integrals")
				sys.exit(1)
			inputFile.Close()

			if len(binnings) == 0:
				outputFileName = fileManager.getIntegralFilePath(binID, eventsType).replace(".root",args.appendString+".root")
				outROOT =  pyRootPwa.ROOT.TFile.Open(outputFileName,"RECREATE")
				integrals[0].Write("integral")
				outROOT.Close()
			else:
				for additionalBinID in range(len(binnings)):
					outputFileName = fileManager.getIntegralFilePath(binID, eventsType, additionalBinID).replace(".root",args.appendString+".root")
					print "am writing",additionalBinID,outputFileName, len(integrals)
					directory = os.path.dirname(outputFileName)
					if not os.path.isdir(directory):
						os.makedirs(directory)
					outROOT = pyRootPwa.ROOT.TFile.Open(outputFileName,"RECREATE")
					integrals[additionalBinID].Write("integral")
					outROOT.Close()
			printSucc("wrote integral matrices to "+outputFileName)
	print datetime.datetime.now()



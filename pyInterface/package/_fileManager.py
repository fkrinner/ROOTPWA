
import glob
import os
import cPickle as pickle

import pyRootPwa
import pyRootPwa.core
ROOT = pyRootPwa.ROOT

def saveFileManager(fileManagerObject, path):
	if not os.path.isfile(path):
		try:
			pickle.dump(fileManagerObject, open(path, "wb"))
		except:
			pyRootPwa.utils.printErr("error saving file manager.")
			return False
		return True
	else:
		pyRootPwa.utils.printErr("cannot open file manager file '" + path + "'. File already exists.")
		return False


def loadFileManager(path):
	try:
		fileManagerObject = pickle.load(open(path, "rb"))
	except IOError:
		pyRootPwa.utils.printErr("error loading file manager. File manager not found at '" + path + "'")
		return None
	except:
		pyRootPwa.utils.printErr("error loading file manager.")
		raise
	if not fileManagerObject.areKeyFilesSynced():
		pyRootPwa.utils.printErr("key files are not the same as in file manager.")
		return None
	if not fileManagerObject.areDataFilesSynced():
		pyRootPwa.utils.printErr("data files are not the same as in file manager.")
		return None
	return fileManagerObject


class InputFile:
	dataFileName = ""
	binningMap = {}
	eventsType = None
	def __init__(self, dataFileName, binningMap, eventsType):
		self.dataFileName = dataFileName
		self.binningMap = binningMap
		self.eventsType = eventsType


class EventsType:
	OTHER     = 0
	REAL      = 1
	GENERATED = 2
	ACCEPTED  = 3


class fileManager:
	dataDirectory      = ""
	keyDirectory       = ""
	amplitudeDirectory = ""
	integralDirectory  = ""
	limitFilesInDir = -1

	dataFiles       = {}
	keyFiles        = {}
	amplitudeFiles  = {}
	intergalFiles   = {}
	globalAxes = {}
	binList = []


	def initialize(self, configObject):
		self.dataDirectory      = configObject.dataDirectory
		self.keyDirectory       = configObject.keyDirectory
		self.amplitudeDirectory = configObject.ampDirectory
		self.integralDirectory  = configObject.intDirectory
		try:
			self.binned = True
			self.additionalBinningFile= configObject.binFile
		except AttributeError:
			self.binned = False

		pyRootPwa.utils.printInfo("data file dir read from config file: '" + self.dataDirectory + "'.")
		pyRootPwa.utils.printInfo("key file dir read from config file: '" + self.keyDirectory + "'.")
		pyRootPwa.utils.printInfo("amplitude file dir read from config file: '" + self.amplitudeDirectory + "'.")
		pyRootPwa.utils.printInfo("integral file dir read from config file: '" + self.integralDirectory + "'.")
		if self.binned:
			pyRootPwa.utils.printInfo("additional binning file is: '"+self.additionalBinningFile+"'.")

		self.dataFiles = self._openDataFiles()
		allAxes = []
		for eventsType in self.dataFiles:
			allAxes.append(self._getBinningAxes(self.dataFiles[eventsType]))
		self.globalAxes = self._combineAxes(allAxes)
		self.binList = self._createBinIDs()
		self.keyFiles = self._openKeyFiles()

		if self.binned:
			self.additionalBinning = self._openAdditionalBinning()

		if len(self.keyFiles) == 0:
			pyRootPwa.utils.printErr("error loading keyfiles.")
			return False

		if int(configObject.limitFilesPerDir) == 0:
			self.limitFilesInDir = len(self.keyFiles) * len(self.dataFiles)
			pyRootPwa.utils.printInfo("limit for files per directory set to " + str(self.limitFilesInDir))
		elif int(configObject.limitFilesPerDir) == -1:
			self.limitFilesInDir = -1
			pyRootPwa.utils.printInfo("limit for files per directory set to infinite")
		else:
			self.limitFilesInDir = int(configObject.limitFilesPerDir)
			pyRootPwa.utils.printInfo("limit for files per directory set to " + str(self.limitFilesInDir))
		self.amplitudeFiles = self._getAmplitudeFilePaths()
		pyRootPwa.utils.printInfo("number of amplitude files: " + str(len(self.amplitudeFiles)))
		self.integralFiles = self._getIntegralFilePaths()
		pyRootPwa.utils.printInfo("number of integral files: " + str(len(self.integralFiles)))
		self.binIDmap = self._initBinIDmap()
		return True


	def getMissingBins(self):
		missingBins = []
		for binID in self.getBinIDList():
			currentBin = self.binList[binID]
			for eventsType in self.dataFiles:
				found = False
				for dataFile in self.dataFiles[eventsType]:
					if dataFile.binningMap == currentBin:
						found = True
						break
				if not found:
					missingBins.append(binID)
		return missingBins


	def getDataFilePaths(self):
		return fileManager.convertDataFilesToPaths(self.dataFiles)


	def getDataFile(self, binID, eventsType):
		internalEventsType = fileManager.eventsTypeFromBpEnum(eventsType)
		if not internalEventsType in self.dataFiles:
			pyRootPwa.utils.printErr("did not find data file with eventsType '" + str(eventsType) + " (internal '"+str(internalEventsType)+").")
			return False
		for dataFile in self.dataFiles[internalEventsType]:
			found = True
			if not dataFile.binningMap.keys() == self.binList[binID].keys():
				continue
			for binningVariable in dataFile.binningMap:
				if not self.binList[binID][binningVariable] == dataFile.binningMap[binningVariable]:
					found = False
					break
			if found: 
				return dataFile
		pyRootPwa.utils.printWarn("no data dataFile found for binID = " + str(binID) + " and eventsType = '" + str(eventsType) + "'.")
		return False


	def getKeyFilePaths(self):
		return fileManager.convertKeyFilesToPaths(self.keyFiles)


	def getKeyFiles(self):
		return self.keyFiles


	def getKeyFile(self, waveName):
		return self.keyFiles[waveName]


	def getAmplitudeFilePaths(self, binID, eventsType):
		ampFileList = {}
		for key in sorted(self.amplitudeFiles):
			if (key[0] == binID and key[2] == eventsType):
				ampFileList[key[1]] = self.amplitudeDirectory + "/" + self.amplitudeFiles[key]
		return ampFileList


	def getAmplitudeFilePath(self, binID, waveName, eventsType):
		return self.amplitudeDirectory + "/" + self.amplitudeFiles[(binID, waveName, fileManager.eventsTypeFromBpEnum(eventsType))]


	def getIntegralFilePath(self, binID, eventsType, additionalBinID = -1):
		if additionalBinID == -1:
			return self.integralFiles[(binID, fileManager.eventsTypeFromBpEnum(eventsType))]
		else:
			return self.integralFiles[(binID, fileManager.eventsTypeFromBpEnum(eventsType), additionalBinID)]


	def getBinID(self, binInformation):
		foundBins = []
		for binID in range(len(self.binList)):
			found = True
			currentBin = self.binList[binID]
			for binningVariable in binInformation:
				lower = currentBin[binningVariable][0]
				upper = currentBin[binningVariable][1]
				given = binInformation[binningVariable]
				if lower > given or upper < given:
					found = False
					break
			if found: foundBins.append(binID)
		return foundBins


	def getBinFromID(self, binID):
		if binID not in self.getBinIDList():
			pyRootPwa.utils.printErr("binID not found: " + str(binID) + ".")
			raise Exception("do this properly")
		return self.binList[binID]

	def getBinIDfromAddBinID(self, addBinID):
		if not self.binned:
			pyRootPwa.utils.printErr("no additional binning available. set in the config file")
			raise Exception
		return self.binIDmap[addBinID]

	def getAddBinIDfromBinID(self, binID):
		if not self.binned:
			pyRootPwa.utils.printErr("no additional binning available. set in the config file")
			raise Exception
		addIDs = []
		for addID in self.binIDmap:
			if binID in self.binIDmap[addID]:
				addIDs.append(addID)
		return addIDs

	def getActiveBins(self):
		active = []
		for binID in self.getBinIDList():
			if len(self.getAddBinIDfromBinID(binID)):
				active.append(binID)
		return active

	def getBinIDList(self):
		return range(len(self.binList))


	def getWaveNameList(self):
		return self.keyFiles.keys()

	def _initBinIDmap(self):
		binIDmap = {}
		for addID, addBin in enumerate(self.additionalBinning):
			binIDmap[addID] = []
			mMin,mMax = addBin["mass"]
			for binID, bin in enumerate(self.binList):
				binMin,binMax = bin["mass"]
				binMin = float(binMin)/1000.
				binMax = float(binMax)/1000.
				if (not binMin >= mMax) and (not binMax <= mMin):
					binIDmap[addID].append(binID)
		return binIDmap

	def _getBinningAxes(self, fileList):
		if not fileList:
			pyRootPwa.utils.printErr("got no input files")
			raise Exception("do this properly")
		binAxes = {}
		for binningVariable in fileList[0].binningMap:
			binAxes[binningVariable] = []
		for inputFile in fileList:
			if inputFile.binningMap.keys() != binAxes.keys():
				pyRootPwa.utils.printErr("data file '" + inputFile.dataFileName +
				                         "' seems to have different binning variables.")
				return {}
			for key in inputFile.binningMap:
				if inputFile.binningMap[key] not in binAxes[key]:
					binAxes[key].append(inputFile.binningMap[key])
		for binningVariable in binAxes:
			binAxes[binningVariable] = sorted(binAxes[binningVariable], key=lambda t: t[0])
			bins = binAxes[binningVariable]
			if not bins:
				pyRootPwa.util.printErr("no bins found for variable '" + binningVariable + "'.")
				return {}
			bound = bins[0]
			for currentBin in bins[1:]:
				if bound[1] < currentBin[0]:
					pyRootPwa.utils.printWarn("gap in bin structure found for binned variable '" +
					                          binningVariable + "' between bin '" + str(bound) +
					                          "' and bin '" + str(currentBin) + "'.")
				bound = currentBin
		return binAxes


	def _combineAxes(self, axes):
		if not axes:
			pyRootPwa.utils.printErr("got no axes.")
			return {}
		binningVariables = axes[0].keys()
		for axis in axes:
			if not axis.keys() == binningVariables:
				pyRootPwa.utils.printWarn("found different binning variables in different axes.")
		globalAxes = {}
		for binningVariable in binningVariables:
			binAxes = [ axis[binningVariable] for axis in axes if binningVariable in axis.keys() ]
			globalBinList = []
			for axis in binAxes:
				for currentBin in axis:
					if currentBin not in globalBinList:
						globalBinList.append(currentBin)
			globalBinList = sorted(globalBinList, key=lambda t: t[0])
			if not globalBinList:
				pyRootPwa.util.printErr("global currentBin list empty for '" + binningVariable + "'.")
				return {}
			bound = globalBinList[0]
			for currentBin in globalBinList[1:]:
				if bound[1] > currentBin[0]:
					pyRootPwa.utils.printErr("overlap in bin structure found for binned variable '" +
					                         binningVariable + "' between bin '" + str(bound) +
					                         "' and bin '" + str(currentBin) + "'.")
					return {}
				bound = currentBin
			globalAxes[binningVariable] = globalBinList
		pyRootPwa.utils.printSucc("combined currentBin axes: " + str(globalAxes.keys()))
		return globalAxes


	def _getAmplitudeFilePaths(self):
		amplitudeFiles = {}
		overLimit = self.isFilePerDirLimitReached()
		fileCount = 0
		dirSet = []
		for binID in self.getBinIDList():
			for waveName in self.keyFiles:
				for eventsType in self.dataFiles:
					if overLimit:
						subDir = fileCount/self.limitFilesInDir
						fullDir = self.amplitudeDirectory + "/" + str(subDir)
						if not fullDir in dirSet:
							dirSet.append(fullDir)
						amplitudeFiles[(binID, waveName, eventsType)] = str(subDir) + "/" + waveName + "_binID-" + str(binID) + "_" + str(eventsType) + ".root"
						fileCount += 1
					elif not overLimit:
						amplitudeFiles[(binID, waveName, eventsType)] = waveName + "_binID-" + str(binID) + "_" + str(eventsType) + ".root"

		# make sure directories exist (create if neccessary) and check if they are empty
		if overLimit:
			for subDir in dirSet:
				if not os.path.isdir(subDir):
					os.mkdir(subDir)
					pyRootPwa.utils.printInfo("created folder for amplitude files: '" + subDir + "'.")
				else:
					if not os.listdir(subDir) == []:
						pyRootPwa.utils.printWarn("directory '" + subDir + "' is not empty.")
		else:
			if not os.listdir(self.amplitudeDirectory) == []:
				pyRootPwa.utils.printWarn("directory '" + self.amplitudeDirectory + "' is not empty.")

		return amplitudeFiles

	def getIntegralFilePathAdditionalID(self, addBinID, eventsType):
		return self.integralDirectory+os.sep+self._getBinningString(addBinID)+"_"+str(eventsType)+".root"


	def _getIntegralFilePaths(self):
		integralFiles = {}
		for binID in self.getBinIDList():
			for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
				if self.binned:
					for additionalBinID in range(len(self.additionalBinning)):
						additionalBin = self._getBinningString(additionalBinID)
						integralFiles[(binID, eventsType, additionalBinID)] = self.integralDirectory + os.sep + "integral_binID-"+str(binID) + os.sep + additionalBin+"_"+str(eventsType)+".root"
				else:
					integralFiles[(binID, eventsType)] = self.integralDirectory + "/integral_binID-" + str(binID) + "_" + str(eventsType) + ".root"
		return integralFiles


	def _getBinningString(self,additionalBinID):
		bin = self.additionalBinning[additionalBinID]
		signleVariables = [key+":"+str(bin[key][0])+"-"+str(bin[key][1]) for key in bin]
		return '_'.join(signleVariables)



	def _openDataFiles(self):
		dataFileNames = glob.glob(self.dataDirectory + "/*.root")
		inputFiles = {}

		for dataFileName in dataFileNames:
			dataFile = ROOT.TFile.Open(dataFileName, "READ")
			if not dataFile:
				pyRootPwa.utils.printErr("could not open event file '" + dataFileName + "'.")
				return {}
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(dataFile)
			if not eventMeta:
				pyRootPwa.utils.printErr("could not find metadata in event file '" + dataFileName + "'.")
				return {}
			inputFile = InputFile(dataFileName, eventMeta.binningMap(), fileManager.eventsTypeFromBpEnum(eventMeta.eventsType()))
			if inputFile.eventsType not in inputFiles:
				inputFiles[inputFile.eventsType] = []
			inputFiles[inputFile.eventsType].append(inputFile)
			dataFile.Close()

		if not inputFiles:
			pyRootPwa.utils.printErr("no binning maps found.")
			return {}
		for eventsType in inputFiles:
			if not inputFiles[eventsType]:
				pyRootPwa.utils.printErr("no binning maps found for eventsType '" + eventsType + "'.")
				return {}

		return inputFiles


	def _openKeyFiles(self):
		keyFileNames = glob.glob(self.keyDirectory + "/*.key")
		keyFiles = {}

		for keyFileID in range(len(keyFileNames)):
			keyFileName = keyFileNames[keyFileID]
			waveDescription = pyRootPwa.core.waveDescription()
			waveDescription.parseKeyFile(keyFileName)
			nmbAmps = waveDescription.nmbAmplitudes() # <<<<<
			for nAmp in range(nmbAmps): # <<<<<
				(success, amplitude) = waveDescription.constructAmplitude(nAmp) # <<<<<
				if not success:
					pyRootPwa.utils.printErr("could not construct decay topology for key file '" + keyFileName + "'.")
					return []
				waveName = waveDescription.waveNameFromTopology(amplitude.decayTopology())
				if waveName in keyFiles:
					pyRootPwa.utils.printErr("duplicate wave name ('" + waveName + "' from files '" + keyFiles[waveName][0] + "' and '" + keyFileName + "'.")
					return []
				keyFiles[waveName] = (keyFileName,nAmp) # <<<<<
		return keyFiles

	def _openAdditionalBinning(self):
		binning = []
		with open(self.additionalBinningFile) as binningIn:
			for line in binningIn.readlines():
				if not line.strip() == "":
					bin = {}
					for chunk in line.split():
						bin[chunk.split(';')[0]] = (float(chunk.split(';')[1]),float(chunk.split(';')[2]))
					binning.append(bin)
		return binning

	def _createBinIDs(self):
		binList = self._iterateBins(0, {}, [])
		return binList


	def _iterateBins(self, dim, currentBin, result):
		keys = self.globalAxes.keys()
		if dim < len(self.globalAxes):
			for newBin in self.globalAxes[keys[dim]]:
				newCurrentBin = currentBin.copy()
				newCurrentBin[keys[dim]] = newBin
				result = self._iterateBins(dim+1, newCurrentBin, result)
		elif dim == len(self.globalAxes):
			result.append(currentBin)
			return result
		else: raise Exception("do this properly")
		return result


	def __repr__(self):
		retStr = "keyfiles:\n"
		for waveName in self.keyFiles:
			retStr += waveName + " >> " + self.keyFiles[waveName] + "\n"
		retStr += "\ndatafiles:\n"
		for eventsType in self.dataFiles:
			for dataFile in self.dataFiles[eventsType]:
				retStr += "eventsType [" + str(eventsType) + "], bin [" + str(dataFile.binningMap) + "] >> " + dataFile.dataFileName + "\n"
		retStr += "\nampfiles:\n"
		for eventsType in self.dataFiles:
			for binID in self.getBinIDList():
				for waveName in self.keyFiles:
					retStr += "eventsType [" + str(eventsType) + "], binID [" + str(binID) + "], wavename [" + waveName + "] >> " + self.amplitudeFiles[(binID, waveName, eventsType)] + "\n"
		retStr += "\nintfiles:\n"
		for eventsType in [EventsType.GENERATED, EventsType.ACCEPTED]:
			for binID in self.getBinIDList():
				retStr += "eventsType [" + str(eventsType) + "], binID [" + str(binID) + "] >> " + self.integralFiles[(binID, eventsType)] + "\n"
		return retStr


	def isFilePerDirLimitReached(self):
		return len(self.binList) * len(self.keyFiles) * len(self.dataFiles) > self.limitFilesInDir and not self.limitFilesInDir == -1


	def areDataFilesSynced(self):
		return set(self.getDataFilePaths()) == set(glob.glob(self.dataDirectory + "/*.root"))


	def areKeyFilesSynced(self):
		return set(self.getKeyFilePaths()) == set(glob.glob(self.keyDirectory + "/*.key"))


	def areFilesSynced(self):
		return self.areDataFilesSynced() and self.areKeyFilesSynced()


	@staticmethod
	def convertDataFilesToPaths(dataFilesDict):
		allDataFiles = []
		for eventsType in dataFilesDict:
			for dataFile in dataFilesDict[eventsType]:
				allDataFiles.append(dataFile.dataFileName)
		return allDataFiles


	@staticmethod
	def convertKeyFilesToPaths(keyFilesList):
		allKeyFiles = []
		for keyFile in keyFilesList:
			allKeyFiles.append(keyFilesList[keyFile][0]) # <<<<<
		return allKeyFiles


	@staticmethod
	def eventsTypeFromBpEnum(eventsTypeBP):
		if eventsTypeBP == pyRootPwa.core.eventMetadata.OTHER:
			return EventsType.OTHER
		elif eventsTypeBP == pyRootPwa.core.eventMetadata.REAL:
			return EventsType.REAL
		elif eventsTypeBP == pyRootPwa.core.eventMetadata.GENERATED:
			return EventsType.GENERATED
		elif eventsTypeBP == pyRootPwa.core.eventMetadata.ACCEPTED:
			return EventsType.ACCEPTED
		else:
			return -1


	@staticmethod
	def eventsTypeToBpEnum(eventsType):
		if eventsType == EventsType.OTHER:
			return pyRootPwa.core.eventMetadata.OTHER
		elif eventsType == EventsType.REAL:
			return pyRootPwa.core.eventMetadata.REAL
		elif eventsType == EventsType.GENERATED:
			return pyRootPwa.core.eventMetadata.GENERATED
		elif eventsType == EventsType.ACCEPTED:
			return pyRootPwa.core.eventMetadata.ACCEPTED
		else:
			return -1

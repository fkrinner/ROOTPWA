import pyRootPwa.utils
ROOT = pyRootPwa.utils.ROOT

def initLikelihood(waveDescThres,
                   massBinCenter,
                   ampFileList,
                   normIntegralFileName,
                   accIntegralFileName,
                   binningMap    = {},
                   eventFileNames = [],
                   accEventsOverride = 0,
                   cauchy = False,
                   cauchyWidth = 0.5,
                   rank = 1,
                   verbose = False
                  ):
	likelihood = pyRootPwa.core.pwaLikelihood()
	likelihood.useNormalizedAmps(True)
	if not verbose:
		likelihood.setQuiet()
	if cauchy:
		likelihood.setPriorType(pyRootPwa.core.HALF_CAUCHY)
		likelihood.setCauchyWidth(cauchyWidth)
	if (not likelihood.init(waveDescThres,
	                        rank,
	                        massBinCenter)):
		pyRootPwa.utils.printErr("could not initialize likelihood. Aborting...")
		return None

	normIntFile = ROOT.TFile.Open(normIntegralFileName, "READ")
	normIntMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(normIntFile)
	normIntMatrix = normIntMeta.getAmpIntegralMatrix()
	if not likelihood.addNormIntegral(normIntMatrix):
		pyRootPwa.utils.printErr("could not add normalization integral. Aborting...")
		return None
	normIntFile.Close()
	accIntFile = ROOT.TFile.Open(accIntegralFileName, "READ")
	accIntMeta = pyRootPwa.core.ampIntegralMatrixMetadata.readIntegralFile(accIntFile)
	accIntMatrix = accIntMeta.getAmpIntegralMatrix()
	if not likelihood.addAccIntegral(accIntMatrix, accEventsOverride):
		pyRootPwa.utils.printErr("could not add acceptance integral. Aborting...")
		return None
	accIntFile.Close()
	eventFiles = []
	if not len(eventFileNames) == 0:
		pyRootPwa.utils.printInfo("On the fly binning activated")
		eventMetadata = []
		for eventFileName in eventFileNames:
			eventFile = ROOT.TFile.Open(eventFileName, "READ")
			if not eventFile:
				pyRootPwa.utils.printErr("could not open event file '" + eventFileName + "'.")
        	                return None
			eventFiles.append(eventFile)
			eventMeta = pyRootPwa.core.eventMetadata.readEventFile(eventFile)
			if not eventMeta:
				pyRootPwa.utils.printErr("Could not get eventMetadata")
				return None
			eventMetadata.append(eventMeta)
		if not likelihood.setOnTheFlyBinning(binningMap, eventMetadata):
			pyRootPwa.utils.printErr("Could not set on the fly binning")
			return None
	for wave in waveDescThres:
		waveName = wave[0]
		ampMetas = []
		ampFiles = []
		ampFileNames = ampFileList[waveName]
		for ampFileName in ampFileNames:
			ampFile = ROOT.TFile.Open(ampFileName, "READ")
			if not ampFile:
				pyRootPwa.utils.printErr("could not open amplitude file '" + ampFileName + "'.")
				return None
			ampMeta = pyRootPwa.core.amplitudeMetadata.readAmplitudeFile(ampFile, waveName)
			if not ampMeta:
				pyRootPwa.utils.printErr("could not get metadata for waveName '" + waveName + "'.")
				return None
			ampFiles.append(ampFile)
			ampMetas.append(ampMeta)
		if not likelihood.addAmplitude(ampMetas):
			pyRootPwa.utils.printErr("could not add amplitude '" + waveName + "'. Aborting...")
			return None
		for ampFile in ampFiles:
			ampFile.Close()
	if not likelihood.finishInit():
		pyRootPwa.utils.printErr("could not finish initialization of likelihood. Aborting...")
		return None

	return likelihood

import FWCore.ParameterSet.Config as cms

isData = True
#skim = 0
year = 2016
period = 'Run2016B'
#sampleName = 'MonteCarlo'
# sampleName = 'TTJets', "QCD", "DYJetsToLL_M50"

#if "@SKIM@".lower()=="0":
#    skim = 0

#sampleName = "@SAMPLE_NAME@"

process = cms.Process("TreeProducer")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.load('Configuration/StandardSequences/Services_cff')
process.load("Configuration.Geometry.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v8', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
# JetHT AOD
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/150/00000/FC972EB3-D819-E611-94F9-02163E0134F4.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/02C276F7-DA19-E611-8189-02163E0141DD.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/065F7CCE-DA19-E611-90AD-02163E0125FC.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/086C3FED-DA19-E611-BC4E-02163E012611.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/087EA43F-DA19-E611-8610-02163E011DF8.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0C69B99A-DA19-E611-9EAB-02163E011808.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/0E75A2A8-DA19-E611-8148-02163E011D55.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/10BEC6A7-DA19-E611-A373-02163E01456A.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/12BCFF67-EB19-E611-B2C4-02163E01417D.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/16DA718F-DA19-E611-BCEE-02163E01376E.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/18788F8B-DA19-E611-BA71-02163E01264D.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1A8F7BFA-E819-E611-910F-02163E014683.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1C0996EA-DA19-E611-827F-02163E014129.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/1E09B0CE-DA19-E611-8411-02163E0125FC.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/20E6C9D5-DA19-E611-888D-02163E012B1F.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/26D91FF0-DF19-E611-93F0-02163E013895.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28177EED-DA19-E611-998B-02163E011FDE.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/28EE6368-DA19-E611-9F3F-02163E01456A.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2A6A2E06-DB19-E611-9127-02163E014501.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2ACD0DBF-DA19-E611-98F9-02163E011FAB.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2C3C32B1-DA19-E611-BE23-02163E014151.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/2CCD4200-DB19-E611-8910-02163E01441D.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/32B30915-DC19-E611-9685-02163E0143CF.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34A253B0-DD19-E611-834E-02163E014614.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/34D1220F-DB19-E611-8212-02163E01456A.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/36E4FFBD-DA19-E611-858B-02163E01376E.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/38406C52-DB19-E611-90DC-02163E0139C3.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A56B173-DA19-E611-8616-02163E01456A.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3A8A29DA-DB19-E611-A5CF-02163E01348D.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/3EA5C90C-DB19-E611-8809-02163E0144F0.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/48F10CC8-DB19-E611-BA67-02163E014501.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4AB035DA-DA19-E611-8C38-02163E013917.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4CA6F7C0-DA19-E611-BA60-02163E0138B2.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E24C919-DB19-E611-A7C4-02163E014540.root',
        '/store/data/Run2016B/JetHT/AOD/PromptReco-v2/000/273/158/00000/4E945EE7-DA19-E611-9084-02163E0139C8.root',
        )		    
)
 
# Electron ID ==========================================================================================

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



### END Electron ID ====================================================================================  


###########################################
######   NTUPLE PRODUCER    ###############
###########################################
process.makeroottree = cms.EDAnalyzer("NTupleMakerAOD",
# data, year, period, skim
IsData = cms.untracked.bool(isData),
Year = cms.untracked.uint32(year),
Period = cms.untracked.string(period),
Skim = cms.untracked.uint32(0),
# switches of collections
GenParticles = cms.untracked.bool(False),
RecPrimVertex = cms.untracked.bool(True),
RecBeamSpot = cms.untracked.bool(True),
RecTrack = cms.untracked.bool(False),
RecPhoton = cms.untracked.bool(False),
RecPiZero = cms.untracked.bool(True),
RecSecVertex = cms.untracked.bool(True),
RecJet = cms.untracked.bool(False),
RecV0 = cms.untracked.bool(True),
# JEC
JECfile = cms.untracked.string("DesyTauAnalyses/NTupleMaker/data/Fall15_25nsV2/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt"),
# collections
JetCollectionTag = cms.InputTag(""),
PFCandidateCollectionTag = cms.InputTag("particleFlow"),
GenParticleCollectionTag = cms.InputTag("genParticles"),
BeamSpotCollectionTag =  cms.InputTag("offlineBeamSpot"),
PVCollectionTag = cms.InputTag("offlinePrimaryVertices"),
KshortCollectionTag = cms.InputTag("generalV0Candidates","Kshort","RECO"),
LambdaCollectionTag = cms.InputTag("generalV0Candidates","Lambda","RECO"),
TauPiZeroCollectionTag = cms.InputTag("hpsPFTauProducer","pizeros","RECO"),
SecVertexCollectionTag = cms.InputTag("inclusiveSecondaryVertices"),
# tracks
RecTrackPtMin = cms.untracked.double(1.0),
RecTrackEtaMax = cms.untracked.double(2.4),
RecTrackDxyMax = cms.untracked.double(1.0),
RecTrackDzMax = cms.untracked.double(1.0),
RecTrackNum = cms.untracked.int32(0),
#
RecPiZeroPtMin = cms.untracked.double(1.0),
RecPiZeroEtaMax = cms.untracked.double(2.5),
RecPizeroNum = cms.untracked.int32(0),
# photons
RecPhotonPtMin = cms.untracked.double(1.),
RecPhotonEtaMax = cms.untracked.double(2.5),
RecPhotonNum = cms.untracked.int32(0),
# electrons
RecElectronPtMin = cms.untracked.double(8.),
RecElectronEtaMax = cms.untracked.double(2.6),
RecElectronNum = cms.untracked.int32(0),
# taus
RecTauPtMin = cms.untracked.double(18),
RecTauEtaMax = cms.untracked.double(2.4),                                      
RecTauNum = cms.untracked.int32(0),
# jets
RecJetPtMin = cms.untracked.double(18.),
RecJetEtaMax = cms.untracked.double(5.2),
RecJetNum = cms.untracked.int32(0),
SampleName = cms.untracked.string("Data") 
)

#process.patJets.addBTagInfo = cms.bool(True)
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.p = cms.Path(
#                     process.particleFlowPtrs +
#                     process.makePatElectrons +
#                     process.makePatMuons + 
#                     process.makePatJets +
#	              process.selectPrimaryVertex *
#                     process.patDefaultSequence * 
#                     process.egmGsfElectronIDSequence *
#                     process.mvaTrigV025nsCSA14 * 
#                     process.mvaNonTrigV025nsCSA14 * 
    process.makeroottree
    )






process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root")
				   )


#process.end = cms.EndPath(process.Out*process.TFileService)

#processDumpFile = open('MyRootMaker.dump', 'w')
#print >> processDumpFile, process.dumpPython()





def customise_for_gc(process):
	import FWCore.ParameterSet.Config as cms
	from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper

	try:
		maxevents = __MAX_EVENTS__
		process.maxEvents = cms.untracked.PSet(
			input = cms.untracked.int32(max(-1, maxevents))
		)
	except:
		pass

	# Dataset related setup
	try:
		primaryFiles = [__FILE_NAMES__]
		process.source = cms.Source('PoolSource',
			skipEvents = cms.untracked.uint32(__SKIP_EVENTS__),
			fileNames = cms.untracked.vstring(primaryFiles)
		)
		try:
			secondaryFiles = [__FILE_NAMES2__]
			process.source.secondaryFileNames = cms.untracked.vstring(secondaryFiles)
		except:
			pass
		try:
			lumirange = [__LUMI_RANGE__]
			if len(lumirange) > 0:
				process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(lumirange)
				process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
		except:
			pass
	except:
		pass

	if hasattr(process, 'RandomNumberGeneratorService'):
		randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
		randSvc.populate()

	process.AdaptorConfig = cms.Service('AdaptorConfig',
		enable = cms.untracked.bool(True),
		stats = cms.untracked.bool(True),
	)

	# Generator related setup
	try:
		if hasattr(process, 'generator') and process.source.type_() != 'PoolSource':
			process.source.firstLuminosityBlock = cms.untracked.uint32(1 + __MY_JOBID__)
			print 'Generator random seed:', process.RandomNumberGeneratorService.generator.initialSeed
	except:
		pass

	return (process)

process = customise_for_gc(process)

# grid-control: https://ekptrac.physik.uni-karlsruhe.de/trac/grid-control




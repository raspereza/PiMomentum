#include "PiMomentumScale/NTupleMaker/plugins/NTupleMakerAOD.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisByIntegration.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisBase.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>
#include <DataFormats/METReco/interface/GenMET.h>
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
//#include "AnalysisDataFormats/TauAnalysis/interface/PFMEtSignCovMatrix.h"
//#include "RecoJets/JetProducers/interface/PileupJetIdentifier.h"
//#include "DataFormats/METReco/interface/PFMEtSignCovMatrix.h"
//#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"

#include <DataFormats/TrackReco/interface/Track.h>

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
//#include "TauAnalysis/CandidateTools/interface/CompositePtrCandidateT1T2MEtProducer.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
//#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "PiMomentumScale/CandidateTools/interface/candidateAuxFunctions.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include <TString.h>

using namespace reco;
using namespace pat;
//typedef std::vector<NSVfitEventHypothesisByIntegration> NSVfitEventHypothesisByIntegrationCollection;
typedef ROOT::Math::XYZVector Vector;

// http://root.cern.ch/root/html/ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___.html#ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___:_SMatrix_double_3_3_ROOT__Math__MatRepSym_double_3___
// http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SVectorDoc.html
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;  // Standard Matrix representation for a general 3x3 matrix of type double.
typedef ROOT::Math::SVector<double, 3> SVector3; //// SVector: vector of size 3 

// static const unsigned int SKIM_MUTAUTAU   = (1 << 0);     //1  : mu+tau+tau
// static const unsigned int SKIM_ETAUTAU    = (1 << 1);     //2  : e+tau+tau  
// static const unsigned int SKIM_MUMU       = (1 << 2);     //4  : mu+mu
// static const unsigned int SKIM_EE         = (1 << 3);     //8  : e+e
// static const unsigned int SKIM_ETAU       = (1 << 4);     //16 : e+tau
// static const unsigned int SKIM_ALL        = (1 << 5);     //32 : all
// static const unsigned int SKIM_EMU        = (1 << 6);     //64 : e+mu
// static const unsigned int SKIM_TAUTAU     = (1 << 7);     //128: tau+tau

#include "PiMomentumScale/NTupleMaker/interface/genMatch.h"


//to set the values from parameter
NTupleMakerAOD::NTupleMakerAOD(const edm::ParameterSet& iConfig) :  
  // data, year, period, skim
  cdata(iConfig.getUntrackedParameter<bool>("IsData", false)),
  cYear(iConfig.getUntrackedParameter<unsigned int>("Year")),
  cPeriod(iConfig.getUntrackedParameter<std::string>("Period")),
  cSkim(iConfig.getUntrackedParameter<unsigned int>("Skim")),
  cJECfile(iConfig.getUntrackedParameter<std::string>("JECfile")),
  // switches (collections)
  cgen(iConfig.getUntrackedParameter<bool>("GenParticles", false)),
  cbeamspot(iConfig.getUntrackedParameter<bool>("BeamSpot", false)),
  crecprimvertex(iConfig.getUntrackedParameter<bool>("RecPrimVertex", false)),
  crectrack(iConfig.getUntrackedParameter<bool>("RecTrack", false)),
  crecphoton(iConfig.getUntrackedParameter<bool>("RecPhoton", false)),
  crecpizero(iConfig.getUntrackedParameter<bool>("RecPiZero", false)),
  crecsv(iConfig.getUntrackedParameter<bool>("RecSecVertex", false)),
  crecpfjet(iConfig.getUntrackedParameter<bool>("RecJet", false)),
  crecv0(iConfig.getUntrackedParameter<bool>("RecV0", false)),
  // tracks
  cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.5)),
  cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 2.4)),
  cTrackDxyMax(iConfig.getUntrackedParameter<double>("RecTrackDxyMax", 2.0)),
  cTrackDzMax(iConfig.getUntrackedParameter<double>("RecTrackDzMax", 2.0)),
  cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
  // photons
  cPhotonPtMin(iConfig.getUntrackedParameter<double>("RecPhotonPtMin", 1.)),
  cPhotonEtaMax(iConfig.getUntrackedParameter<double>("RecPhotonEtaMax", 2.5)),
  cPhotonNum(iConfig.getUntrackedParameter<int>("RecPhotonNum", 0)),
  // pzeros
  cPizeroPtMin(iConfig.getUntrackedParameter<double>("RecPiZeroPtMin", 30.)),
  cPizeroEtaMax(iConfig.getUntrackedParameter<double>("RecPiZeroEtaMax", 4.5)),
  cPizeroNum(iConfig.getUntrackedParameter<int>("RecPizeroNum", 0)),
  // jets
  cJetPtMin(iConfig.getUntrackedParameter<double>("RecJetPtMin", 30.)),
  cJetEtaMax(iConfig.getUntrackedParameter<double>("RecJetEtaMax", 4.5)),
  cJetNum(iConfig.getUntrackedParameter<int>("RecJetNum", 0)),
  // sample name
  sampleName(iConfig.getUntrackedParameter<std::string>("SampleName", "Higgs"))
{

  KshortCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("KshortCollectionTag"));
  LambdaCollectionToken_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LambdaCollectionTag"));
  PFCandidateCollectionToken_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandidateCollectionTag"));
  TauPiZeroCollectionToken_ = consumes<reco::RecoTauPiZeroCollection>(iConfig.getParameter<edm::InputTag>("TauPiZeroCollectionTag"));
  SVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("SecVertexCollectionTag"));
  JetCollectionToken_ = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("JetCollectionTag"));
  GenParticleCollectionToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag"));
  BeamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag"));
  PVToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PVCollectionTag"));

  if(cYear != 2011 && cYear != 2012 && cYear != 2015 && cYear!= 2016)
    throw cms::Exception("NTupleMakerAOD") << "Invalid Year, only 2011, 2012, 2015 and 2016 are allowed!";
  if(cPeriod != "Summer11" && cPeriod != "Fall11" && cPeriod != "Summer12" && cPeriod != "PHYS14" && cPeriod != "Spring15" && cPeriod != "Run2015B" && cPeriod != "Run2015C" && cPeriod != "Run2015D" && cPeriod != "Run2016B" && cPeriod != "Spring16")
    throw cms::Exception("NTupleMakerAOD") << "Invalid period, only Summer11, Fall11, Summer12, PHYS14, Spring15, Spring16, Run2015B, Run2015C, Run2015D and Run2016B are allowed!";
  
  double barrelRadius = 129.;  //p81, p50, ECAL TDR
  double endcapZ      = 320.5; // fig 3.26, p81, ECAL TDR
  Surface::RotationType rot;
  ecalBarrel         = Cylinder::build(Surface::PositionType(0, 0, 0), rot, barrelRadius);
  ecalNegativeEtaEndcap = Plane::build(Surface::PositionType(0, 0, -endcapZ), rot);
  ecalPositiveEtaEndcap = Plane::build(Surface::PositionType(0, 0, endcapZ), rot);
  
  const char* cmssw_base = getenv("CMSSW_BASE");
  if(!cmssw_base) throw cms::Exception("No CMSSW runtime environment found");
  std::string prefix = std::string(cmssw_base) + "/src/";
}

//destructor
NTupleMakerAOD::~NTupleMakerAOD(){
  
}


void NTupleMakerAOD::beginJob(){
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("AC1B", "AC1B", 1);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);
  
  tree->Branch("event_nr", &event_nr, "event_nr/i");
  tree->Branch("event_run", &event_run, "event_run/i");
  tree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
  tree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
  tree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
  tree->Branch("errors", &errors, "errors/i");
  tree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
  tree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
  tree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");
  
  // beam spot
  if (cbeamspot) {
    tree->Branch("beamspot_x", &beamspot_x, "beamspot_x/F");
    tree->Branch("beamspot_y", &beamspot_y, "beamspot_y/F");
    tree->Branch("beamspot_z", &beamspot_z, "beamspot_z/F");
    tree->Branch("beamspot_xwidth", &beamspot_xwidth, "beamspot_xwidth/F");
    tree->Branch("beamspot_ywidth", &beamspot_ywidth, "beamspot_ywidth/F");
    tree->Branch("beamspot_zsigma", &beamspot_zsigma, "beamspot_zsigma/F");
    tree->Branch("beamspot_cov", &beamspot_cov, "beamspot_cov[6]/F");
  }  

  // primary vertex
  if (crecprimvertex) {
    tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i"); 
    tree->Branch("goodprimvertex_count", &goodprimvertex_count, "goodprimvertex_count/i"); 
    tree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
    tree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
    tree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
    tree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
    tree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
    tree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_pdf/F");
    tree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
    tree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");
  }  

  // pf jets
  if (crecpfjet) {
    tree->Branch("pfjet_count", &pfjet_count, "pfjet_count/i");
    tree->Branch("pfjet_e", pfjet_e, "pfjet_e[pfjet_count]/F");
    tree->Branch("pfjet_px", pfjet_px, "pfjet_px[pfjet_count]/F");
    tree->Branch("pfjet_py", pfjet_py, "pfjet_py[pfjet_count]/F");
    tree->Branch("pfjet_pz", pfjet_pz, "pfjet_pz[pfjet_count]/F");
    tree->Branch("pfjet_pt", pfjet_pt, "pfjet_pt[pfjet_count]/F");
    tree->Branch("pfjet_eta", pfjet_eta, "pfjet_eta[pfjet_count]/F");
    tree->Branch("pfjet_phi", pfjet_phi, "pfjet_phi[pfjet_count]/F");
    tree->Branch("pfjet_neutralhadronicenergy", pfjet_neutralhadronicenergy, "pfjet_neutralhadronicenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedhadronicenergy", pfjet_chargedhadronicenergy, "pfjet_chargedhadronicenergy[pfjet_count]/F");
    tree->Branch("pfjet_neutralemenergy", pfjet_neutralemenergy, "pfjet_neutralemenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedemenergy", pfjet_chargedemenergy, "pfjet_chargedemenergy[pfjet_count]/F");
    tree->Branch("pfjet_muonenergy", pfjet_muonenergy, "pfjet_muonenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedmuonenergy", pfjet_chargedmuonenergy, "pfjet_chargedmuonenergy[pfjet_count]/F");
    tree->Branch("pfjet_chargedmulti", pfjet_chargedmulti, "pfjet_chargedmulti[pfjet_count]/i");	
    tree->Branch("pfjet_neutralmulti", pfjet_neutralmulti, "pfjet_neutralmulti[pfjet_count]/i");	
    tree->Branch("pfjet_chargedhadronmulti", pfjet_chargedhadronmulti, "pfjet_chargedhadronmulti[pfjet_count]/i");
    //    tree->Branch("pfjet_energycorr", pfjet_energycorr, "pfjet_energycorr[pfjet_count]/F");
    //    tree->Branch("pfjet_energycorr_l1fastjet", pfjet_energycorr_l1fastjet, "pfjet_energycorr_l1fastjet[pfjet_count]/F");
    //    tree->Branch("pfjet_energycorr_l2relative", pfjet_energycorr_l2relative, "pfjet_energycorr_l2relative[pfjet_count]/F");
    //    tree->Branch("pfjet_energycorr_l3absolute", pfjet_energycorr_l3absolute, "pfjet_energycorr_l3absolute[pfjet_count]/F");
    //    tree->Branch("pfjet_energycorr_l2l3residual", pfjet_energycorr_l2l3residual, "pfjet_energycorr_l2l3residual[pfjet_count]/F");
    //    tree->Branch("pfjet_pu_jet_cut_loose", pfjet_pu_jet_cut_loose, "pfjet_pu_jet_cut_loose[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_cut_medium", pfjet_pu_jet_cut_medium, "pfjet_pu_jet_cut_medium[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_cut_tight", pfjet_pu_jet_cut_tight, "pfjet_pu_jet_cut_tight[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_cut_mva", pfjet_pu_jet_cut_mva, "pfjet_pu_jet_cut_mva[pfjet_count]/F");
    //    tree->Branch("pfjet_pu_jet_simple_loose", pfjet_pu_jet_simple_loose, "pfjet_pu_jet_simple_loose[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_simple_medium", pfjet_pu_jet_simple_medium, "pfjet_pu_jet_simple_medium[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_simple_tight", pfjet_pu_jet_simple_tight, "pfjet_pu_jet_simple_tight[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_simple_mva", pfjet_pu_jet_simple_mva, "pfjet_pu_jet_simple_mva[pfjet_count]/F");
    //    tree->Branch("pfjet_pu_jet_full_loose", pfjet_pu_jet_full_loose, "pfjet_pu_jet_full_loose[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_full_medium", pfjet_pu_jet_full_medium, "pfjet_pu_jet_full_medium[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_full_tight", pfjet_pu_jet_full_tight, "pfjet_pu_jet_full_tight[pfjet_count]/O");
    //    tree->Branch("pfjet_pu_jet_full_mva", pfjet_pu_jet_full_mva, "pfjet_pu_jet_full_mva[pfjet_count]/F");
    //    tree->Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour[pfjet_count]/I");
    //    tree->Branch("pfjet_btag", pfjet_btag,"pfjet_btag[pfjet_count][10]/F");
    //    tree->Branch("pfjet_jecUncertainty",pfjet_jecUncertainty,"pfjet_jecUncertainty[pfjet_count]/F");
  }    

  if (crecphoton) {
    tree->Branch("photon_count", &photon_count, "photon_count/i");
    tree->Branch("photon_px", photon_px, "photon_px[photon_count]/F");
    tree->Branch("photon_py", photon_py, "photon_py[photon_count]/F");
    tree->Branch("photon_pz", photon_pz, "photon_pz[photon_count]/F");
    tree->Branch("photon_pt", photon_pt, "photon_pt[photon_count]/F");
    tree->Branch("photon_e" , photon_e,  "photon_e[photon_count]/F");
    tree->Branch("photon_eta", photon_eta, "photon_eta[photon_count]/F");
    tree->Branch("photon_phi", photon_phi, "photon_phi[photon_count]/F");
    /*
    tree->Branch("photon_e1x5", photon_e1x5, "photon_e1x5[photon_count]/F");
    tree->Branch("photon_e2x5", photon_e2x5, "photon_e2x5[photon_count]/F");
    tree->Branch("photon_e3x3", photon_e3x3, "photon_e3x3[photon_count]/F");
    tree->Branch("photon_e5x5", photon_e5x5, "photon_e5x5[photon_count]/F");
    tree->Branch("photon_sigmaetaeta", photon_sigmaetaeta, "photon_sigmaetaeta[photon_count]/F");
    tree->Branch("photon_sigmaietaieta", photon_sigmaietaieta, "photon_sigmaietaieta[photon_count]/F");
    tree->Branch("photon_ehcaloverecal", photon_ehcaloverecal, "photon_ehcaloverecal[photon_count]/F");
    tree->Branch("photon_ehcaloverecaldepth1", photon_ehcaloverecaldepth1, "photon_ehcaloverecaldepth1[photon_count]/F");
    tree->Branch("photon_ehcaloverecaldepth2", photon_ehcaloverecaldepth2, "photon_ehcaloverecaldepth2[photon_count]/F");
    tree->Branch("photon_maxenergyxtal", photon_maxenergyxtal, "photon_maxenergyxtal[photon_count]/F");
    tree->Branch("photon_isolationr3track", photon_isolationr3track, "photon_isolationr3track[photon_count]/F");
    tree->Branch("photon_isolationr3trackhollow", photon_isolationr3trackhollow, "photon_isolationr3trackhollow[photon_count]/F");
    tree->Branch("photon_isolationr3ntrack", photon_isolationr3ntrack, "photon_isolationr3ntrack[photon_count]/i");
    tree->Branch("photon_isolationr3ntrackhollow", photon_isolationr3ntrackhollow, "photon_isolationr3ntrackhollow[photon_count]/i");
    tree->Branch("photon_isolationr3ecal", photon_isolationr3ecal, "photon_isolationr3ecal[photon_count]/F");
    tree->Branch("photon_isolationr3hcal", photon_isolationr3hcal, "photon_isolationr3hcal[photon_count]/F");
    tree->Branch("photon_isolationr4track", photon_isolationr4track, "photon_isolationr4track[photon_count]/F");
    tree->Branch("photon_isolationr4trackhollow", photon_isolationr4trackhollow, "photon_isolationr4trackhollow[photon_count]/F");
    tree->Branch("photon_isolationr4ntrack", photon_isolationr4ntrack, "photon_isolationr4ntrack[photon_count]/i");
    tree->Branch("photon_isolationr4ntrackhollow", photon_isolationr4ntrackhollow, "photon_isolationr4ntrackhollow[photon_count]/i");
    tree->Branch("photon_isolationr4ecal", photon_isolationr4ecal, "photon_isolationr4ecal[photon_count]/F");
    tree->Branch("photon_isolationr4hcal", photon_isolationr4hcal, "photon_isolationr4hcal[photon_count]/F");
    tree->Branch("photon_superclusterindex", photon_superclusterindex, "photon_superclusterindex[photon_count]/I");
    tree->Branch("photon_info", photon_info, "photon_info[photon_count]/b");
    tree->Branch("photon_gapinfo", photon_gapinfo, "photon_gapinfo[photon_count]/i");
    tree->Branch("photon_conversionbegin", photon_conversionbegin, "photon_conversionbegin[photon_count]/i");
    tree->Branch("photon_superClusterX", photon_superClusterX, "photon_superClusterX[photon_count]/F");
    tree->Branch("photon_superClusterY", photon_superClusterY, "photon_superClusterY[photon_count]/F");
    tree->Branch("photon_superClusterZ", photon_superClusterZ, "photon_superClusterZ[photon_count]/F");
    */
  }  

  if (crectrack) { 
    tree->Branch("track_count", &track_count, "track_count/i");
    tree->Branch("track_px", track_px, "track_px[track_count]/F");
    tree->Branch("track_py", track_py, "track_py[track_count]/F");
    tree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
    tree->Branch("track_pt", track_pt, "track_pt[track_count]/F");
    tree->Branch("track_eta", track_eta, "track_eta[track_count]/F");
    tree->Branch("track_phi", track_phi, "track_phi[track_count]/F");
    tree->Branch("track_charge", track_charge, "track_charge[track_count]/F");
    tree->Branch("track_mass", track_mass, "track_mass[track_count]/F");
    tree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
    tree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
    tree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
    tree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
    tree->Branch("track_ID", track_ID, "track_ID[track_count]/I");
    tree->Branch("track_highPurity", track_highPurity, "track_highPurity[track_count]/O");
  }  

  // V0s
  if (crecv0){
    tree->Branch("v0_count", &v0_count, "v0_count/i");
    tree->Branch("v0_ID", v0_ID, "v0_ID[v0_count]/I");
    tree->Branch("v0_vx", v0_vx, "v0_vx[v0_count]/F");
    tree->Branch("v0_vy", v0_vy, "v0_vy[v0_count]/F");
    tree->Branch("v0_vz", v0_vz, "v0_vz[v0_count]/F");
    tree->Branch("v0_px", v0_px, "v0_px[v0_count]/F");
    tree->Branch("v0_py", v0_py, "v0_py[v0_count]/F");
    tree->Branch("v0_pz", v0_pz, "v0_pz[v0_count]/F");
    tree->Branch("v0_pt", v0_pt, "v0_pt[v0_count]/F");
    tree->Branch("v0_eta", v0_eta, "v0_eta[v0_count]/F");
    tree->Branch("v0_phi", v0_phi, "v0_phi[v0_count]/F");
    tree->Branch("v0_mass", v0_mass, "v0_mass[v0_count]/F");
    tree->Branch("v0_chi2", v0_chi2, "v0_chi2[v0_count]/F");
    tree->Branch("v0_ndof", v0_ndof, "v0_ndof[v0_count]/F");
    tree->Branch("v0_ndaughters", v0_ndaughters, "v0_ndaughters[v0_count]/i");

    tree->Branch("v0_pos_px", v0_pos_px, "v0_pos_px[v0_count]/F");
    tree->Branch("v0_pos_py", v0_pos_py, "v0_pos_py[v0_count]/F");    
    tree->Branch("v0_pos_pz", v0_pos_px, "v0_pos_pz[v0_count]/F");
    tree->Branch("v0_pos_mass", v0_pos_mass, "v0_pos_mass[v0_count]/F");    
    tree->Branch("v0_pos_ID", v0_pos_ID, "v0_pos_ID[v0_count]/I");    
    tree->Branch("v0_pos_pt", v0_pos_pt, "v0_pos_pt[v0_count]/F");
    tree->Branch("v0_pos_eta",v0_pos_eta,"v0_pos_eta[v0_count]/F");
    tree->Branch("v0_pos_phi",v0_pos_phi,"v0_pos_phi[v0_count]/F");

    tree->Branch("v0_neg_px", v0_neg_px, "v0_neg_px[v0_count]/F");
    tree->Branch("v0_neg_py", v0_neg_py, "v0_neg_py[v0_count]/F");    
    tree->Branch("v0_neg_pz", v0_neg_px, "v0_neg_pz[v0_count]/F");
    tree->Branch("v0_neg_mass", v0_neg_mass, "v0_neg_mass[v0_count]/F");    
    tree->Branch("v0_neg_ID", v0_neg_ID, "v0_neg_ID[v0_count]/I");    
    tree->Branch("v0_neg_pt", v0_neg_pt, "v0_neg_pt[v0_count]/F");
    tree->Branch("v0_neg_eta",v0_neg_eta,"v0_neg_eta[v0_count]/F");
    tree->Branch("v0_neg_phi",v0_neg_phi,"v0_neg_phi[v0_count]/F");

  }

  // pizeros
  if (crecpizero) {
    tree->Branch("pizero_count",&pizero_count,"pizero_count/i");
    tree->Branch("pizero_px",pizero_px,"pizero_px[pizero_count]/F");
    tree->Branch("pizero_py",pizero_py,"pizero_py[pizero_count]/F");
    tree->Branch("pizero_pz",pizero_pz,"pizero_pz[pizero_count]/F");
    tree->Branch("pizero_pt",pizero_pt,"pizero_pt[pizero_count]/F");
    tree->Branch("pizero_eta",pizero_eta,"pizero_eta[pizero_count]/F");
    tree->Branch("pizero_phi",pizero_phi,"pizero_phi[pizero_count]/F");
    tree->Branch("pizero_e",pizero_e,"pizero_e[pizero_count]/F");
    tree->Branch("pizero_x",pizero_x,"pizero_x[pizero_count]/F");
    tree->Branch("pizero_y",pizero_y,"pizero_y[pizero_count]/F");
    tree->Branch("pizero_z",pizero_z,"pizero_z[pizero_count]/F");
  }
  
  // SVs
  if (crecsv) {
    tree->Branch("sv_count",&sv_count,"sv_count/i");
    tree->Branch("sv_covxx",sv_covxx,"sv_covxx[sv_count]/F");
    tree->Branch("sv_covxy",sv_covxy,"sv_covxy[sv_count]/F");
    tree->Branch("sv_covxz",sv_covxz,"sv_covxz[sv_count]/F");
    tree->Branch("sv_covyy",sv_covyy,"sv_covyy[sv_count]/F");
    tree->Branch("sv_covyz",sv_covyz,"sv_covyz[sv_count]/F");
    tree->Branch("sv_covzz",sv_covzz,"sv_covzz[sv_count]/F");
    tree->Branch("sv_vx",sv_vx,"sv_vx[sv_count]/F");
    tree->Branch("sv_vy",sv_vy,"sv_vy[sv_count]/F");
    tree->Branch("sv_vz",sv_vz,"sv_vz[sv_count]/F");
    tree->Branch("sv_ntracks",sv_ntracks,"sv_ntracks[sv_count]/i");
    tree->Branch("sv_chi2",sv_chi2,"sv_chi2[sv_count]/F");
    tree->Branch("sv_ndof",sv_ndof,"sv_ndof[sv_count]/i");
    tree->Branch("sv_exceedTrkSize",sv_exceedTrkSize,"sv_exceedTrkSize[sv_count]/O");
    tree->Branch("sv_track_px",sv_track_px,"sv_track_px[sv_count][20]/F");
    tree->Branch("sv_track_py",sv_track_py,"sv_track_py[sv_count][20]/F");
    tree->Branch("sv_track_pz",sv_track_pz,"sv_track_pz[sv_count][20]/F");
    tree->Branch("sv_track_pt",sv_track_px,"sv_track_pt[sv_count][20]/F");
    tree->Branch("sv_track_eta",sv_track_px,"sv_track_eta[sv_count][20]/F");
    tree->Branch("sv_track_phi",sv_track_px,"sv_track_phi[sv_count][20]/F");
  }

  // generator info
  if (cgen && !cdata) {
    tree->Branch("genweight", &genweight, "genweight/F");
    tree->Branch("genid1", &genid1, "genid1/F");
    tree->Branch("genx1", &genx1, "genx1/F");
    tree->Branch("genid2", &genid2, "genid2/F");
    tree->Branch("genx2", &genx2, "genx2/F");
    tree->Branch("genScale", &genScale, "genScale/F");
    
    tree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
    tree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
    tree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
    tree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");

    // generated taus
    tree->Branch("gentau_count", &gentau_count, "gentau_count/i");
    tree->Branch("gentau_e",  gentau_e,  "tau_e[gentau_count]/F");
    tree->Branch("gentau_px", gentau_px, "tau_px[gentau_count]/F");
    tree->Branch("gentau_py", gentau_py, "tau_py[gentau_count]/F");
    tree->Branch("gentau_pz", gentau_pz, "tau_pz[gentau_count]/F");

    tree->Branch("gentau_visible_e",  gentau_visible_e,  "tau_visible_e[gentau_count]/F");
    tree->Branch("gentau_visible_px", gentau_visible_px, "tau_visible_px[gentau_count]/F");
    tree->Branch("gentau_visible_py", gentau_visible_py, "tau_visible_py[gentau_count]/F");
    tree->Branch("gentau_visible_pz", gentau_visible_pz, "tau_visible_pz[gentau_count]/F");

    tree->Branch("gentau_visible_pt",  gentau_visible_pt,  "tau_visible_pt[gentau_count]/F");
    tree->Branch("gentau_visible_eta", gentau_visible_eta, "tau_visible_eta[gentau_count]/F");
    tree->Branch("gentau_visible_phi", gentau_visible_phi, "tau_visible_phi[gentau_count]/F");
    tree->Branch("gentau_visible_mass", gentau_visible_mass, "tau_visible_mass[gentau_count]/F");

    tree->Branch("gentau_visibleNoLep_e",  gentau_visibleNoLep_e,  "tau_visibleNoLep_e[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_px", gentau_visibleNoLep_px, "tau_visibleNoLep_px[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_py", gentau_visibleNoLep_py, "tau_visibleNoLep_py[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_pz", gentau_visibleNoLep_pz, "tau_visibleNoLep_pz[gentau_count]/F");

    tree->Branch("gentau_visibleNoLep_pt",  gentau_visibleNoLep_pt,  "tau_visibleNoLep_pt[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_eta", gentau_visibleNoLep_eta, "tau_visibleNoLep_eta[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_phi", gentau_visibleNoLep_phi, "tau_visibleNoLep_phi[gentau_count]/F");
    tree->Branch("gentau_visibleNoLep_mass", gentau_visibleNoLep_mass, "tau_visibleNoLep_mass[gentau_count]/F");
    
    tree->Branch("gentau_status", gentau_status, "gentau_status[gentau_count]/I");
    tree->Branch("gentau_fromHardProcess", gentau_fromHardProcess, "gentau_fromHardProcess[gentau_count]/I");
    tree->Branch("gentau_fromHardProcessBeforeFSR", gentau_fromHardProcessBeforeFSR, "gentau_fromHardProcessBeforeFSR[gentau_count]/I");
    tree->Branch("gentau_isDecayedLeptonHadron", gentau_isDecayedLeptonHadron, "gentau_isDecayedLeptonHadron[gentau_count]/I");
    tree->Branch("gentau_isDirectHadronDecayProduct", gentau_isDirectHadronDecayProduct, "gentau_isDirectHadronDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isDirectHardProcessTauDecayProduct", gentau_isDirectHardProcessTauDecayProduct, "gentau_isDirectHardProcessTauDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isDirectPromptTauDecayProduct", gentau_isDirectPromptTauDecayProduct, "gentau_isDirectPromptTauDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isDirectTauDecayProduct", gentau_isDirectTauDecayProduct, "gentau_isDirectTauDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isFirstCopy", gentau_isFirstCopy, "gentau_isFirstCopy[gentau_count]/I");
    tree->Branch("gentau_isHardProcess", gentau_isHardProcess, "gentau_isHardProcess[gentau_count]/I");
    tree->Branch("gentau_isHardProcessTauDecayProduct", gentau_isHardProcessTauDecayProduct, "gentau_isHardProcessTauDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isLastCopy", gentau_isLastCopy, "gentau_isLastCopy[gentau_count]/I");
    tree->Branch("gentau_isLastCopyBeforeFSR", gentau_isLastCopyBeforeFSR, "gentau_isLastCopyBeforeFSR[gentau_count]/I");
    tree->Branch("gentau_isPrompt", gentau_isPrompt, "gentau_isPrompt[gentau_count]/I");
    tree->Branch("gentau_isPromptTauDecayProduct", gentau_isPromptTauDecayProduct, "gentau_isPromptTauDecayProduct[gentau_count]/I");
    tree->Branch("gentau_isTauDecayProduct", gentau_isTauDecayProduct, "gentau_isTauDecayProduct[gentau_count]/I");

    tree->Branch("gentau_decayMode",  gentau_decayMode,  "tau_decayMode[gentau_count]/I");
    tree->Branch("gentau_decayMode_name",  gentau_decayMode_name,  "tau_decayMode_name[gentau_count]/C");
    tree->Branch("gentau_mother",gentau_mother,"gentau_mother[gentau_count]/b");

    // generated particles
    tree->Branch("genparticles_lheHt", &genparticles_lheHt, "genparticles_lheHt/F");
    tree->Branch("genparticles_noutgoing", &genparticles_noutgoing, "genparticles_noutgoing/i");
    tree->Branch("genparticles_count", &genparticles_count, "genparticles_count/i");
    tree->Branch("genparticles_e", genparticles_e, "genparticles_e[genparticles_count]/F");
    tree->Branch("genparticles_px", genparticles_px, "genparticles_px[genparticles_count]/F");
    tree->Branch("genparticles_py", genparticles_py, "genparticles_py[genparticles_count]/F");
    tree->Branch("genparticles_pz", genparticles_pz, "genparticles_pz[genparticles_count]/F");
    tree->Branch("genparticles_vx", genparticles_vx, "genparticles_vx[genparticles_count]/F");
    tree->Branch("genparticles_vy", genparticles_vy, "genparticles_vy[genparticles_count]/F");
    tree->Branch("genparticles_vz", genparticles_vz, "genparticles_vz[genparticles_count]/F");
    tree->Branch("genparticles_pdgid", genparticles_pdgid, "genparticles_pdgid[genparticles_count]/I");
    tree->Branch("genparticles_status", genparticles_status, "genparticles_status[genparticles_count]/I");
    tree->Branch("genparticles_info", genparticles_info, "genparticles_info[genparticles_count]/i");

    tree->Branch("genparticles_fromHardProcess", genparticles_fromHardProcess, "genparticles_fromHardProcess[genparticles_count]/I");
    tree->Branch("genparticles_fromHardProcessBeforeFSR", genparticles_fromHardProcessBeforeFSR, "genparticles_fromHardProcessBeforeFSR[genparticles_count]/I");
    tree->Branch("genparticles_isDecayedLeptonHadron", genparticles_isDecayedLeptonHadron, "genparticles_isDecayedLeptonHadron[genparticles_count]/I");
    tree->Branch("genparticles_isDirectHadronDecayProduct", genparticles_isDirectHadronDecayProduct, "genparticles_isDirectHadronDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isDirectHardProcessTauDecayProduct", genparticles_isDirectHardProcessTauDecayProduct, "genparticles_isDirectHardProcessTauDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isDirectPromptTauDecayProduct", genparticles_isDirectPromptTauDecayProduct, "genparticles_isDirectPromptTauDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isDirectTauDecayProduct", genparticles_isDirectTauDecayProduct, "genparticles_isDirectTauDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isFirstCopy", genparticles_isFirstCopy, "genparticles_isFirstCopy[genparticles_count]/I");
    tree->Branch("genparticles_isHardProcess", genparticles_isHardProcess, "genparticles_isHardProcess[genparticles_count]/I");
    tree->Branch("genparticles_isHardProcessTauDecayProduct", genparticles_isHardProcessTauDecayProduct, "genparticles_isHardProcessTauDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isLastCopy", genparticles_isLastCopy, "genparticles_isLastCopy[genparticles_count]/I");
    tree->Branch("genparticles_isLastCopyBeforeFSR", genparticles_isLastCopyBeforeFSR, "genparticles_isLastCopyBeforeFSR[genparticles_count]/I");
    tree->Branch("genparticles_isPrompt", genparticles_isPrompt, "genparticles_isPrompt[genparticles_count]/I");
    tree->Branch("genparticles_isPromptTauDecayProduct", genparticles_isPromptTauDecayProduct, "genparticles_isPromptTauDecayProduct[genparticles_count]/I");
    tree->Branch("genparticles_isTauDecayProduct", genparticles_isTauDecayProduct, "genparticles_isTauDecayProduct[genparticles_count]/I");
 
    tree->Branch("genparticles_mother", genparticles_mother, "genparticles_mother[genparticles_count]/b");
  }    

} //void NTupleMakerAOD::beginJob()

void NTupleMakerAOD::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  run_number = iRun.run();
}

void NTupleMakerAOD::beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumi_run = iLumiBlock.run();
  lumi_block = iLumiBlock.luminosityBlock();
  
  if(cdata)
    {
      //edm::Handle<LumiSummary> lumiSummary;
      //iLumiBlock.getByLabel(edm::InputTag("lumiProducer"), lumiSummary);
      lumi_value = 0;//lumiSummary->avgInsDelLumi();
      lumi_valueerr = 0;//lumiSummary->avgInsDelLumiErr();
      lumi_livefrac = 0;//lumiSummary->lumiSecQual();
      lumi_deadfrac = 0;//lumiSummary->deadFrac();
      lumi_quality = 0;//lumiSummary->liveFrac();
      lumi_eventsprocessed = 0;
      lumi_eventsfiltered = 0;
    }
}

void NTupleMakerAOD::endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{

}


void NTupleMakerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebugAOD)  cout<<"inside the analyze function"<< endl; 

  track_count = 0;
  goodprimvertex_count = 0;
  primvertex_count = 0;
  muon_count = 0;
  dimuon_count = 0;
  v0_count = 0;
  sv_count = 0;
  gentau_count = 0;
  pfjet_count = 0;
  electron_count = 0;
  photon_count = 0;
  genparticles_count = 0;
  pizero_count = 0;
  errors = 0;

  bool takeevent = true;

  nEvents->Fill(0);
  pv_position = math::XYZPoint(0.,0.,0.);
  
  lumi_eventsprocessed++;
  
  event_nr      = iEvent.id().event();
  event_run      = iEvent.id().run();
  event_timeunix = iEvent.time().unixTime();
  event_timemicrosec = iEvent.time().microsecondOffset();
  event_luminosityblock = iEvent.getLuminosityBlock().luminosityBlock();

  if(cbeamspot)
    {
      edm::Handle<BeamSpot> TheBeamSpot;
      iEvent.getByToken(BeamSpotToken_, TheBeamSpot);
      if(TheBeamSpot.isValid())
	{
	  beamspot_x = TheBeamSpot->x0();
	  beamspot_y = TheBeamSpot->y0();
	  beamspot_z = TheBeamSpot->z0();
	  beamspot_xwidth = TheBeamSpot->BeamWidthX();
	  beamspot_ywidth = TheBeamSpot->BeamWidthY();
	  beamspot_zsigma = TheBeamSpot->sigmaZ();
	  beamspot_cov[0] = TheBeamSpot->covariance(0,0); 
	  beamspot_cov[1] = TheBeamSpot->covariance(0,1); 
	  beamspot_cov[2] = TheBeamSpot->covariance(0,2); 
	  beamspot_cov[3] = TheBeamSpot->covariance(1,1); 
	  beamspot_cov[4] = TheBeamSpot->covariance(1,2); 
	  beamspot_cov[5] = TheBeamSpot->covariance(2,2); 
	  pv_position = math::XYZPoint(TheBeamSpot->x0(), TheBeamSpot->y0(), TheBeamSpot->z0());
	}
      else
	{
	  beamspot_x = 0.;
	  beamspot_y = 0.;
	  beamspot_z = 0.;
	  beamspot_xwidth = 0.;
	  beamspot_ywidth = 0.;
	  beamspot_zsigma = 0.;
	  beamspot_cov[0] = 0.;
	  beamspot_cov[1] = 0.;
	  beamspot_cov[2] = 0.;
	  beamspot_cov[3] = 0.;
	  beamspot_cov[4] = 0.;
	  beamspot_cov[5] = 0.;
	}
    }
  
  if(crecprimvertex)
    {
      edm::Handle<VertexCollection> Vertex;
      iEvent.getByToken(PVToken_, Vertex);
      if(Vertex.isValid()) {
	for(unsigned i = 0 ; i < Vertex->size(); i++) {
	  primvertex_count++;
	  if(i == 0) {
	    primvertex_x = (*Vertex)[i].x();
	    primvertex_y = (*Vertex)[i].y();
	    primvertex_z = (*Vertex)[i].z();
	    primvertex_chi2 = (*Vertex)[i].chi2();
	    primvertex_ndof = (*Vertex)[i].ndof();
	    primvertex_ntracks = (*Vertex)[i].tracksSize();
	    primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
	    primvertex_cov[1] = (*Vertex)[i].covariance(0,1); 
	    primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
	    primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
	    primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
	    primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()
	    Float_t ptq = 0.;
	    for(Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin() ; it != (*Vertex)[i].tracks_end() ; ++it)
	      {
		ptq += (*it)->pt() * (*it)->pt();
	      }
	    primvertex_ptq = ptq;
	    
	    pv_position = (*Vertex)[i].position();
	    primvertex = (*Vertex)[i];
	  }
	  if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake() && (*Vertex)[i].ndof() >= 4 && (*Vertex)[i].z() > -24 && (*Vertex)[i].z() < 24 && (*Vertex)[i].position().Rho() < 2.)
	    goodprimvertex_count++;
	}
      }
    }

  if (crectrack)
    AddTracks(iEvent, iSetup);
  if (crecphoton)
    AddGammas(iEvent, iSetup);
  if (crecv0) {
    AddV0s(iEvent, iSetup, true);
    AddV0s(iEvent, iSetup, false);
    //    std::cout << std::endl;
  }
  if (crecsv) {
    AddSVs(iEvent, iSetup);
  }
  if (crecpizero)
    AddPiZeros(iEvent, iSetup);
  

  if (crecpfjet) 
    {
      if(doDebugAOD)  cout<<"add PF jets"<< endl; 
      int numberOfJets = int(AddPFJets(iEvent,iSetup));
      // if (cSkim) {
      // 	bool goodJets = false;
      // 	for (unsigned int i=0; i<pfjet_count; ++i) {
      // 	  if (pfjet_pt[i]>50. && fabs(pfjet_eta[i])<2.3) {
      // 	    goodJets = true;
      // 	    // std::cout << "Good jet : pt = " << pfjet_pt[i]
      // 	    //  	      << "  eta = " << pfjet_eta[i] 
      // 	    // 	      << "  phi = " << pfjet_phi[i] << std::endl;
      // 	    break;
      // 	  }
      // 	}
      // 	if (!goodJets) return;
      // }
    } // crecpfjet



  genweight = 1.;
  numpileupinteractionsminus = -1;
  numpileupinteractions      = -1;
  numpileupinteractionsplus  = -1;
  numtruepileupinteractions  = -1.0f;
  hepNUP_ = -1;

  genparticles_lheHt = 0.;
  genparticles_noutgoing = 0;

  // generator info and generated particles 
  if(doDebugAOD)  cout<<"add gen info"<< endl; 
  if(cgen && !cdata)
    {
      AddGenHt(iEvent);

      bool haveGenParticles = AddGenParticles(iEvent);

      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByLabel(edm::InputTag("generator"), HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  //	  cout << "Event weight from HEPMC : " << genweight << endl;
	  genid1 = HEPMC->pdf()->id.first;
	  genx1 = HEPMC->pdf()->x.first;
	  genid2 = HEPMC->pdf()->id.second;
	  genx2 = HEPMC->pdf()->x.second;
	  genScale = HEPMC->qScale();
	}

      edm::Handle<vector<PileupSummaryInfo> > PUInfo;
      iEvent.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PUInfo);
      if(PUInfo.isValid())
	{
	  for(vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	    {
	      int BX = PVI->getBunchCrossing();
	      if(BX == -1)
		{ 
		  numpileupinteractionsminus = PVI->getPU_NumInteractions();
		}
	      else if(BX == 0)
		{ 
		  numpileupinteractions = PVI->getPU_NumInteractions();
		}
	      else if(BX == 1)
		{ 
		  numpileupinteractionsplus = PVI->getPU_NumInteractions();
		}
	      
	      numtruepileupinteractions = PVI->getTrueNumInteractions();
	    }
	}
    } // cgen

  tree->Fill();

  //Store Tau embedding information
  // edm::Handle<GenFilterInfo> embeddingWeightHandle;
  // iEvent.getByLabel(edm::InputTag("generator","minVisPtFilter",""), embeddingWeightHandle);
  // embeddingWeight_ = embeddingWeightHandle.isValid() ? embeddingWeightHandle->filterEfficiency() : 1.0;
  //cout << "--- EMBEDDING WEIGHT : " << embeddingWeight_ << endl;
  /*
    if (isRhEmb_){
    iEvent.getByLabel("genZdecayToTausForEmbeddingKineReweight", genDiTauHandle);      
    iEvent.getByLabel(edm::InputTag("TauSpinnerReco","TauSpinnerWT"), TauSpinnerHandle);
    iEvent.getByLabel(edm::InputTag("ZmumuEvtSelEffCorrWeightProducer","weight"), ZmumuEffHandle);
    if(isMC_){
    iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genDiTauMassVsGenDiTauPt"), diTauMassVSdiTauPtHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genTau2EtaVsGenTau1Eta"), tau2EtaVStau1EtaHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightGENembedding","genTau2PtVsGenTau1Pt"), tau2PtVStau1PtHandle);
    }
    else{
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genDiTauMassVsGenDiTauPt"), diTauMassVSdiTauPtHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genTau2EtaVsGenTau1Eta"), tau2EtaVStau1EtaHandle);
      iEvent.getByLabel(edm::InputTag("embeddingKineReweightRECembedding","genTau2PtVsGenTau1Pt"), tau2PtVStau1PtHandle);
    }
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weight"), muonRadiationHandle);
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weightDown"), muonRadiationDownHandle);
    iEvent.getByLabel(edm::InputTag("muonRadiationCorrWeightProducer","weightUp"), muonRadiationUpHandle);

    TauSpinnerWeight = TauSpinnerHandle.isValid() ? (*TauSpinnerHandle) : 1.0;
    ZmumuEffWeight   = ZmumuEffHandle.isValid() ? (*ZmumuEffHandle) : 1.0;
    diTauMassVSdiTauPtWeight = diTauMassVSdiTauPtHandle.isValid() ? (*diTauMassVSdiTauPtHandle) : 1.0;
    tau2EtaVStau1EtaWeight = tau2EtaVStau1EtaHandle.isValid() ? (*tau2EtaVStau1EtaHandle) : 1.0;
    tau2PtVStau1PtWeight = tau2PtVStau1PtHandle.isValid() ? (*tau2PtVStau1PtHandle) : 1.0;
    muonRadiationWeight = muonRadiationHandle.isValid() ? (*muonRadiationHandle) : 1.0;
    muonRadiationDownWeight = muonRadiationDownHandle.isValid() ? (*muonRadiationDownHandle) : 1.0;
    muonRadiationUpWeight = muonRadiationUpHandle.isValid() ? (*muonRadiationUpHandle) : 1.0;
    genDiTauMass_ = genDiTauHandle.isValid() && genDiTauHandle->size()>0 ? genDiTauHandle->at(0).mass() : 9999;
  }
  embeddingWeights_->push_back(TauSpinnerWeight);
  embeddingWeights_->push_back(ZmumuEffWeight);
  embeddingWeights_->push_back(diTauMassVSdiTauPtWeight);
  embeddingWeights_->push_back(tau2EtaVStau1EtaWeight);
  embeddingWeights_->push_back(tau2PtVStau1PtWeight);
  embeddingWeights_->push_back(muonRadiationWeight);
  embeddingWeights_->push_back(muonRadiationDownWeight);
  embeddingWeights_->push_back(muonRadiationUpWeight);
  */
  
 
  
} //void NTupleMakerAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)



Int_t NTupleMakerAOD::HasAnyMother(const GenParticle* particle, int id)
{
  vector<unsigned> bknummother;
  vector<const GenParticle*> bkparticle;
  bknummother.reserve(10);
  bkparticle.reserve(10);
  int level = 0;
  bkparticle.push_back(particle);
  bknummother.push_back(0);
  
  unsigned j = 0;
  while(true)
    {
      if(j == bkparticle[level]->numberOfMothers())
	{
	  level--;
	  if(level == -1){return(0);}
	  j = bknummother[level];
	  bkparticle.resize(level+1);
	  bknummother.resize(level+1);
	  continue;
	}
      
      if(bkparticle[level]->mother(j)->pdgId() == id) return(2);
      if(abs(bkparticle[level]->mother(j)->pdgId()) == abs(id)) return(1);
      
      if(bkparticle[level]->mother(j)->numberOfMothers() > 0)
	{
	  bknummother[level] = j+1;
	  bkparticle.push_back(dynamic_cast<const GenParticle*>(bkparticle[level]->mother(j)));
	  bknummother.push_back(0);
	  j = 0;
	  level++;
	  continue;
	}
      j++;
    }
  return(0);
} // Int_t NTupleMakerAOD::HasAnyMother(const GenParticle* particle, int id)

void NTupleMakerAOD::endJob()
{
}

bool NTupleMakerAOD::AddGenHt(const edm::Event& iEvent) {
  genparticles_lheHt = 0.;
  genparticles_noutgoing = 0;

  edm::Handle<LHEEventProduct> lheEventProduct;  
  iEvent.getByLabel( "externalLHEProducer", lheEventProduct);
  
  if(!lheEventProduct.isValid())
    return false;

  const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

  size_t numParticles = lheParticles.size();
  for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
   int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
   int status = lheEvent.ISTUP[idxParticle];
   if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
       genparticles_lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
       ++genparticles_noutgoing;
   } 
  }

  return true;
}

bool NTupleMakerAOD::AddGenParticles(const edm::Event& iEvent) {

  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByToken(GenParticleCollectionToken_, GenParticles);

  bool passed = false;
  
  if(GenParticles.isValid())
    {
      bool count_partons = false;
      passed = true;
      for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	{
	  bool fill = false;
	  UInt_t info = 0;
	  UInt_t mother = 100;

	  if(abs((*GenParticles)[i].pdgId()) == 13 /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 || 
		 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenMuon : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta() 
	      //			<< "   phi = " << (*GenParticles)[i].phi() 
	      //			<< "   status = " << (*GenParticles)[i].status() 
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 11 /*&& (*GenParticles)[i].pt() > 8.*/ && (*GenParticles)[i].status()==1)
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
                 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenElectron : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta() 
	      //			<< "   phi = " << (*GenParticles)[i].phi() 
	      //			<< "   status = " << (*GenParticles)[i].status() 
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 15 /*&& (*GenParticles)[i].pt() > 10.*/)
	    {
	      fill = false;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother=ZBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother=WBOSON;}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother=TAU;}
	      if(HasAnyMother(&(*GenParticles)[i], 25) > 0 || HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
                 HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3; mother=HIGGS;}
	      //	      std::cout << "GenTau : "  
	      //	       		<< "   pt = " << (*GenParticles)[i].pt() 
	      //	       		<< "   eta = " << (*GenParticles)[i].eta()
	      //	       		<< "   phi = " << (*GenParticles)[i].phi() 
	      //	       		<< "   status = " << (*GenParticles)[i].status() 
	      //	       		<< "   mother = " << mother << std::endl;
	      reco::Candidate::LorentzVector tau_visible_p4 = getVisMomentum(&(*GenParticles)[i]);
	      reco::Candidate::LorentzVector tau_visibleNoLep_p4 = utils_genMatch::getVisMomentumNoLep(&(*GenParticles)[i]);
	      // std::cout << "   visible pt = " << tau_visible_p4.pt() 
	      // 		<< "   eta = " << tau_visible_p4.eta() 
	      // 		<< "   phi = " << tau_visible_p4.phi() 
	      // 		<< "   mode = " << getGenTauDecayMode(&(*GenParticles)[i]) << std::endl;
	      
	      std::string genTauDecayMode = getGenTauDecayMode(&(*GenParticles)[i]);

	      gentau_px[gentau_count] = (*GenParticles)[i].px();
	      gentau_py[gentau_count] = (*GenParticles)[i].py();
	      gentau_pz[gentau_count] = (*GenParticles)[i].pz();
	      gentau_e[gentau_count]  = (*GenParticles)[i].energy();
	      gentau_status[gentau_count] = (*GenParticles)[i].status();

	      gentau_visible_px[gentau_count] = tau_visible_p4.px();
	      gentau_visible_py[gentau_count] = tau_visible_p4.py();
	      gentau_visible_pz[gentau_count] = tau_visible_p4.pz();
	      gentau_visible_e[gentau_count]  = tau_visible_p4.energy();
	      
	      gentau_visible_pt[gentau_count]   = tau_visible_p4.pt();
	      gentau_visible_eta[gentau_count]  = tau_visible_p4.eta();
	      gentau_visible_phi[gentau_count]  = tau_visible_p4.phi();
	      gentau_visible_mass[gentau_count] = tau_visible_p4.mass();

	      gentau_visibleNoLep_px[gentau_count] = tau_visibleNoLep_p4.px();
	      gentau_visibleNoLep_py[gentau_count] = tau_visibleNoLep_p4.py();
	      gentau_visibleNoLep_pz[gentau_count] = tau_visibleNoLep_p4.pz();
	      gentau_visibleNoLep_e[gentau_count]  = tau_visibleNoLep_p4.energy();
	      
	      gentau_visibleNoLep_pt[gentau_count]   = tau_visibleNoLep_p4.pt();
	      gentau_visibleNoLep_eta[gentau_count]  = tau_visibleNoLep_p4.eta();
	      gentau_visibleNoLep_phi[gentau_count]  = tau_visibleNoLep_p4.phi();
	      gentau_visibleNoLep_mass[gentau_count] = tau_visibleNoLep_p4.mass();
	      
	      const GenStatusFlags statusFlags = (*GenParticles)[i].statusFlags();
	      gentau_fromHardProcess[gentau_count] = statusFlags.fromHardProcess();
	      gentau_fromHardProcessBeforeFSR[gentau_count] = statusFlags.fromHardProcessBeforeFSR();
	      gentau_isDecayedLeptonHadron[gentau_count] = statusFlags.isDecayedLeptonHadron();
	      gentau_isDirectHadronDecayProduct[gentau_count] = statusFlags.isDirectHadronDecayProduct();
	      gentau_isDirectHardProcessTauDecayProduct[gentau_count] = statusFlags.isDirectHardProcessTauDecayProduct();
	      gentau_isDirectPromptTauDecayProduct[gentau_count] = statusFlags.isDirectPromptTauDecayProduct();
	      gentau_isDirectTauDecayProduct[gentau_count] = statusFlags.isDirectTauDecayProduct();
	      gentau_isFirstCopy[gentau_count] = statusFlags.isFirstCopy();
	      gentau_isHardProcess[gentau_count] = statusFlags.isHardProcess();
	      gentau_isHardProcessTauDecayProduct[gentau_count] = statusFlags.isHardProcessTauDecayProduct();
	      gentau_isLastCopy[gentau_count] = statusFlags.isLastCopy();
	      gentau_isLastCopyBeforeFSR[gentau_count] = statusFlags.isLastCopyBeforeFSR();
	      gentau_isPrompt[gentau_count] = statusFlags.isPrompt();
	      gentau_isPromptTauDecayProduct[gentau_count] = statusFlags.isPromptTauDecayProduct();
	      gentau_isTauDecayProduct[gentau_count] = statusFlags.isTauDecayProduct();

	      gentau_decayMode_name[gentau_count] = genTauDecayMode;

	      if( genTauDecayMode.find("oneProng0Pi0")!=string::npos )          gentau_decayMode[gentau_count] = 0;
	      else if( genTauDecayMode.find("oneProng1Pi0")!=string::npos )     gentau_decayMode[gentau_count] = 1;
	      else if( genTauDecayMode.find("oneProng2Pi0")!=string::npos )     gentau_decayMode[gentau_count] = 2;
	      else if( genTauDecayMode.find("oneProngOther")!=string::npos )    gentau_decayMode[gentau_count] = 3;
	      else if( genTauDecayMode.find("threeProng0Pi0")!=string::npos )   gentau_decayMode[gentau_count] = 4;
	      else if( genTauDecayMode.find("threeProng1Pi0")!=string::npos )   gentau_decayMode[gentau_count] = 5;
	      else if( genTauDecayMode.find("threeProngOther")!=string::npos )  gentau_decayMode[gentau_count] = 6;
	      else if( genTauDecayMode.find("rare")!=string::npos )             gentau_decayMode[gentau_count] = 7;
	      else if( genTauDecayMode.find("muon")!=string::npos )             gentau_decayMode[gentau_count] = 8;
	      else if( genTauDecayMode.find("electron")!=string::npos )         gentau_decayMode[gentau_count] = 9;
	      else gentau_decayMode[gentau_count] = -99;

	      gentau_mother[gentau_count] = mother;

	      gentau_count++;

	    }
	  else if( (abs((*GenParticles)[i].pdgId()) == 16 || abs((*GenParticles)[i].pdgId()) == 14 || abs((*GenParticles)[i].pdgId()) == 12))
	    {
	      if ((*GenParticles)[i].status()==1) {
		fill = true;
		if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;mother=ZBOSON;}
		if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;mother=WBOSON;}
		if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2;mother=TAU;}
		//		std::cout << "GenNeutrino : " << (*GenParticles)[i].pdgId() 
		//		 	  << "   pt = " << (*GenParticles)[i].pt() 
		//		 	  << "   eta = " << (*GenParticles)[i].eta() 
		//		 	  << "   phi = " << (*GenParticles)[i].phi()
		//			  << "   status = " << (*GenParticles)[i].status() 
		//		 	  << "   mother =  " << mother << std::endl;

	      }
	    }
	  // Save partons (quarks)
	  else if(abs((*GenParticles)[i].pdgId()) < 6)
	    {
	      if ((*GenParticles)[i].status()==3 && count_partons)
		fill = true;
	    }
	  // Save all tops from Madgraph
	  else if ( abs((*GenParticles)[i].pdgId()) == 6 && (*GenParticles)[i].status()==62 ) 
	    {
	      fill = true;
	      //	      	std::cout << "GenTop : " << (*GenParticles)[i].pdgId() 
	      //			  << "   pt = " << (*GenParticles)[i].pt() 
	      //		 	  << "   eta = " << (*GenParticles)[i].eta()
	      //		 	  << "   phi = " << (*GenParticles)[i].phi() 
	      //		 	  << "   status = " << (*GenParticles)[i].status() << std::endl;
	    }
	  // Save partons (gluons) 
	  else if(abs((*GenParticles)[i].pdgId()) == 21)
	    {
	      if ((*GenParticles)[i].status()==3 && count_partons)
		fill = true;
	    }
	  // Save all W/Z bosons from Madgraph
	  else if(abs((*GenParticles)[i].pdgId()) == 23 || abs((*GenParticles)[i].pdgId()) == 24 )
	    {
	      count_partons = true;
	      fill = true;
	      int nDaughters = (*GenParticles)[i].numberOfDaughters();
	      bool posElectronFound = false;
	      bool negElectronFound = false;
	      bool posMuonFound = false;
	      bool negMuonFound = false;
	      bool posTauFound = false;
	      bool negTauFound = false;
	      for (int iD=0 ; iD<nDaughters; ++iD) {
		const reco::Candidate * kid = (*GenParticles)[i].daughter(iD);
		int pdgId = kid->pdgId();
		if (pdgId==11) negElectronFound = true;
		if (pdgId==-11) posElectronFound = true;
		if (pdgId==13) negMuonFound = true;
		if (pdgId==-13) posMuonFound = true;
		if (pdgId==15) negTauFound = true;
		if (pdgId==-15) posTauFound = true;
	      }
	      if (posElectronFound&&negElectronFound) info = 11;
	      else if (posMuonFound&&negMuonFound) info = 13;
	      else if (posTauFound&&negTauFound) info = 15;
	      //	      std::cout << "GenBoson : " << (*GenParticles)[i].pdgId() 
	      //			<< "   pt = " << (*GenParticles)[i].pt() 
	      //			<< "   eta = " << (*GenParticles)[i].eta()
	      //			<< "   phi = " << (*GenParticles)[i].phi() 
	      //			<< "   status = " << (*GenParticles)[i].status() 
	      //			<< "   info = " << info << std::endl;
	    }
	  //Save Higgs bosons
	  else if(abs((*GenParticles)[i].pdgId()) == 25 || abs((*GenParticles)[i].pdgId()) == 35 ||
		  abs((*GenParticles)[i].pdgId()) == 36){
	    fill = true;
	  }
	  if((HasAnyMother(&(*GenParticles)[i], 23) > 0 ||
	      HasAnyMother(&(*GenParticles)[i], 25) > 0 ||
	      HasAnyMother(&(*GenParticles)[i], 35) > 0 ||
	      HasAnyMother(&(*GenParticles)[i], 36) > 0)&&
	     HasAnyMother(&(*GenParticles)[i], 15) > 0 && 
	     (*GenParticles)[i].status()==1 && abs((*GenParticles)[i].pdgId())!=15) {
	    if(HasAnyMother(&(*GenParticles)[i], 23) > 0) {info |= 1<<0;mother=ZBOSON;}
	    if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;mother=WBOSON;}
	    if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2;mother=TAU;}
	    if(HasAnyMother(&(*GenParticles)[i], 25) > 0||
	       HasAnyMother(&(*GenParticles)[i], 35) > 0||
	       HasAnyMother(&(*GenParticles)[i], 36) > 0) {info |= 1<<3;mother=HIGGS;}
	    fill = true;
	  }

	  if(fill)
	    {
	      genparticles_e[genparticles_count] = (*GenParticles)[i].energy();
	      genparticles_px[genparticles_count] = (*GenParticles)[i].px();
	      genparticles_py[genparticles_count] = (*GenParticles)[i].py();
	      genparticles_pz[genparticles_count] = (*GenParticles)[i].pz();
	      genparticles_vx[genparticles_count] = (*GenParticles)[i].vx();
	      genparticles_vy[genparticles_count] = (*GenParticles)[i].vy();
	      genparticles_vz[genparticles_count] = (*GenParticles)[i].vz();
	      genparticles_pdgid[genparticles_count] = (*GenParticles)[i].pdgId();
	      genparticles_status[genparticles_count] = (*GenParticles)[i].status();
	      genparticles_info[genparticles_count] = info;
	      genparticles_mother[genparticles_count] = mother;

	      const GenStatusFlags statusFlags = (*GenParticles)[i].statusFlags();
	      genparticles_fromHardProcess[genparticles_count] = statusFlags.fromHardProcess();
	      genparticles_fromHardProcessBeforeFSR[genparticles_count] = statusFlags.fromHardProcessBeforeFSR();
	      genparticles_isDecayedLeptonHadron[genparticles_count] = statusFlags.isDecayedLeptonHadron();
	      genparticles_isDirectHadronDecayProduct[genparticles_count] = statusFlags.isDirectHadronDecayProduct();
	      genparticles_isDirectHardProcessTauDecayProduct[genparticles_count] = statusFlags.isDirectHardProcessTauDecayProduct();
	      genparticles_isDirectPromptTauDecayProduct[genparticles_count] = statusFlags.isDirectPromptTauDecayProduct();
	      genparticles_isDirectTauDecayProduct[genparticles_count] = statusFlags.isDirectTauDecayProduct();
	      genparticles_isFirstCopy[genparticles_count] = statusFlags.isFirstCopy();
	      genparticles_isHardProcess[genparticles_count] = statusFlags.isHardProcess();
	      genparticles_isHardProcessTauDecayProduct[genparticles_count] = statusFlags.isHardProcessTauDecayProduct();
	      genparticles_isLastCopy[genparticles_count] = statusFlags.isLastCopy();
	      genparticles_isLastCopyBeforeFSR[genparticles_count] = statusFlags.isLastCopyBeforeFSR();
	      genparticles_isPrompt[genparticles_count] = statusFlags.isPrompt();
	      genparticles_isPromptTauDecayProduct[genparticles_count] = statusFlags.isPromptTauDecayProduct();
	      genparticles_isTauDecayProduct[genparticles_count] = statusFlags.isTauDecayProduct();
	      
	      genparticles_count++;
	    }
	} // for(unsigned i = 0 ; i < GenParticles->size() ; i++)
    } // if(GenParticles.isValid())

  return passed;

} // bool NTupleMakerAOD::AddGenParticles(const edm::Event& iEvent) 

unsigned int NTupleMakerAOD::AddV0s(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool isKshort) {

  edm::Handle<reco::VertexCompositeCandidateCollection> Vertices;
  if (isKshort)
    iEvent.getByToken( KshortCollectionToken_, Vertices);
  else
    iEvent.getByToken( LambdaCollectionToken_, Vertices);
  
  

  if(Vertices.isValid())
    {
      //      std::cout << "Here we are" << std::endl;
      for(unsigned i = 0 ; i < Vertices->size() ; i++)
	{
	  v0_px[v0_count] = (*Vertices)[i].px();
	  v0_py[v0_count] = (*Vertices)[i].py();
	  v0_pz[v0_count] = (*Vertices)[i].pz();
	  v0_vx[v0_count] = (*Vertices)[i].vx();
	  v0_vy[v0_count] = (*Vertices)[i].vy();
	  v0_vz[v0_count] = (*Vertices)[i].vz();
	  v0_ID[v0_count] = (*Vertices)[i].pdgId();
	  v0_mass[v0_count] = (*Vertices)[i].mass();
	  v0_chi2[v0_count] = (*Vertices)[i].vertexChi2();
	  v0_ndof[v0_count] = (*Vertices)[i].vertexNdof();

	  TLorentzVector v0LV; v0LV.SetXYZM(v0_px[v0_count],v0_py[v0_count],v0_pz[v0_count],
					      v0_mass[v0_count]);

	  v0_pt[v0_count]  = v0LV.Pt();
	  v0_eta[v0_count] = v0LV.Eta();
	  v0_phi[v0_count] = v0LV.Phi();

	  float deltax = v0_vx[v0_count] - primvertex_x;
	  float deltay = v0_vy[v0_count] - primvertex_y;
	  float deltaz = v0_vz[v0_count] - primvertex_z;

	  v0_decay[v0_count] = TMath::Sqrt(deltax*deltax+deltay*deltay+deltaz*deltaz);

	  float vp = deltax*v0_px[v0_count]+deltay*v0_py[v0_count]+deltaz*v0_pz[v0_count];
	  float p2 = TMath::Sqrt(v0_px[v0_count]*v0_px[v0_count]+v0_py[v0_count]*v0_py[v0_count]+v0_pz[v0_count]*v0_pz[v0_count]);
	  float time = -vp/p2;
	  float dcax = v0_vx[v0_count]+time*v0_px[v0_count];
	  float dcay = v0_vy[v0_count]+time*v0_py[v0_count];
	  float dcaz = v0_vz[v0_count]+time*v0_pz[v0_count];

	  v0_ip[v0_count] = TMath::Sqrt(dcax*dcax+dcay*dcay+dcaz*dcaz);

	  unsigned int numDaughters = (*Vertices)[i].numberOfDaughters();
	  v0_ndaughters[v0_count] = numDaughters;
	  if (numDaughters==2) {
	    for (unsigned int iD=0; iD<2; ++iD) {
	      const reco::Candidate * daughter = (*Vertices)[i].daughter(iD);
	      if (daughter->charge()>0) {
		v0_pos_px[v0_count] = daughter->px();
		v0_pos_py[v0_count] = daughter->py();
		v0_pos_pz[v0_count] = daughter->pz();
		v0_pos_mass[v0_count] = daughter->mass();
		v0_pos_ID[v0_count] = daughter->pdgId();
		TLorentzVector partLV; partLV.SetXYZM(v0_pos_px[v0_count],v0_pos_py[v0_count],v0_pos_pz[v0_count],
						      v0_pos_mass[v0_count]);
		v0_pos_pt[v0_count]  = partLV.Pt();
		v0_pos_eta[v0_count] = partLV.Eta();
		v0_pos_phi[v0_count] = partLV.Phi();
	      }
	      else {
		v0_neg_px[v0_count] = daughter->px();
                v0_neg_py[v0_count] = daughter->py();
                v0_neg_pz[v0_count] = daughter->pz();
                v0_neg_mass[v0_count] = daughter->mass();
                v0_neg_ID[v0_count] = daughter->pdgId();
		TLorentzVector partLV; partLV.SetXYZM(v0_neg_px[v0_count],v0_neg_py[v0_count],v0_neg_pz[v0_count],
						      v0_neg_mass[v0_count]);
		v0_neg_pt[v0_count]  = partLV.Pt();
		v0_neg_eta[v0_count] = partLV.Eta();
		v0_neg_phi[v0_count] = partLV.Phi();
	      }

	    }
	  }

	  /*
	  std::cout << "Vertex " << v0_count 
		    << "  pdgId = " << v0_ID[v0_count]
		    << "  mass = " << v0_mass[v0_count] 
		    << "  ip = " << v0_ip[v0_count] 
		    << "  dL = " << v0_decay[v0_count]
		    << "  number of daughters = " << numDaughters
		    << std::endl;
	  std::cout << "   neg  pdgId = " << v0_neg_ID[v0_count]
		    << "  mass = " << v0_neg_mass[v0_count]
		    << "  vx = " << v0_neg_vx[v0_count]
		    << "  vy = " << v0_neg_vy[v0_count]
		    << "  vz = " << v0_neg_vz[v0_count] 
		    << "  px = " << v0_neg_px[v0_count]
		    << "  py = " << v0_neg_py[v0_count]
		    << "  pz = " << v0_neg_pz[v0_count] 
		    << std::endl;
	  std::cout << "   pos  pdgId = " << v0_pos_ID[v0_count]
		    << "  mass = " << v0_pos_mass[v0_count]
		    << "  vx = " << v0_pos_vx[v0_count]
		    << "  vy = " << v0_pos_vy[v0_count]
		    << "  vz = " << v0_pos_vz[v0_count] 
		    << "  px = " << v0_pos_px[v0_count]
		    << "  py = " << v0_pos_py[v0_count]
		    << "  pz = " << v0_pos_pz[v0_count] 
		    << std::endl;
	  */
	  v0_count++;

	  if (v0_count>=M_v0maxcount) {
	    cerr << "number of V0s > M_v0maxcount. They are missing." << endl; errors |= 1<<1;
	    break;
	  }

	}
    }      

  return v0_count;

}

unsigned int NTupleMakerAOD::AddSVs(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::VertexCollection> SecVertices;
  iEvent.getByToken( SVToken_, SecVertices);

  if (SecVertices.isValid()) {
    for(unsigned i = 0 ; i < SecVertices->size(); i++) {
      sv_vx[sv_count] = (*SecVertices)[i].x();
      sv_vy[sv_count] = (*SecVertices)[i].y();
      sv_vz[sv_count] = (*SecVertices)[i].z();
      sv_chi2[sv_count] = (*SecVertices)[i].chi2();
      sv_ndof[sv_count] = (*SecVertices)[i].ndof();
      sv_covxx[sv_count] = (*SecVertices)[i].covariance(0,0);
      sv_covxy[sv_count] = (*SecVertices)[i].covariance(0,1);
      sv_covxz[sv_count] = (*SecVertices)[i].covariance(0,2);
      sv_covyy[sv_count] = (*SecVertices)[i].covariance(1,1);
      sv_covyz[sv_count] = (*SecVertices)[i].covariance(1,2);
      sv_covzz[sv_count] = (*SecVertices)[i].covariance(2,2);
      sv_ntracks[sv_count] = (*SecVertices)[i].tracksSize();
      unsigned int nTrk = 0;
      bool exceedVectorSize = false;
      for(Vertex::trackRef_iterator it = (*SecVertices)[i].tracks_begin() ; it != (*SecVertices)[i].tracks_end() ; ++it)
	{
	  sv_track_px[sv_count][nTrk] = (*it)->px();
	  sv_track_py[sv_count][nTrk] = (*it)->py();
	  sv_track_pz[sv_count][nTrk] = (*it)->pz();
	  sv_track_pt[sv_count][nTrk] = (*it)->pt();
	  sv_track_eta[sv_count][nTrk] = (*it)->eta();
	  sv_track_phi[sv_count][nTrk] = (*it)->phi();
	  sv_track_charge[sv_count][nTrk] = (*it)->charge();
	  nTrk++;
	  if (nTrk>=20) { 
	    break;
	    exceedVectorSize = true;
	  }
	}
      sv_exceedTrkSize[sv_count] = exceedVectorSize;

    }
  }

  return sv_count;

}

unsigned int NTupleMakerAOD::AddPiZeros(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::RecoTauPiZeroCollection> Strips;
  iEvent.getByToken( TauPiZeroCollectionToken_, Strips);

  if (Strips.isValid()) {
    if (Strips->size()>0)
      std::cout << "Number of Pizeros = " << Strips->size() << std::endl;
    for (unsigned int i=0; i<Strips->size(); ++i) {
      pizero_px[pizero_count] = (*Strips)[i].px();
      pizero_py[pizero_count] = (*Strips)[i].py();
      pizero_pz[pizero_count] = (*Strips)[i].pz();
      pizero_e[pizero_count] = (*Strips)[i].p();
      pizero_pt[pizero_count] = (*Strips)[i].pt();
      pizero_eta[pizero_count] = (*Strips)[i].eta();
      pizero_phi[pizero_count] = (*Strips)[i].phi();
      pizero_x[pizero_count] = (*Strips)[i].vx();
      pizero_y[pizero_count] = (*Strips)[i].vy();
      pizero_z[pizero_count] = (*Strips)[i].vz();
      pizero_count++;
      std::cout << " " << i 
		<< "  px = " << pizero_px[pizero_count]
		<< "  py = " << pizero_py[pizero_count]
		<< "  pz = " << pizero_pz[pizero_count]
		<< "  vx = " << pizero_x[pizero_count]
		<< "  vy = " << pizero_y[pizero_count]
		<< "  vz = " << pizero_z[pizero_count] << std::endl;

      if (pizero_count>=M_pizeromaxcount) {
	cerr << "number of pizeros > M_pizeromaxcount. They are missing." << endl; errors |= 1<<1;
	break;
      }
    }
  } 

  return pizero_count;

}

unsigned int NTupleMakerAOD::AddTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::PFCandidateCollection> Tracks;
  iEvent.getByToken( PFCandidateCollectionToken_, Tracks);

  if(Tracks.isValid())
    {
      for(unsigned i = 0 ; i < Tracks->size() ; i++){
	const reco::Track * trkRef = (*Tracks)[i].bestTrack();
	if (trkRef==NULL) continue;
	if ((*Tracks)[i].pt() < cTrackPtMin) continue;
        if (fabs((*Tracks)[i].eta()) > cTrackEtaMax) continue;
        if (fabs((*Tracks)[i].charge()) < 0.5) continue;
	if (fabs(trkRef->dxy()) > cTrackDxyMax) continue;
	if (fabs(trkRef->dz()) > cTrackDzMax) continue;
        track_px[track_count] = (*Tracks)[i].px();
        track_py[track_count] = (*Tracks)[i].py();
        track_pz[track_count] = (*Tracks)[i].pz();
        track_pt[track_count] = (*Tracks)[i].pt();
        track_eta[track_count] = (*Tracks)[i].eta();
        track_phi[track_count] = (*Tracks)[i].phi();
        track_charge[track_count] = (*Tracks)[i].charge();
        track_mass[track_count] = (*Tracks)[i].mass();
        track_dxy[track_count] = trkRef->dxy();
        track_dz[track_count] = trkRef->dz();
        track_ID[track_count] = (*Tracks)[i].pdgId();
	track_highPurity[track_count] = trkRef->quality(reco::Track::highPurity);

        track_count++;

        if (track_count==M_trackmaxcount) {
          cerr << "number of tracks > M_trackmaxcount. They are missing." << endl; errors |= 1<<1;
          break;
        }
      }

    }

  return track_count;

}

unsigned int NTupleMakerAOD::AddGammas(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::PFCandidateCollection> Tracks;
  iEvent.getByToken( PFCandidateCollectionToken_, Tracks);

  if(Tracks.isValid())
    {
      for(unsigned i = 0 ; i < Tracks->size() ; i++){
	if ((*Tracks)[i].pdgId()!=22) continue;
        if ((*Tracks)[i].pt() < cPhotonPtMin) continue;
        if (fabs((*Tracks)[i].eta()) > cPhotonEtaMax) continue;
        photon_px[photon_count] = (*Tracks)[i].px();
        photon_py[photon_count] = (*Tracks)[i].py();
        photon_pz[photon_count] = (*Tracks)[i].pz();
        photon_pt[photon_count] = (*Tracks)[i].pt();
        photon_eta[photon_count] = (*Tracks)[i].eta();
        photon_phi[photon_count] = (*Tracks)[i].phi();
	photon_e[photon_count] = (*Tracks)[i].p();
	const reco::SuperClusterRef superCluster = (*Tracks)[i].superClusterRef();
	photon_superClusterX[photon_count] = superCluster->x();
	photon_superClusterY[photon_count] = superCluster->y();
	photon_superClusterZ[photon_count] = superCluster->z();
        photon_count++;

        if (photon_count==M_photonmaxcount) {
          cerr << "number of tracks > M_trackmaxcount. They are missing." << endl; errors |= 1<<1;
          break;
        }
      }

    }

  return photon_count;

}

// #if 0
// template<typename TCollection>
// const NSVfitEventHypothesisBase* matchSVfitHypothesis(const NSVfitEventHypothesisBase& hypo, const TCollection& coll) //const NSVfitEventHypothesisBaseCollection& coll)
// {
//   for(unsigned int i = 0; i < coll.size(); ++i)
//     {
//       const NSVfitEventHypothesisBase& cand = coll[i];
//       assert(hypo.numResonances() == cand.numResonances());
      
//       bool particle_match = true;
//       for(unsigned int j = 0; j < hypo.numResonances(); ++j)
// 	{
// 	  const NSVfitResonanceHypothesisBase* hypo_res = hypo.resonance(j);
// 	  const NSVfitResonanceHypothesisBase* cand_res = cand.resonance(j);
// 	  assert(hypo_res && cand_res);
	  
// 	  assert(hypo_res->numDaughters() == cand_res->numDaughters());
// 	  for(unsigned int k = 0; k < hypo_res->numDaughters(); ++k)
// 	    {
// 	      if(hypo_res->daughter(k)->particle() != cand_res->daughter(k)->particle())
// 		particle_match = false;
// 	    }
// 	}
      
//       if(particle_match) return &cand;
//     }
  
//   assert(false);
//   return NULL;
// }
// #endif

unsigned int NTupleMakerAOD::AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::PFJetCollection> pfjets;
  iEvent.getByToken(JetCollectionToken_, pfjets);
  
  //  edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > svInfos;
  //  iEvent.getByLabel(edm::InputTag("secondaryVertexTagInfosEI"), svInfos);
  //  assert(svInfos.isValid());
  
  //  edm::Handle<edm::ValueMap<float> > puJetIdMVAFull;
  //iEvent.getByLabel(edm::InputTag("puJetIdForPFMVAMEt","fullDiscriminant"), puJetIdMVAFull);
  //  iEvent.getByLabel(edm::InputTag("pileupJetIdFull","full53xDiscriminant"), puJetIdMVAFull);

  //  edm::Handle<reco::PFJetCollection> ak4jets;
  //iEvent.getByLabel(edm::InputTag("calibratedAK4PFJetsForPFMVAMEt"), ak4jets);
  //iEvent.getByLabel(edm::InputTag("ak4PFJets"), ak4jets);
  //iEvent.getByLabel(edm::InputTag("AK4PFCHS"), ak4jets);
  //  iEvent.getByLabel(edm::InputTag("slimmedJets"), ak4jets);
  
  //	edm::Handle<edm::ValueMap<int> > puJetIdFlagFull;
  //	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullId"), puJetIdFlagFull);

  if(pfjets.isValid())
    {
      // if (pfjets->size()>0) {
      // 	const std::vector< std::pair< std::string, float > > pairDiscriVector = (*pfjets)[0].getPairDiscri();
      // 	int nDiscri = pairDiscriVector.size();
      // 	std::cout << "Number of discriminators = " << nDiscri << std::endl;
      // 	for (int iD=0;iD<nDiscri;++iD) {
      // 	  std::pair<std::string, float> pairDiscri = pairDiscriVector[iD];
      // 	  std::cout << "Dicsr = " << pairDiscriVector[iD].first << std::endl;
      // 	}
      // }
      // std::cout << std::endl;
      
      
      for(unsigned i = 0 ; i < pfjets->size() ; i++)
	{

	  if(pfjet_count == M_jetmaxcount){
	    cerr << "number of pfjet > M_jetmaxcount. They are missing." << endl; 
	    errors |= 1<<4; 
	    break;
	  }

	  if((*pfjets)[i].pt() < cJetPtMin) continue;
	  if(fabs((*pfjets)[i].eta()) > cJetEtaMax) continue;
	  //	  std::cout << "Jet  " << i <<  ", pT=" <<  (*pfjets)[i].pt() << std::endl;
	  
	  pfjet_e[pfjet_count] = (*pfjets)[i].energy();
	  pfjet_px[pfjet_count] = (*pfjets)[i].px();
	  pfjet_py[pfjet_count] = (*pfjets)[i].py();
	  pfjet_pz[pfjet_count] = (*pfjets)[i].pz();
          pfjet_pt[pfjet_count] = (*pfjets)[i].pt();
          pfjet_eta[pfjet_count] = (*pfjets)[i].eta();
          pfjet_phi[pfjet_count] = (*pfjets)[i].phi();
	  pfjet_neutralhadronicenergy[pfjet_count] = (*pfjets)[i].neutralHadronEnergy();
	  pfjet_chargedhadronicenergy[pfjet_count] = (*pfjets)[i].chargedHadronEnergy();
	  pfjet_neutralemenergy[pfjet_count] = (*pfjets)[i].neutralEmEnergy();
	  pfjet_chargedemenergy[pfjet_count] = (*pfjets)[i].chargedEmEnergy();
	  pfjet_muonenergy[pfjet_count] = (*pfjets)[i].muonEnergy();
	  pfjet_chargedmuonenergy[pfjet_count] = (*pfjets)[i].chargedMuEnergy();
	  pfjet_chargedmulti[pfjet_count] = (*pfjets)[i].chargedMultiplicity();
	  pfjet_neutralmulti[pfjet_count] = (*pfjets)[i].neutralMultiplicity();
	  pfjet_chargedhadronmulti[pfjet_count] = (*pfjets)[i].chargedHadronMultiplicity();

	  pfjet_energycorr[pfjet_count] = -1.;
	  pfjet_energycorr_l1fastjet[pfjet_count] = -1.;
	  pfjet_energycorr_l2relative[pfjet_count] = -1.;
	  pfjet_energycorr_l3absolute[pfjet_count] = -1.;
	  pfjet_energycorr_l2l3residual[pfjet_count] = -1.;
		
	  // std::cout << "Jet Energy corrections : " << std::endl;
	  // std::cout << "    L1FastJet    = " << pfjet_energycorr_l1fastjet[pfjet_count] << std::endl;
	  // std::cout << "    L2Relative   = " << pfjet_energycorr_l2relative[pfjet_count] << std::endl;
	  // std::cout << "    L3Absolute   = " << pfjet_energycorr_l3absolute[pfjet_count] << std::endl;
	  // std::cout << "    L2L3Residual = " << pfjet_energycorr_l2l3residual[pfjet_count] << std::endl;
	  // std::cout << "    Total (Uncor)= " << pfjet_energycorr[pfjet_count] << std::endl;
	  
	  
	  //pfjet_pu_jet_simple_loose[pfjet_count] = false;
	  //pfjet_pu_jet_simple_medium[pfjet_count] = false;
	  //pfjet_pu_jet_simple_tight[pfjet_count] = false;
	  //pfjet_pu_jet_simple_mva[pfjet_count] = (*puJetIdMVAFull)[(*pfjets)[i].originalObjectRef()];
	  
	  
	  //		    pfjet_pu_jet_full_loose[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kLoose);
	  //		    pfjet_pu_jet_full_medium[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kMedium);
	  //		    pfjet_pu_jet_full_tight[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kTight);
	  //                pfjet_pu_jet_full_mva[pfjet_count] = (*puJetIdMVAFull)[(*pfjets)[i].originalObjectRef()];
	  
	  //get MVA Id
          //for(reco::PFJetCollection::const_iterator iak4jets = ak4jets->begin(); iak4jets != ak4jets->end(); iak4jets++){
	  //	  pfjet_pu_jet_full_mva[pfjet_count] = -9999;
	  //	  if (puJetIdMVAFull.isValid()&&ak4jets.isValid()) {
	  //	    for(size_t ij = 0; ij < ak4jets->size(); ij++){
	  //	      reco::PFJetRef jetRef (ak4jets, ij);
	  //	      if(deltaR((*pfjets)[i].p4(), jetRef->p4()) < 0.3){
	  //		if(deltaR((*pfjets)[i].p4(), jetRef->p4()) > 0.1)
	  //		  std::cout<<"original jet pt "<<(*pfjets)[i].pt()<<" re-recoed jet pt "<<jetRef->pt()<<" pu mva value "<<(*puJetIdMVAFull)[jetRef]<<std::endl;
	  //		pfjet_pu_jet_full_mva[pfjet_count] = (*puJetIdMVAFull)[jetRef];
	  //	      }
	  //	    }
	  //	  }	  
	  //	  pfjet_pu_jet_full_mva[pfjet_count] = (*pfjets)[i].userFloat("pileupJetId:fullDiscriminant");
	  jecUnc->setJetEta(pfjet_eta[pfjet_count]);
	  jecUnc->setJetPt(pfjet_pt[pfjet_count]);
	  pfjet_jecUncertainty[pfjet_count] = jecUnc->getUncertainty(true);
		
	  //	  for(unsigned n = 0 ; n < cBtagDiscriminators.size() ; n++)
	  //	    {
	  //	      pfjet_btag[pfjet_count][n] = -1000;
	  //	      if(cBtagDiscriminators[n] != "F"){
		//		std::cout << " " << cBtagDiscriminators.at(n) << "  : " <<  (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) << std::endl;
	  //		pfjet_btag[pfjet_count][n] = (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) ;
	  //	      }
	  //	    }
	  pfjet_count++;
	}
    }
  
  return  pfjet_count;
}


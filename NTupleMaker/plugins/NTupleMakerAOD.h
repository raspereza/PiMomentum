#ifndef NTupleTauIdEfficiencyW_h
#define NTupleTauIdEfficiencyW_h
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <memory>

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
//#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h" 
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
//#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;
using namespace reco;

#define M_trackmaxcount 1000
#define M_superclustermaxcount 1000
#define M_superclustermembermaxcount 1000
#define M_superclusterhitmaxcount 5000
#define M_primvertexmaxcount 1000
#define M_muonmaxcount 1000
#define M_taumaxcount 1000
#define M_electronmaxcount 1000
#define M_photonmaxcount 1000
#define M_conversionmaxcount 1000
#define M_jetmaxcount 1000
#define M_mvametmaxcount 2000
#define M_genparticlesmaxcount 1000
#define M_trigobjectmaxcount 1000
#define M_pzeromaxcount 1000
#define M_v0maxcount 1000
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::SMatrix<double, 2, 2, ROOT::Math::MatRepSym<double, 2> > CovMatrix2D;

bool doDebugAOD = false;
class NTupleMakerAOD : public edm::EDAnalyzer{ 
 public:
  explicit NTupleMakerAOD( const edm::ParameterSet& iConfig );
  ~NTupleMakerAOD();

 private:
  enum MotherNames{HIGGS=1, WBOSON, ZBOSON, TAU};
  enum MvaMetChannel{EMU=1, ETAU, MUTAU, TAUTAU, MUMU, EE, UNKNOWN};

  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup );

  bool AddGenHt(const edm::Event& iEvent);
  bool AddGenParticles(const edm::Event& iEvent);
  //  unsigned int AddTaus(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddGammas(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddV0s(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool isKshort);
  //  unsigned int AddPiZeros(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  Int_t HasAnyMother(const GenParticle* particle, int id);
  math::XYZPoint PositionOnECalSurface(reco::TransientTrack&);

  //  DCA calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2);

  TTree* tree;
  TH1D*  nEvents;

  //Configuration (steering cards)

  bool cdata;
  unsigned int cYear;
  std::string cPeriod;
  unsigned int cSkim;
  std::string cJECfile;
  
  bool cgen;
  bool cbeamspot;
  bool crecprimvertex;
  bool crectrack;
  bool crecphoton;
  bool crectau;
  bool crecpizero;
  bool crecpfjet;
  bool crecv0;

  double cMuPtMin;
  double cMuEtaMax;
  int cMuNum;

  // electrons
  double cElPtMin;
  double cElEtaMax;
  int cElNum;

  // taus
  double cTauPtMin;
  double cTauEtaMax;
  int cTauNum;

  // tracks
  double cTrackPtMin;
  double cTrackEtaMax;
  double cTrackDxyMax;
  double cTrackDzMax;
  int cTrackNum;

  // photons
  double cPhotonPtMin;
  double cPhotonEtaMax;
  int cPhotonNum;

  // pizeros
  double cPizeroPtMin;
  double cPizeroEtaMax;
  int cPizeroNum;

  // jets
  double cJetPtMin;
  double cJetEtaMax;
  int cJetNum;

  std::string sampleName;

  // Collections
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> KshortCollectionToken_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> LambdaCollectionToken_;
  edm::EDGetTokenT<reco::PFTauCollection> TauCollectionToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> PFCandidateCollectionToken_;
  edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;
  edm::EDGetTokenT<reco::PFJetCollection> JetCollectionToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> GenParticleCollectionToken_;
  edm::EDGetTokenT<BeamSpot> BeamSpotToken_;
  edm::EDGetTokenT<VertexCollection> PVToken_;


  //Variables
  edm::ESHandle<TransientTrackBuilder>  TTrackBuilder         ;
  edm::ESHandle<MagneticField>          magneticField         ; 
  Cylinder::ConstCylinderPointer        ecalBarrel            ;
  Plane::ConstPlanePointer              ecalNegativeEtaEndcap ;
  Plane::ConstPlanePointer              ecalPositiveEtaEndcap ;
  //  MEtRecoilCorrection *                 metRecCorr; 
  //  RecoilCorrector *                     corrector_ ;  
  
  math::XYZPoint pv_position;
  Vertex primvertex;

  //Data		
  Int_t njets4RC;

  UInt_t errors;
  UInt_t event_nr;
  UInt_t event_luminosityblock;
  UInt_t event_run;
  UInt_t event_timeunix;
  UInt_t event_timemicrosec;
  UChar_t trigger_level1bits[8];
  UChar_t trigger_level1[128];
  UChar_t trigger_HLT[128];

  // beam spot   
  Float_t beamspot_x;
  Float_t beamspot_y;
  Float_t beamspot_z;
  Float_t beamspot_xwidth;
  Float_t beamspot_ywidth;
  Float_t beamspot_zsigma;
  Float_t beamspot_cov[6];

  // primary vertex
  UInt_t  primvertex_count;
  UInt_t  goodprimvertex_count;
  Float_t primvertex_x;
  Float_t primvertex_y;
  Float_t primvertex_z;
  Float_t primvertex_chi2;
  Float_t primvertex_ndof;
  Float_t primvertex_ptq;
  Int_t   primvertex_ntracks;
  Float_t primvertex_cov[6];

  // tracks
  //tracks                                                                                                                                  
  UInt_t track_count;
  Float_t track_px[M_trackmaxcount];
  Float_t track_py[M_trackmaxcount];
  Float_t track_pz[M_trackmaxcount];
  Float_t track_pt[M_trackmaxcount];
  Float_t track_eta[M_trackmaxcount];
  Float_t track_phi[M_trackmaxcount];
  Float_t track_mass[M_trackmaxcount];
  Float_t track_charge[M_trackmaxcount];
  Float_t track_dxy[M_trackmaxcount];
  Float_t track_dxyerr[M_trackmaxcount];
  Float_t track_dz[M_trackmaxcount];
  Float_t track_dzerr[M_trackmaxcount];
  Int_t track_ID[M_trackmaxcount];
  Bool_t track_highPurity[M_trackmaxcount];

  // muons
  UInt_t muon_count;
  Float_t muon_px[M_muonmaxcount];
  Float_t muon_py[M_muonmaxcount];
  Float_t muon_pz[M_muonmaxcount];
  Float_t muon_pt[M_muonmaxcount];
  Float_t muon_eta[M_muonmaxcount];
  Float_t muon_phi[M_muonmaxcount];
  Float_t muon_pterror[M_muonmaxcount];
  Float_t muon_chi2[M_muonmaxcount];
  Float_t muon_normChi2[M_muonmaxcount];
  Float_t muon_ndof[M_muonmaxcount];
  Float_t muon_charge[M_muonmaxcount];

  // needed for medium muon Id
  Float_t muon_combQ_chi2LocalPosition[M_muonmaxcount];
  Float_t muon_combQ_trkKink[M_muonmaxcount];
  Float_t muon_validFraction[M_muonmaxcount];
  Float_t muon_segmentComp[M_muonmaxcount];

  UInt_t muon_nMuonStations[M_muonmaxcount];
  UInt_t muon_nMuonHits[M_muonmaxcount];
  UInt_t muon_nPixelHits[M_muonmaxcount];
  UInt_t muon_nTrackerHits[M_muonmaxcount];

  Float_t muon_dz[M_muonmaxcount];
  Float_t muon_dzerr[M_muonmaxcount];
  Float_t muon_dxy[M_muonmaxcount];
  Float_t muon_dxyerr[M_muonmaxcount];

  Float_t muon_chargedHadIso[M_muonmaxcount];
  Float_t muon_neutralHadIso[M_muonmaxcount];
  Float_t muon_photonIso[M_muonmaxcount];
  Float_t muon_puIso[M_muonmaxcount];

  Float_t muon_r03_sumChargedHadronPt[M_muonmaxcount];
  Float_t muon_r03_sumChargedParticlePt[M_muonmaxcount];
  Float_t muon_r03_sumNeutralHadronEt[M_muonmaxcount];
  Float_t muon_r03_sumPhotonEt[M_muonmaxcount];
  Float_t muon_r03_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
  Float_t muon_r03_sumPhotonEtHighThreshold[M_muonmaxcount];
  Float_t muon_r03_sumPUPt[M_muonmaxcount];

  Float_t muon_r04_sumChargedHadronPt[M_muonmaxcount];
  Float_t muon_r04_sumChargedParticlePt[M_muonmaxcount];
  Float_t muon_r04_sumNeutralHadronEt[M_muonmaxcount];
  Float_t muon_r04_sumPhotonEt[M_muonmaxcount];
  Float_t muon_r04_sumNeutralHadronEtHighThreshold[M_muonmaxcount];
  Float_t muon_r04_sumPhotonEtHighThreshold[M_muonmaxcount];
  Float_t muon_r04_sumPUPt[M_muonmaxcount];

  Bool_t muon_isPF[M_muonmaxcount];
  Bool_t muon_isGlobal[M_muonmaxcount];
  Bool_t muon_isTracker[M_muonmaxcount];
  Bool_t muon_isTight[M_muonmaxcount];
  Bool_t muon_isLoose[M_muonmaxcount];
  Bool_t muon_isMedium[M_muonmaxcount];

  Bool_t muon_globalTrack[M_muonmaxcount];
  Bool_t muon_innerTrack[M_muonmaxcount];
  
  Int_t muon_genmatch[M_muonmaxcount];

  UInt_t dimuon_count;
  UInt_t dimuon_leading[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  UInt_t dimuon_trailing[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist2D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist2DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist3D[M_muonmaxcount*(M_muonmaxcount - 1)/2];
  Float_t dimuon_dist3DE[M_muonmaxcount*(M_muonmaxcount - 1)/2];

  // pat jets 
  UInt_t pfjet_count;
  Float_t pfjet_e[M_jetmaxcount];
  Float_t pfjet_px[M_jetmaxcount];
  Float_t pfjet_py[M_jetmaxcount];
  Float_t pfjet_pz[M_jetmaxcount];
  Float_t pfjet_pt[M_jetmaxcount];
  Float_t pfjet_eta[M_jetmaxcount];
  Float_t pfjet_phi[M_jetmaxcount];

  Float_t pfjet_neutralhadronicenergy[M_jetmaxcount];
  Float_t pfjet_chargedhadronicenergy[M_jetmaxcount];
  Float_t pfjet_neutralemenergy[M_jetmaxcount];
  Float_t pfjet_chargedemenergy[M_jetmaxcount];
  Float_t pfjet_muonenergy[M_jetmaxcount];
  Float_t pfjet_chargedmuonenergy[M_jetmaxcount];

  UInt_t pfjet_chargedmulti[M_jetmaxcount];	
  UInt_t pfjet_neutralmulti[M_jetmaxcount];	
  UInt_t pfjet_chargedhadronmulti[M_jetmaxcount];
  Float_t pfjet_energycorr[M_jetmaxcount];
  Float_t pfjet_energycorr_l1fastjet[M_jetmaxcount];
  Float_t pfjet_energycorr_l2relative[M_jetmaxcount];
  Float_t pfjet_energycorr_l3absolute[M_jetmaxcount];
  Float_t pfjet_energycorr_l2l3residual[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_cut_mva[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_simple_mva[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_full_mva[M_jetmaxcount];
  Int_t pfjet_flavour[M_jetmaxcount];
  Float_t pfjet_btag[M_jetmaxcount][10];
  Float_t pfjet_jecUncertainty[M_jetmaxcount];

  // pat electrons 
  UInt_t electron_count;
  Float_t electron_px[M_electronmaxcount];
  Float_t electron_py[M_electronmaxcount];
  Float_t electron_pz[M_electronmaxcount];
  Float_t electron_pt[M_electronmaxcount];
  Float_t electron_eta[M_electronmaxcount];
  Float_t electron_phi[M_electronmaxcount];
  Float_t electron_trackchi2[M_electronmaxcount];
  Float_t electron_trackndof[M_electronmaxcount];
  Float_t electron_outerx[M_electronmaxcount];
  Float_t electron_outery[M_electronmaxcount];
  Float_t electron_outerz[M_electronmaxcount];
  Float_t electron_closestpointx[M_electronmaxcount];
  Float_t electron_closestpointy[M_electronmaxcount];
  Float_t electron_closestpointz[M_electronmaxcount];
  Float_t electron_esuperclusterovertrack[M_electronmaxcount];
  Float_t electron_eseedclusterovertrack[M_electronmaxcount];
  Float_t electron_deltaetasuperclustertrack[M_electronmaxcount];
  Float_t electron_deltaphisuperclustertrack[M_electronmaxcount];
  Float_t electron_e1x5[M_electronmaxcount];
  Float_t electron_e2x5[M_electronmaxcount];
  Float_t electron_e5x5[M_electronmaxcount];
  Float_t electron_sigmaetaeta[M_electronmaxcount];
  Float_t electron_sigmaietaieta[M_electronmaxcount];
  Float_t electron_ehcaloverecal[M_electronmaxcount];
  Float_t electron_ehcaloverecaldepth1[M_electronmaxcount];
  Float_t electron_ehcaloverecaldepth2[M_electronmaxcount];
  Float_t electron_full5x5_sigmaietaieta[M_electronmaxcount];
  Float_t electron_ooemoop[M_electronmaxcount];

  Float_t electron_superClusterEta[M_electronmaxcount];
  Float_t electron_superClusterPhi[M_electronmaxcount];
  
  Float_t electron_superClusterX[M_electronmaxcount];
  Float_t electron_superClusterY[M_electronmaxcount];
  Float_t electron_superClusterZ[M_electronmaxcount];

  Float_t electron_chargedHadIso[M_electronmaxcount];
  Float_t electron_neutralHadIso[M_electronmaxcount];
  Float_t electron_photonIso[M_electronmaxcount];
  Float_t electron_puIso[M_electronmaxcount];
  Float_t electron_miniISO[M_electronmaxcount];

  Float_t electron_r03_sumChargedHadronPt[M_electronmaxcount];
  Float_t electron_r03_sumChargedParticlePt[M_electronmaxcount];
  Float_t electron_r03_sumNeutralHadronEt[M_electronmaxcount];
  Float_t electron_r03_sumPhotonEt[M_electronmaxcount];
  Float_t electron_r03_sumNeutralHadronEtHighThreshold[M_electronmaxcount];
  Float_t electron_r03_sumPhotonEtHighThreshold[M_electronmaxcount];
  Float_t electron_r03_sumPUPt[M_electronmaxcount];

  Float_t electron_charge[M_electronmaxcount];
  UChar_t electron_nhits[M_electronmaxcount];
  UChar_t electron_nmissinghits[M_electronmaxcount];
  UChar_t electron_npixelhits[M_electronmaxcount];
  UChar_t electron_nmissinginnerhits[M_electronmaxcount];
  UChar_t electron_npixellayers[M_electronmaxcount];
  UChar_t electron_nstriplayers[M_electronmaxcount];
  Float_t electron_dxy[M_electronmaxcount];
  Float_t electron_dxyerr[M_electronmaxcount];
  Float_t electron_dz[M_electronmaxcount];
  Float_t electron_dzerr[M_electronmaxcount];
  Float_t electron_convdist[M_electronmaxcount];
  Float_t electron_convdcot[M_electronmaxcount];
  Float_t electron_convradius[M_electronmaxcount];
  UInt_t electron_gapinfo[M_electronmaxcount];
  UInt_t electron_chargeinfo[M_electronmaxcount];
  Float_t electron_fbrems[M_electronmaxcount];
  Int_t electron_numbrems[M_electronmaxcount];
  Int_t electron_superclusterindex[M_electronmaxcount];
  UChar_t electron_info[M_electronmaxcount];

  Float_t electron_mva_id_nontrigPhys14[M_electronmaxcount];
  Float_t electron_mva_value_nontrig_Spring15_v1[M_electronmaxcount];
  Float_t electron_mva_value_trig_Spring15_v1[M_electronmaxcount];
  Int_t electron_mva_category_nontrig_Spring15_v1[M_electronmaxcount];
  Int_t electron_mva_category_trig_Spring15_v1[M_electronmaxcount];

  Bool_t electron_mva_wp80_nontrig_Spring15_v1[M_electronmaxcount];
  Bool_t electron_mva_wp90_nontrig_Spring15_v1[M_electronmaxcount];
  Bool_t electron_mva_wp80_trig_Spring15_v1[M_electronmaxcount];
  Bool_t electron_mva_wp90_trig_Spring15_v1[M_electronmaxcount];

  Bool_t electron_cutId_veto_Spring15[M_electronmaxcount];
  Bool_t electron_cutId_loose_Spring15[M_electronmaxcount];
  Bool_t electron_cutId_medium_Spring15[M_electronmaxcount];
  Bool_t electron_cutId_tight_Spring15[M_electronmaxcount];

  Bool_t electron_pass_conversion[M_electronmaxcount];

  Int_t electron_genmatch[M_electronmaxcount];
  
  UInt_t photon_count;
  Float_t photon_px[M_photonmaxcount];
  Float_t photon_py[M_photonmaxcount];
  Float_t photon_pz[M_photonmaxcount];
  Float_t photon_pt[M_photonmaxcount];
  Float_t photon_eta[M_photonmaxcount];
  Float_t photon_phi[M_photonmaxcount];
  Float_t photon_e1x5[M_photonmaxcount];
  Float_t photon_e2x5[M_photonmaxcount];
  Float_t photon_e3x3[M_photonmaxcount];
  Float_t photon_e5x5[M_photonmaxcount];
  Float_t photon_sigmaetaeta[M_photonmaxcount];
  Float_t photon_sigmaietaieta[M_photonmaxcount];
  Float_t photon_ehcaloverecal[M_photonmaxcount];
  Float_t photon_ehcaloverecaldepth1[M_photonmaxcount];
  Float_t photon_ehcaloverecaldepth2[M_photonmaxcount];
  Float_t photon_maxenergyxtal[M_photonmaxcount];
  Float_t photon_isolationr3track[M_photonmaxcount];
  Float_t photon_isolationr3trackhollow[M_photonmaxcount];
  UInt_t photon_isolationr3ntrack[M_photonmaxcount];
  UInt_t photon_isolationr3ntrackhollow[M_photonmaxcount];
  Float_t photon_isolationr3ecal[M_photonmaxcount];
  Float_t photon_isolationr3hcal[M_photonmaxcount];
  Float_t photon_isolationr4track[M_photonmaxcount];
  Float_t photon_isolationr4trackhollow[M_photonmaxcount];
  UInt_t photon_isolationr4ntrack[M_photonmaxcount];
  UInt_t photon_isolationr4ntrackhollow[M_photonmaxcount];
  Float_t photon_isolationr4ecal[M_photonmaxcount];
  Float_t photon_isolationr4hcal[M_photonmaxcount];
  Int_t photon_superclusterindex[M_photonmaxcount];
  UChar_t photon_info[M_photonmaxcount];
  UInt_t photon_gapinfo[M_photonmaxcount];
  UInt_t photon_conversionbegin[M_photonmaxcount];

  // taus  
  UInt_t tau_count;
  Float_t tau_e[M_taumaxcount];
  Float_t tau_px[M_taumaxcount];
  Float_t tau_py[M_taumaxcount];
  Float_t tau_pz[M_taumaxcount];
  Float_t tau_pt[M_taumaxcount];
  Float_t tau_eta[M_taumaxcount];
  Float_t tau_phi[M_taumaxcount];
  Float_t tau_mass[M_taumaxcount];

  Float_t tau_leadchargedhadrcand_px[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_py[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_pz[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_mass[M_taumaxcount];
  Int_t   tau_leadchargedhadrcand_id[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_dxy[M_taumaxcount];
  Float_t tau_leadchargedhadrcand_dz[M_taumaxcount];


  Float_t tau_vertexx[M_taumaxcount];
  Float_t tau_vertexy[M_taumaxcount];
  Float_t tau_vertexz[M_taumaxcount];

  Float_t tau_dxy[M_taumaxcount];
  Float_t tau_dz[M_taumaxcount];
  Float_t tau_ip3d[M_taumaxcount];
  Float_t tau_ip3dSig[M_taumaxcount];
  Float_t tau_charge[M_taumaxcount];
  Float_t tau_genjet_e[M_taumaxcount];
  Float_t tau_genjet_px[M_taumaxcount];
  Float_t tau_genjet_py[M_taumaxcount];
  Float_t tau_genjet_pz[M_taumaxcount];
  Int_t tau_genmatch[M_taumaxcount];

  // main tau discriminators
  bool setTauBranches;
  std::vector<std::pair<std::string, unsigned int> >tauIdIndx;
  Float_t tau_ids[100][M_taumaxcount];
  
  // number of tracks around 
  UInt_t tau_ntracks_pt05[M_taumaxcount];
  UInt_t tau_ntracks_pt08[M_taumaxcount];
  UInt_t tau_ntracks_pt1[M_taumaxcount];
  
  // l1 match
  Bool_t  tau_L1trigger_match[M_taumaxcount]; 

  UInt_t tau_signalChargedHadrCands_size[M_taumaxcount];
  UInt_t tau_signalNeutralHadrCands_size[M_taumaxcount];
  UInt_t tau_signalGammaCands_size[M_taumaxcount];
  UInt_t tau_isolationChargedHadrCands_size[M_taumaxcount];
  UInt_t tau_isolationNeutralHadrCands_size[M_taumaxcount];
  UInt_t tau_isolationGammaCands_size[M_taumaxcount];

  string tau_genDecayMode_name[M_taumaxcount];
  Int_t tau_genDecayMode[M_taumaxcount];

  string tau_decayMode_name[M_taumaxcount];
  Int_t tau_decayMode[M_taumaxcount];

  // PiZeros
  UInt_t pizero_count;
  Float_t pizero_px[M_pzeromaxcount];
  Float_t pizero_py[M_pzeromaxcount];
  Float_t pizero_pz[M_pzeromaxcount];
  Float_t pizero_pt[M_pzeromaxcount];
  Float_t pizero_eta[M_pzeromaxcount];
  Float_t pizero_phi[M_pzeromaxcount];

  // V0s
  UInt_t v0_count;
  Int_t v0_ID[M_v0maxcount];
  Float_t v0_vx[M_v0maxcount];
  Float_t v0_vy[M_v0maxcount];
  Float_t v0_vz[M_v0maxcount];
  Float_t v0_px[M_v0maxcount];
  Float_t v0_py[M_v0maxcount];
  Float_t v0_pz[M_v0maxcount];
  Float_t v0_pt[M_v0maxcount];
  Float_t v0_eta[M_v0maxcount];
  Float_t v0_phi[M_v0maxcount];
  Float_t v0_mass[M_v0maxcount];
  Float_t v0_chi2[M_v0maxcount];
  Float_t v0_ndof[M_v0maxcount];
  Float_t v0_decay[M_v0maxcount];
  Float_t v0_ip[M_v0maxcount];
  UInt_t  v0_ndaughters[M_v0maxcount];

  Float_t v0_pos_px[M_v0maxcount];
  Float_t v0_pos_py[M_v0maxcount];
  Float_t v0_pos_pz[M_v0maxcount];
  Float_t v0_pos_pt[M_v0maxcount];
  Float_t v0_pos_eta[M_v0maxcount];
  Float_t v0_pos_phi[M_v0maxcount];
  Float_t v0_pos_mass[M_v0maxcount];
  Float_t v0_pos_vx[M_v0maxcount];
  Float_t v0_pos_vy[M_v0maxcount];
  Float_t v0_pos_vz[M_v0maxcount];
  Int_t   v0_pos_ID[M_v0maxcount]; 

  Float_t v0_neg_px[M_v0maxcount];
  Float_t v0_neg_py[M_v0maxcount];
  Float_t v0_neg_pz[M_v0maxcount];
  Float_t v0_neg_pt[M_v0maxcount];
  Float_t v0_neg_eta[M_v0maxcount];
  Float_t v0_neg_phi[M_v0maxcount];
  Float_t v0_neg_mass[M_v0maxcount];
  Float_t v0_neg_vx[M_v0maxcount];
  Float_t v0_neg_vy[M_v0maxcount];
  Float_t v0_neg_vz[M_v0maxcount];
  Int_t   v0_neg_ID[M_v0maxcount]; 

  // generated tau
  UInt_t gentau_count;
  Float_t gentau_px[M_taumaxcount];
  Float_t gentau_py[M_taumaxcount];
  Float_t gentau_pz[M_taumaxcount];
  Float_t gentau_e[M_taumaxcount];

  Float_t gentau_visible_px[M_taumaxcount];
  Float_t gentau_visible_py[M_taumaxcount];
  Float_t gentau_visible_pz[M_taumaxcount];
  Float_t gentau_visible_e[M_taumaxcount];

  Float_t gentau_visible_pt[M_taumaxcount];
  Float_t gentau_visible_eta[M_taumaxcount];
  Float_t gentau_visible_phi[M_taumaxcount];
  Float_t gentau_visible_mass[M_taumaxcount];

  Float_t gentau_visibleNoLep_px[M_taumaxcount];
  Float_t gentau_visibleNoLep_py[M_taumaxcount];
  Float_t gentau_visibleNoLep_pz[M_taumaxcount];
  Float_t gentau_visibleNoLep_e[M_taumaxcount];

  Float_t gentau_visibleNoLep_pt[M_taumaxcount];
  Float_t gentau_visibleNoLep_eta[M_taumaxcount];
  Float_t gentau_visibleNoLep_phi[M_taumaxcount];
  Float_t gentau_visibleNoLep_mass[M_taumaxcount];
  
  Int_t gentau_status[M_taumaxcount];
  Int_t gentau_fromHardProcess[M_taumaxcount];
  Int_t gentau_fromHardProcessBeforeFSR[M_taumaxcount];
  Int_t gentau_isDecayedLeptonHadron[M_taumaxcount];
  Int_t gentau_isDirectHadronDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectHardProcessTauDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectPromptTauDecayProduct[M_taumaxcount];
  Int_t gentau_isDirectTauDecayProduct[M_taumaxcount];
  Int_t gentau_isFirstCopy[M_taumaxcount];
  Int_t gentau_isHardProcess[M_taumaxcount];
  Int_t gentau_isHardProcessTauDecayProduct[M_taumaxcount];
  Int_t gentau_isLastCopy[M_taumaxcount];
  Int_t gentau_isLastCopyBeforeFSR[M_taumaxcount];
  Int_t gentau_isPrompt[M_taumaxcount];
  Int_t gentau_isPromptTauDecayProduct[M_taumaxcount];
  Int_t gentau_isTauDecayProduct[M_taumaxcount];

  Int_t gentau_decayMode[M_taumaxcount];
  string gentau_decayMode_name[M_taumaxcount];
  UChar_t gentau_mother[M_taumaxcount];

  // L1 Iso Tau
  UInt_t l1isotau_count;
  Float_t l1isotau_e[M_taumaxcount];
  Float_t l1isotau_px[M_taumaxcount];
  Float_t l1isotau_py[M_taumaxcount];
  Float_t l1isotau_pz[M_taumaxcount];
  Float_t l1isotau_pt[M_taumaxcount];
  Float_t l1isotau_eta[M_taumaxcount];
  Float_t l1isotau_phi[M_taumaxcount];
  Float_t l1isotau_mass[M_taumaxcount];  
  Float_t l1isotau_charge[M_taumaxcount]; 
  
  // rho neutral
  Float_t rhoNeutral;

  // met
  Float_t pfmet_ex;
  Float_t pfmet_ey;
  Float_t pfmet_ez;
  Float_t pfmet_pt;
  Float_t pfmet_phi;
  Float_t pfmet_sumet;

  Float_t pfmet_sig;
  Float_t pfmet_sigxx;
  Float_t pfmet_sigxy;
  Float_t pfmet_sigyx;
  Float_t pfmet_sigyy;

  Float_t pfmet_ex_JetEnUp;
  Float_t pfmet_ey_JetEnUp;

  Float_t pfmet_ex_JetEnDown;
  Float_t pfmet_ey_JetEnDown;

  Float_t pfmet_ex_UnclusteredEnUp;
  Float_t pfmet_ey_UnclusteredEnUp;

  Float_t pfmet_ex_UnclusteredEnDown;
  Float_t pfmet_ey_UnclusteredEnDown;

  Float_t pfmetcorr_ex;
  Float_t pfmetcorr_ey;
  Float_t pfmetcorr_ez;
  Float_t pfmetcorr_pt;
  Float_t pfmetcorr_phi;
  Float_t pfmetcorr_sumet;

  Float_t pfmetcorr_sig;
  Float_t pfmetcorr_sigxx;
  Float_t pfmetcorr_sigxy;
  Float_t pfmetcorr_sigyx;
  Float_t pfmetcorr_sigyy;

  Float_t pfmetcorr_ex_JetEnUp;
  Float_t pfmetcorr_ey_JetEnUp;

  Float_t pfmetcorr_ex_JetEnDown;
  Float_t pfmetcorr_ey_JetEnDown;

  Float_t pfmetcorr_ex_UnclusteredEnUp;
  Float_t pfmetcorr_ey_UnclusteredEnUp;

  Float_t pfmetcorr_ex_UnclusteredEnDown;
  Float_t pfmetcorr_ey_UnclusteredEnDown;

  Float_t puppimet_ex;
  Float_t puppimet_ey;
  Float_t puppimet_ez;
  Float_t puppimet_pt;
  Float_t puppimet_phi;
  Float_t puppimet_sumet;

  Float_t puppimet_ex_JetEnUp;
  Float_t puppimet_ey_JetEnUp;

  Float_t puppimet_ex_JetEnDown;
  Float_t puppimet_ey_JetEnDown;

  Float_t puppimet_ex_UnclusteredEnUp;
  Float_t puppimet_ey_UnclusteredEnUp;

  Float_t puppimet_ex_UnclusteredEnDown;
  Float_t puppimet_ey_UnclusteredEnDown;

  Float_t puppimet_sigxx;
  Float_t puppimet_sigxy;
  Float_t puppimet_sigyx;
  Float_t puppimet_sigyy;

  UInt_t mvamet_count;
  Float_t mvamet_ex[M_mvametmaxcount];
  Float_t mvamet_ey[M_mvametmaxcount];

  Float_t mvamet_sigxx[M_mvametmaxcount];
  Float_t mvamet_sigxy[M_mvametmaxcount];
  Float_t mvamet_sigyx[M_mvametmaxcount];
  Float_t mvamet_sigyy[M_mvametmaxcount];
  UChar_t mvamet_channel[M_mvametmaxcount];
  UInt_t mvamet_lep1[M_mvametmaxcount];
  UInt_t mvamet_lep2[M_mvametmaxcount];
  Float_t mvamet_lep1_pt[M_mvametmaxcount];
  Float_t mvamet_lep2_pt[M_mvametmaxcount]; 

  Float_t genmet_ex;
  Float_t genmet_ey;

  //Generator Information
  Float_t genweight;
  Float_t genid1;
  Float_t genx1;
  Float_t genid2;
  Float_t genx2;
  Float_t genScale;

  Int_t numpileupinteractionsminus;
  Int_t numpileupinteractions;
  Int_t numpileupinteractionsplus;
  Float_t numtruepileupinteractions;
  Int_t hepNUP_;
  
  Float_t genparticles_lheHt;
  UInt_t genparticles_noutgoing;
  UInt_t genparticles_count;
  Float_t genparticles_e[M_genparticlesmaxcount];
  Float_t genparticles_px[M_genparticlesmaxcount];
  Float_t genparticles_py[M_genparticlesmaxcount];
  Float_t genparticles_pz[M_genparticlesmaxcount];
  Float_t genparticles_vx[M_genparticlesmaxcount];
  Float_t genparticles_vy[M_genparticlesmaxcount];
  Float_t genparticles_vz[M_genparticlesmaxcount];
  Int_t genparticles_pdgid[M_genparticlesmaxcount];
  Int_t genparticles_status[M_genparticlesmaxcount];
  UInt_t genparticles_info[M_genparticlesmaxcount];
  UChar_t genparticles_mother[M_genparticlesmaxcount];

  Int_t genparticles_fromHardProcess[M_genparticlesmaxcount];
  Int_t genparticles_fromHardProcessBeforeFSR[M_genparticlesmaxcount];
  Int_t genparticles_isDecayedLeptonHadron[M_genparticlesmaxcount];
  Int_t genparticles_isDirectHadronDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectHardProcessTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectPromptTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isDirectTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isFirstCopy[M_genparticlesmaxcount];
  Int_t genparticles_isHardProcess[M_genparticlesmaxcount];
  Int_t genparticles_isHardProcessTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isLastCopy[M_genparticlesmaxcount];
  Int_t genparticles_isLastCopyBeforeFSR[M_genparticlesmaxcount];
  Int_t genparticles_isPrompt[M_genparticlesmaxcount];
  Int_t genparticles_isPromptTauDecayProduct[M_genparticlesmaxcount];
  Int_t genparticles_isTauDecayProduct[M_genparticlesmaxcount];

  // trigger objects
  UInt_t trigobject_count;
  Float_t trigobject_px[M_trigobjectmaxcount];
  Float_t trigobject_py[M_trigobjectmaxcount];
  Float_t trigobject_pz[M_trigobjectmaxcount];
  Float_t trigobject_pt[M_trigobjectmaxcount];
  Float_t trigobject_eta[M_trigobjectmaxcount];
  Float_t trigobject_phi[M_trigobjectmaxcount];
  Bool_t  trigobject_filters[M_trigobjectmaxcount][50];
  Bool_t trigobject_isMuon[M_trigobjectmaxcount];
  Bool_t trigobject_isElectron[M_trigobjectmaxcount];
  Bool_t trigobject_isTau[M_trigobjectmaxcount];
  Bool_t trigobject_isJet[M_trigobjectmaxcount];
  Bool_t trigobject_isMET[M_trigobjectmaxcount];

  //lumitree
  UInt_t lumi_run;
  UInt_t lumi_block;
  Float_t lumi_value;
  Float_t lumi_valueerr;
  Float_t lumi_livefrac;
  Float_t lumi_deadfrac;
  UInt_t lumi_quality;
  UInt_t lumi_eventsprocessed;
  UInt_t lumi_eventsfiltered;
  UInt_t lumi_hltprescaletable;
  UInt_t lumi_l1algoprescaletable;
  UInt_t lumi_l1techprescaletable;

  //runtree
  UInt_t run_number;
  UInt_t run_hltcount;
  vector<string> run_hltnames;
  vector<string> run_hltfilters;
  vector<string> run_hltmufilters;
  vector<string> run_hltelectronfilters;
  vector<string> run_hlttaufilters;
  vector<string> run_hltphotonfilters;
  vector<string> run_hltjetfilters;
  vector<string> run_floattaudiscriminators;
  vector<string> run_binarytaudiscriminators;
  vector<string> run_btagdiscriminators;
  UInt_t run_hltprescaletablescount;
  UInt_t run_hltprescaletables[10000];
  UInt_t run_l1algocount;
  UInt_t run_l1algoprescaletablescount;
  UInt_t run_l1algoprescaletables[10000];
  UInt_t run_l1techcount;
  UInt_t run_l1techprescaletablescount;
  UInt_t run_l1techprescaletables[10000];		

  std::map<std::string, int>* hltriggerresults_;
  std::map<std::string, int>* hltriggerprescales_;
  std::vector<std::string>hltriggerresultsV_;
  float embeddingWeight_;
  //std::vector< double > embeddingWeights_; //for RhEmb
  //float TauSpinnerWeight_;
  EGammaMvaEleEstimatorCSA14* myMVAnonTrigPhys14;
  JetCorrectionUncertainty *jecUnc;

};

DEFINE_FWK_MODULE(NTupleMakerAOD);

#endif

#ifndef MyDemoAnalyzer_h
#define MyDemoAnalyzer_h 1

#include "marlin/Processor.h"
#include <marlin/Global.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>
#include <EVENT/LCEvent.h>
#include <EVENT/Cluster.h>
#include <EVENT/CalorimeterHit.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gear/TPCParameters.h>
#include <gear/TPCModule.h>
#include <gear/CalorimeterParameters.h>
#include "IMPL/ClusterImpl.h"
#include "lcio.h"
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Api/PandoraApi.h"
#include "Objects/Helix.h"
#include "Objects/CartesianVector.h"
#include "Pandora/PdgTable.h"
#include "LCHelpers/ClusterHelper.h"

//#include "CalorimeterHitType.h"


using namespace lcio ;
using namespace marlin ;
using namespace std;


/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: MyDemoAnalyzer.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class MyDemoAnalyzer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new MyDemoAnalyzer ; }
  
  
  MyDemoAnalyzer() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;

  
 protected:
  //below two functions are a one to one copy of marlinpandora
    /**
     *  @brief  Obtain track states at start and end of track and the momentum at the dca
     * 
     *  @param  pTrack the lcio track
     *  @param  trackParameters the track parameters
     */
    void GetTrackStatesOld(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const;

    /**
     *  @brief  Project helix to the surface of the ecal
     * 
     *  @param  pHelix helix fit to be projected to ecal surface
     *  @param  signPz sign w.r.t. increasing z direction
     *  @param  trackParameters the track parameters
     */
    void GetECalProjectionOld(const pandora::Helix *const pHelix, const int signPz, PandoraApi::Track::Parameters &trackParameters) const;
    void GetTrackClusterDistance (const PandoraApi::Track::Parameters trackParameters, const EVENT::Cluster *const pCluster, const unsigned int maxSearchLayer, const float parallelDistanceCut, const float minTrackClusterCosAngle, float &trackClusterDistance) const;
 
 

    //needed for above calculations
    float             m_bField;                       ///< The bfield
    int               m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
    float             m_eCalBarrelInnerPhi0;          ///< ECal barrel inner phi 0
    float             m_eCalBarrelInnerR;             ///< ECal barrel inner radius
    float             m_eCalEndCapInnerZ;             ///< ECal endcap inner z
    unsigned int maxSearchLayer;
    float parallelDistanceCut;  
    float minTrackClusterCosAngle;


  /** Input collection name.
   */
  std::string _mcColName ;
  std::string _mcColNameReReco ;
  std::string _jetColName ;
  std::string _jetColNameReReco ;
  std::string _fileName;
  std::string _recoPartName;
  std::string _recoPartNameReReco;
  std::string _recoPartNameReRecoNoTrackSigmaPOverP;
  std::string _trackColName;
  std::string _seltrackColName;
  std::string _v0ColName;
  bool _hasrereco_information;
  bool _hasrerecoNoTrackSigmaPOverP_information;

  float _lepEMin;
  float _phEMin;
  float _jetEMin;
  float _mcPartEMin;
  float _lepGENEMin;
  float _phGENEMin;
  float _jetGENEMin;
  float _nuGENEMin;

  TFile* _outputFile;
  TTree* _tree;
  TTree* _processedEvents;

  int _nRun ;
  int _nEvt ;
  int _numRun ;
  int _numEvt ;
  double _weight;

  std::vector<float> *_jetPx,*_jetPy,*_jetPz,*_jetE,*_jetPhi,*_jetTheta,*_jetMass,*_jetNF,*_jetLF,*_jetKF,*_jetCHF,*_jetPhF,*_jetElF,*_jetMuF;
  std::vector<int> *_jetNMult,*_jetLMult,*_jetKMult,*_jetCHMult,*_jetPhMult,*_jetMuMult,*_jetElMult;
  std::vector<float> *_elPx,*_elPy,*_elPz,*_elE,*_elPhi,*_elTheta,*_elMass,*_elEoverTot,*_elECAL,*_elHCAL,*_elClusterE,*_el0D0,*_el0Z0,*_el0Phi,*_el0Chi2_NDOF,*_elCHIsoDR03,*_elNHIsoDR03,*_elPhIsoDR03,*_elMuIsoDR03,*_elElIsoDR03,*_elCHIsoDR04,*_elNHIsoDR04,*_elPhIsoDR04,*_elMuIsoDR04,*_elElIsoDR04,*_el1D0,*_el1Z0,*_el1Phi,*_el1Chi2_NDOF,*_elCHIsoCosTheta995,*_elNHIsoCosTheta995,*_elPhIsoCosTheta995,*_elMuIsoCosTheta995,*_elElIsoCosTheta995,/*_elLongShowerProfile,*_elTransShowerProfile,*_elTrackBasedTransProfile,*/*_elEnergyEE,*_elEnergyEB,*_elEnergyEElse,*_elEnergyHE,*_elEnergyHB,*_elEnergyHElse,*_elEnergyElseB,*_elEnergyElseE,*_elEnergyElseElse,*_elEnergyEM,*_elEnergyHAD,*_elEnergyMuon,/*_elPeakEnergy,*_elNRadiationLength,*/*_elMaxEInECALLayer,*_elMaxEInHCALLayer,*_elOmega,*_elMinTrackClustDiff,*_elSigmaPOverPTrack,*_elTrackP,*_elTrackPt,*_elTrackClusterOpeningAngle;
  std::vector<int> *_elNHitsVTX,*_elNHitsFTD,*_elNHitsSIT,*_elNHitsTPC,*_elNHitsSET,*_elNHitsETD,*_elNHitsVTXVeto,*_elNHitsFTDVeto,*_elNHitsSITVeto,*_elNHitsTPCVeto,*_elNHitsSETVeto,*_elNHitsETDVeto,*_elID,*_el1NHits,*_el1NHitsVeto,*_elFirstLayerECAL, *_elLastLayerECAL,*_elLayerMaxHitsECAL,*_elHitsECAL,*_elMaxHitsPerLayerECAL,*_elFirstLayerHCAL,*_elLastLayerHCAL,*_elLayerMaxEECAL,*_elMaxHitsPerLayerHCAL,*_elHitsHCAL;
  std::vector<float> *_muPx,*_muPy,*_muPz,*_muE,*_muPhi,*_muTheta,*_muMass,*_muEoverTot,*_muECAL,*_muHCAL,*_muClusterE,*_mu0D0,*_mu0Z0,*_mu0Phi,*_mu0Chi2_NDOF,*_muCHIsoDR03,*_muNHIsoDR03,*_muPhIsoDR03,*_muMuIsoDR03,*_muElIsoDR03,*_muCHIsoDR04,*_muNHIsoDR04,*_muPhIsoDR04,*_muMuIsoDR04,*_muElIsoDR04,*_mu1D0,*_mu1Z0,*_mu1Phi,*_mu1Chi2_NDOF,*_muCHIsoCosTheta995,*_muNHIsoCosTheta995,*_muPhIsoCosTheta995,*_muMuIsoCosTheta995,*_muElIsoCosTheta995,*_muOmega,*_muTrackP,*_muTrackPt,*_muSigmaPOverPTrack;
  std::vector<int> *_muNHitsVTX,*_muNHitsFTD,*_muNHitsSIT,*_muNHitsTPC,*_muNHitsSET,*_muNHitsETD,*_muNHitsVTXVeto,*_muNHitsFTDVeto,*_muNHitsSITVeto,*_muNHitsTPCVeto,*_muNHitsSETVeto,*_muNHitsETDVeto,*_muID,*_mu1NHits,*_mu1NHitsVeto;

  //f1 fraction of first cluster in total energies(for converted photons we have 0-2 clusters), for converted photons track 1/2 phi values
  //try more ee based isolation using an opening angle cos theta > 0.995 degrees, as \vec(p1)*\vec(p2)/(|p1|*|p2|
  std::vector<float> *_phPx,*_phPy,*_phPz,*_phE,*_phPhi,*_phTheta,*_phMass,*_phEoverTot,*_phECAL,*_phHCAL,*_phClusterE,*_phCHIsoDR03,*_phNHIsoDR03,*_phPhIsoDR03,*_phCHIsoDR04,*_phMuIsoDR03,*_phElIsoDR03,*_phNHIsoDR04,*_phPhIsoDR04,*_phMuIsoDR04,*_phElIsoDR04,*_phf1ECAL,*_phf1HCAL,*_phf1ClusterE,*_phT1Chi2_NDOF,*_phT1Phi,*_phT2Chi2_NDOF,*_phT2Phi,*_phCHIsoCosTheta995,*_phNHIsoCosTheta995,*_phPhIsoCosTheta995,*_phMuIsoCosTheta995,*_phElIsoCosTheta995,/*_phLongShowerProfile,*_phTransShowerProfile,*_phTrackBasedTransProfile,*/*_phEnergyEE,*_phEnergyEB,*_phEnergyEElse,*_phEnergyHE,*_phEnergyHB,*_phEnergyHElse,*_phEnergyElseB,*_phEnergyElseE,*_phEnergyElseElse,*_phEnergyEM,*_phEnergyHAD,*_phEnergyMuon,/*_phPeakEnergy,*_phNRadiationLength,*/*_phMaxEInECALLayer,*_phMaxEInHCALLayer,*_phMinTrackClustDiff,*_phSigmaPOverPClosestTrack,*_phClosestTrackP,*_phClosestTrackPt,*_phMinSelTrackClustDiff,*_phSigmaPOverPClosestSelTrack,*_phClosestSelTrackP,*_phClosestSelTrackPt,*_phClosestTrackD0,*_phClosestTrackZ0,*_phClosestTrackPhi,*_phClosestTrackChi2_NDOF,*_phClosestTrackOmega,*_phClosestTrackClusterOpeningAngle,*_phClosestSelTrackClusterOpeningAngle;
  //for converted photons: number of hits used in the fit/total hits of track1/2 phi for track1/2
  std::vector<int> *_phT1NHitsFitted,*_phT1NHitsVeto,*_phT2NHitsFitted,*_phT2NHitsVeto,*_phFirstLayerECAL, *_phLastLayerECAL,*_phLayerMaxHitsECAL,*_phHitsECAL,*_phMaxHitsPerLayerECAL,*_phFirstLayerHCAL,*_phLastLayerHCAL,*_phLayerMaxEECAL,*_phMaxHitsPerLayerHCAL,*_phHitsHCAL,*_phClosestTrackNHitsVTX,*_phClosestTrackNHitsFTD,*_phClosestTrackNHitsSIT,*_phClosestTrackNHitsTPC,*_phClosestTrackNHitsSET,*_phClosestTrackNHitsETD,*_phClosestTrackNHitsVTXVeto,*_phClosestTrackNHitsFTDVeto,*_phClosestTrackNHitsSITVeto,*_phClosestTrackNHitsTPCVeto,*_phClosestTrackNHitsSETVeto,*_phClosestTrackNHitsETDVeto;

  //MET etc MET, px,py,E,pz
  Float_t _MEx,_MEy,_MEz,_MET,_SumE,_SumPt,_SumChargedHadronE,_SumPhotonE,_SumNeutralHadronE,_SumChargedHadronPt,_SumPhotonPt,_SumNeutralHadronPt,_SumMuonE,_SumMuonPt,_SumElectronE,_SumElectronPt;

  //testing new photon reconstruction -->involves also new calibration constants
  //f1 fraction of first cluster in total energies(for converted photons we have 0-2 clusters), for converted photons track 1/2 phi values
  //try more ee based isolation using an opening angle cos theta > 0.995 degrees, as \vec(p1)*\vec(p2)/(|p1|*|p2|

  //now check ReReco stuff --> calibration changed, SHOULD have an effect on any object
  std::vector<float> *_jetPxReReco,*_jetPyReReco,*_jetPzReReco,*_jetEReReco,*_jetPhiReReco,*_jetThetaReReco,*_jetMassReReco,*_jetNFReReco,*_jetLFReReco,*_jetKFReReco,*_jetCHFReReco,*_jetPhFReReco,*_jetElFReReco,*_jetMuFReReco;
  std::vector<int> *_jetNMultReReco,*_jetLMultReReco,*_jetKMultReReco,*_jetCHMultReReco,*_jetPhMultReReco,*_jetMuMultReReco,*_jetElMultReReco;
  std::vector<float> *_elPxReReco,*_elPyReReco,*_elPzReReco,*_elEReReco,*_elPhiReReco,*_elThetaReReco,*_elMassReReco,*_elEoverTotReReco,*_elECALReReco,*_elHCALReReco,*_elClusterEReReco,*_el0D0ReReco,*_el0Z0ReReco,*_el0PhiReReco,*_el0Chi2_NDOFReReco,*_elCHIsoDR03ReReco,*_elNHIsoDR03ReReco,*_elPhIsoDR03ReReco,*_elMuIsoDR03ReReco,*_elElIsoDR03ReReco,*_elCHIsoDR04ReReco,*_elNHIsoDR04ReReco,*_elPhIsoDR04ReReco,*_elMuIsoDR04ReReco,*_elElIsoDR04ReReco,*_el1D0ReReco,*_el1Z0ReReco,*_el1PhiReReco,*_el1Chi2_NDOFReReco,*_elCHIsoCosTheta995ReReco,*_elNHIsoCosTheta995ReReco,*_elPhIsoCosTheta995ReReco,*_elMuIsoCosTheta995ReReco,*_elElIsoCosTheta995ReReco,/**_elLongShowerProfileReReco,*_elTransShowerProfileReReco,*_elTrackBasedTransProfileReReco,*/*_elEnergyEEReReco,*_elEnergyEBReReco,*_elEnergyEElseReReco,*_elEnergyHEReReco,*_elEnergyHBReReco,*_elEnergyHElseReReco,*_elEnergyElseBReReco,*_elEnergyElseEReReco,*_elEnergyElseElseReReco,*_elEnergyEMReReco,*_elEnergyHADReReco,*_elEnergyMuonReReco,/**_elPeakEnergyReReco,*_elNRadiationLengthReReco,*/*_elMaxEInECALLayerReReco,*_elMaxEInHCALLayerReReco,*_elOmegaReReco,*_elMinTrackClustDiffReReco,*_elSigmaPOverPTrackReReco,*_elTrackPReReco,*_elTrackPtReReco,*_elTrackClusterOpeningAngleReReco;
  std::vector<int> *_elNHitsVTXReReco,*_elNHitsFTDReReco,*_elNHitsSITReReco,*_elNHitsTPCReReco,*_elNHitsSETReReco,*_elNHitsETDReReco,*_elNHitsVTXVetoReReco,*_elNHitsFTDVetoReReco,*_elNHitsSITVetoReReco,*_elNHitsTPCVetoReReco,*_elNHitsSETVetoReReco,*_elNHitsETDVetoReReco,*_elIDReReco,*_el1NHitsReReco,*_el1NHitsVetoReReco,*_elFirstLayerECALReReco, *_elLastLayerECALReReco,*_elLayerMaxHitsECALReReco,*_elHitsECALReReco,*_elMaxHitsPerLayerECALReReco,*_elFirstLayerHCALReReco, *_elLastLayerHCALReReco,*_elLayerMaxEECALReReco,*_elMaxHitsPerLayerHCALReReco,*_elHitsHCALReReco;
  std::vector<float> *_muPxReReco,*_muPyReReco,*_muPzReReco,*_muEReReco,*_muPhiReReco,*_muThetaReReco,*_muMassReReco,*_muEoverTotReReco,*_muECALReReco,*_muHCALReReco,*_muClusterEReReco,*_mu0D0ReReco,*_mu0Z0ReReco,*_mu0PhiReReco,*_mu0Chi2_NDOFReReco,*_muCHIsoDR03ReReco,*_muNHIsoDR03ReReco,*_muPhIsoDR03ReReco,*_muMuIsoDR03ReReco,*_muElIsoDR03ReReco,*_muCHIsoDR04ReReco,*_muNHIsoDR04ReReco,*_muPhIsoDR04ReReco,*_muMuIsoDR04ReReco,*_muElIsoDR04ReReco,*_mu1D0ReReco,*_mu1Z0ReReco,*_mu1PhiReReco,*_mu1Chi2_NDOFReReco,*_muCHIsoCosTheta995ReReco,*_muNHIsoCosTheta995ReReco,*_muPhIsoCosTheta995ReReco,*_muMuIsoCosTheta995ReReco,*_muElIsoCosTheta995ReReco,*_muOmegaReReco,*_muTrackPReReco,*_muTrackPtReReco,*_muSigmaPOverPTrackReReco;
  std::vector<int> *_muNHitsVTXReReco,*_muNHitsFTDReReco,*_muNHitsSITReReco,*_muNHitsTPCReReco,*_muNHitsSETReReco,*_muNHitsETDReReco,*_muNHitsVTXVetoReReco,*_muNHitsFTDVetoReReco,*_muNHitsSITVetoReReco,*_muNHitsTPCVetoReReco,*_muNHitsSETVetoReReco,*_muNHitsETDVetoReReco,*_muIDReReco,*_mu1NHitsReReco,*_mu1NHitsVetoReReco;

  //f1 fraction of first cluster in total energies(for converted photons we have 0-2 clusters), for converted photons track 1/2 phi values
  //try more ee based isolation using an opening angle cos theta > 0.995 degrees, as \vec(p1)*\vec(p2)/(|p1|*|p2|
  std::vector<float> *_phPxReReco,*_phPyReReco,*_phPzReReco,*_phEReReco,*_phPhiReReco,*_phThetaReReco,*_phMassReReco,*_phEoverTotReReco,*_phECALReReco,*_phHCALReReco,*_phClusterEReReco,*_phCHIsoDR03ReReco,*_phNHIsoDR03ReReco,*_phPhIsoDR03ReReco,*_phMuIsoDR03ReReco,*_phElIsoDR03ReReco,*_phCHIsoDR04ReReco,*_phNHIsoDR04ReReco,*_phPhIsoDR04ReReco,*_phMuIsoDR04ReReco,*_phElIsoDR04ReReco,*_phf1ECALReReco,*_phf1HCALReReco,*_phf1ClusterEReReco,*_phT1Chi2_NDOFReReco,*_phT1PhiReReco,*_phT2Chi2_NDOFReReco,*_phT2PhiReReco,*_phCHIsoCosTheta995ReReco,*_phNHIsoCosTheta995ReReco,*_phPhIsoCosTheta995ReReco,*_phMuIsoCosTheta995ReReco,*_phElIsoCosTheta995ReReco,/*_phLongShowerProfileReReco,*_phTransShowerProfileReReco,*_phTrackBasedTransProfileReReco,*/*_phEnergyEEReReco,*_phEnergyEBReReco,*_phEnergyEElseReReco,*_phEnergyHEReReco,*_phEnergyHBReReco,*_phEnergyHElseReReco,*_phEnergyElseBReReco,*_phEnergyElseEReReco,*_phEnergyElseElseReReco,*_phEnergyEMReReco,*_phEnergyHADReReco,*_phEnergyMuonReReco,/*_phPeakEnergyReReco,*_phNRadiationLengthReReco,*/*_phMaxEInECALLayerReReco,*_phMaxEInHCALLayerReReco,*_phMinTrackClustDiffReReco,*_phSigmaPOverPClosestTrackReReco,*_phClosestTrackPReReco,*_phClosestTrackPtReReco,*_phMinSelTrackClustDiffReReco,*_phSigmaPOverPClosestSelTrackReReco,*_phClosestSelTrackPReReco,*_phClosestSelTrackPtReReco,*_phClosestTrackD0ReReco,*_phClosestTrackZ0ReReco,*_phClosestTrackPhiReReco,*_phClosestTrackChi2_NDOFReReco,*_phClosestTrackOmegaReReco,*_phClosestTrackClusterOpeningAngleReReco,*_phClosestSelTrackClusterOpeningAngleReReco;
  //for converted photons: number of hits used in the fit/total hits of track1/2 phi for track1/2
  std::vector<int> *_phT1NHitsFittedReReco,*_phT1NHitsVetoReReco,*_phT2NHitsFittedReReco,*_phT2NHitsVetoReReco,*_phFirstLayerECALReReco, *_phLastLayerECALReReco,*_phLayerMaxHitsECALReReco,*_phHitsECALReReco,*_phLayerMaxEECALReReco,*_phMaxHitsPerLayerECALReReco,*_phFirstLayerHCALReReco, *_phLastLayerHCALReReco,*_phMaxHitsPerLayerHCALReReco,*_phHitsHCALReReco,*_phClosestTrackNHitsVTXReReco,*_phClosestTrackNHitsFTDReReco,*_phClosestTrackNHitsSITReReco,*_phClosestTrackNHitsTPCReReco,*_phClosestTrackNHitsSETReReco,*_phClosestTrackNHitsETDReReco,*_phClosestTrackNHitsVTXVetoReReco,*_phClosestTrackNHitsFTDVetoReReco,*_phClosestTrackNHitsSITVetoReReco,*_phClosestTrackNHitsTPCVetoReReco,*_phClosestTrackNHitsSETVetoReReco,*_phClosestTrackNHitsETDVetoReReco;

  /////////////////////////////////////// now THIRD option check Pandora without any cut on Sigma P/P for tracks in checks electron vs photon
  std::vector<float> *_elPxReRecoNoTrkSigmaPOverP,*_elPyReRecoNoTrkSigmaPOverP,*_elPzReRecoNoTrkSigmaPOverP,*_elEReRecoNoTrkSigmaPOverP,*_elPhiReRecoNoTrkSigmaPOverP,*_elThetaReRecoNoTrkSigmaPOverP,*_elMassReRecoNoTrkSigmaPOverP,*_elEoverTotReRecoNoTrkSigmaPOverP,*_elECALReRecoNoTrkSigmaPOverP,*_elHCALReRecoNoTrkSigmaPOverP,*_elClusterEReRecoNoTrkSigmaPOverP,*_el0D0ReRecoNoTrkSigmaPOverP,*_el0Z0ReRecoNoTrkSigmaPOverP,*_el0PhiReRecoNoTrkSigmaPOverP,*_el0Chi2_NDOFReRecoNoTrkSigmaPOverP,*_elCHIsoDR03ReRecoNoTrkSigmaPOverP,*_elNHIsoDR03ReRecoNoTrkSigmaPOverP,*_elPhIsoDR03ReRecoNoTrkSigmaPOverP,*_elMuIsoDR03ReRecoNoTrkSigmaPOverP,*_elElIsoDR03ReRecoNoTrkSigmaPOverP,*_elCHIsoDR04ReRecoNoTrkSigmaPOverP,*_elNHIsoDR04ReRecoNoTrkSigmaPOverP,*_elPhIsoDR04ReRecoNoTrkSigmaPOverP,*_elMuIsoDR04ReRecoNoTrkSigmaPOverP,*_elElIsoDR04ReRecoNoTrkSigmaPOverP,*_el1D0ReRecoNoTrkSigmaPOverP,*_el1Z0ReRecoNoTrkSigmaPOverP,*_el1PhiReRecoNoTrkSigmaPOverP,*_el1Chi2_NDOFReRecoNoTrkSigmaPOverP,*_elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_elMuIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_elElIsoCosTheta995ReRecoNoTrkSigmaPOverP,/**_elLongShowerProfileReRecoNoTrkSigmaPOverP,*_elTransShowerProfileReRecoNoTrkSigmaPOverP,*_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*/*_elEnergyEEReRecoNoTrkSigmaPOverP,*_elEnergyEBReRecoNoTrkSigmaPOverP,*_elEnergyEElseReRecoNoTrkSigmaPOverP,*_elEnergyHEReRecoNoTrkSigmaPOverP,*_elEnergyHBReRecoNoTrkSigmaPOverP,*_elEnergyHElseReRecoNoTrkSigmaPOverP,*_elEnergyElseBReRecoNoTrkSigmaPOverP,*_elEnergyElseEReRecoNoTrkSigmaPOverP,*_elEnergyElseElseReRecoNoTrkSigmaPOverP,*_elEnergyEMReRecoNoTrkSigmaPOverP,*_elEnergyHADReRecoNoTrkSigmaPOverP,*_elEnergyMuonReRecoNoTrkSigmaPOverP,/**_elPeakEnergyReRecoNoTrkSigmaPOverP,*_elNRadiationLengthReRecoNoTrkSigmaPOverP,*/*_elMaxEInECALLayerReRecoNoTrkSigmaPOverP,*_elMaxEInHCALLayerReRecoNoTrkSigmaPOverP,*_elOmegaReRecoNoTrkSigmaPOverP,*_elMinTrackClustDiffReRecoNoTrkSigmaPOverP,*_elSigmaPOverPTrackReRecoNoTrkSigmaPOverP,*_elTrackPReRecoNoTrkSigmaPOverP,*_elTrackPtReRecoNoTrkSigmaPOverP;
 std::vector<int> *_elNHitsVTXReRecoNoTrkSigmaPOverP,*_elNHitsFTDReRecoNoTrkSigmaPOverP,*_elNHitsSITReRecoNoTrkSigmaPOverP,*_elNHitsTPCReRecoNoTrkSigmaPOverP,*_elNHitsSETReRecoNoTrkSigmaPOverP,*_elNHitsETDReRecoNoTrkSigmaPOverP,*_elNHitsVTXVetoReRecoNoTrkSigmaPOverP,*_elNHitsFTDVetoReRecoNoTrkSigmaPOverP,*_elNHitsSITVetoReRecoNoTrkSigmaPOverP,*_elNHitsTPCVetoReRecoNoTrkSigmaPOverP,*_elNHitsSETVetoReRecoNoTrkSigmaPOverP,*_elNHitsETDVetoReRecoNoTrkSigmaPOverP,*_elIDReRecoNoTrkSigmaPOverP,*_el1NHitsReRecoNoTrkSigmaPOverP,*_el1NHitsVetoReRecoNoTrkSigmaPOverP,*_elFirstLayerECALReRecoNoTrkSigmaPOverP, *_elLastLayerECALReRecoNoTrkSigmaPOverP,*_elLayerMaxHitsECALReRecoNoTrkSigmaPOverP,*_elHitsECALReRecoNoTrkSigmaPOverP,*_elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP,*_elFirstLayerHCALReRecoNoTrkSigmaPOverP, *_elLastLayerHCALReRecoNoTrkSigmaPOverP,*_elLayerMaxEECALReRecoNoTrkSigmaPOverP,*_elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP,*_elHitsHCALReRecoNoTrkSigmaPOverP;
  std::vector<float> *_phPxReRecoNoTrkSigmaPOverP,*_phPyReRecoNoTrkSigmaPOverP,*_phPzReRecoNoTrkSigmaPOverP,*_phEReRecoNoTrkSigmaPOverP,*_phPhiReRecoNoTrkSigmaPOverP,*_phThetaReRecoNoTrkSigmaPOverP,*_phMassReRecoNoTrkSigmaPOverP,*_phEoverTotReRecoNoTrkSigmaPOverP,*_phECALReRecoNoTrkSigmaPOverP,*_phHCALReRecoNoTrkSigmaPOverP,*_phClusterEReRecoNoTrkSigmaPOverP,*_phCHIsoDR03ReRecoNoTrkSigmaPOverP,*_phNHIsoDR03ReRecoNoTrkSigmaPOverP,*_phMuIsoDR03ReRecoNoTrkSigmaPOverP,*_phElIsoDR03ReRecoNoTrkSigmaPOverP,*_phPhIsoDR03ReRecoNoTrkSigmaPOverP,*_phCHIsoDR04ReRecoNoTrkSigmaPOverP,*_phNHIsoDR04ReRecoNoTrkSigmaPOverP,*_phPhIsoDR04ReRecoNoTrkSigmaPOverP,*_phMuIsoDR04ReRecoNoTrkSigmaPOverP,*_phElIsoDR04ReRecoNoTrkSigmaPOverP,*_phf1ECALReRecoNoTrkSigmaPOverP,*_phf1HCALReRecoNoTrkSigmaPOverP,*_phf1ClusterEReRecoNoTrkSigmaPOverP,*_phT1Chi2_NDOFReRecoNoTrkSigmaPOverP,*_phT1PhiReRecoNoTrkSigmaPOverP,*_phT2Chi2_NDOFReRecoNoTrkSigmaPOverP,*_phT2PhiReRecoNoTrkSigmaPOverP,*_phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_phMuIsoCosTheta995ReRecoNoTrkSigmaPOverP,*_phElIsoCosTheta995ReRecoNoTrkSigmaPOverP,/*_phLongShowerProfileReRecoNoTrkSigmaPOverP,*_phTransShowerProfileReRecoNoTrkSigmaPOverP,*_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*/*_phEnergyEEReRecoNoTrkSigmaPOverP,*_phEnergyEBReRecoNoTrkSigmaPOverP,*_phEnergyEElseReRecoNoTrkSigmaPOverP,*_phEnergyHEReRecoNoTrkSigmaPOverP,*_phEnergyHBReRecoNoTrkSigmaPOverP,*_phEnergyHElseReRecoNoTrkSigmaPOverP,*_phEnergyElseBReRecoNoTrkSigmaPOverP,*_phEnergyElseEReRecoNoTrkSigmaPOverP,*_phEnergyElseElseReRecoNoTrkSigmaPOverP,*_phEnergyEMReRecoNoTrkSigmaPOverP,*_phEnergyHADReRecoNoTrkSigmaPOverP,*_phEnergyMuonReRecoNoTrkSigmaPOverP,/*_phPeakEnergyReRecoNoTrkSigmaPOverP,*_phNRadiationLengthReRecoNoTrkSigmaPOverP,*/*_phMaxEInECALLayerReRecoNoTrkSigmaPOverP,*_phMaxEInHCALLayerReRecoNoTrkSigmaPOverP,*_phMinTrackClustDiffReRecoNoTrkSigmaPOverP,*_phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP,*_phClosestTrackPReRecoNoTrkSigmaPOverP,*_phClosestTrackPtReRecoNoTrkSigmaPOverP,*_phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP,*_phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP,*_phClosestSelTrackPReRecoNoTrkSigmaPOverP,*_phClosestSelTrackPtReRecoNoTrkSigmaPOverP;
  //for converted photons: number of hits used in the fit/total hits of track1/2 phi for track1/2
  std::vector<int> *_phT1NHitsFittedReRecoNoTrkSigmaPOverP,*_phT1NHitsVetoReRecoNoTrkSigmaPOverP,*_phT2NHitsFittedReRecoNoTrkSigmaPOverP,*_phT2NHitsVetoReRecoNoTrkSigmaPOverP,*_phFirstLayerECALReRecoNoTrkSigmaPOverP, *_phLastLayerECALReRecoNoTrkSigmaPOverP,*_phLayerMaxHitsECALReRecoNoTrkSigmaPOverP,*_phHitsECALReRecoNoTrkSigmaPOverP,*_phLayerMaxEECALReRecoNoTrkSigmaPOverP,*_phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP,*_phFirstLayerHCALReRecoNoTrkSigmaPOverP, *_phLastLayerHCALReRecoNoTrkSigmaPOverP,*_phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP,*_phHitsHCALReRecoNoTrkSigmaPOverP;

 //MET etc MET, px,py,E,pz
  Float_t _MExReReco,_MEyReReco,_MEzReReco,_METReReco,_SumEReReco,_SumPtReReco,_SumChargedHadronEReReco,_SumPhotonEReReco,_SumNeutralHadronEReReco,_SumChargedHadronPtReReco,_SumPhotonPtReReco,_SumNeutralHadronPtReReco,_SumMuonEReReco,_SumMuonPtReReco,_SumElectronEReReco,_SumElectronPtReReco;


  //MC based information
  //
  std::vector<float> *_mcpartStatus,*_mcpartPx,*_mcpartPy,*_mcpartPz,*_mcpartE,*_mcpartPhi,*_mcpartTheta,*_mcpartMass;
  //save status of mothers -> seems to be maybe an issue for photon conversions, also check if daughter is identical with particle, then try not to save the daughter??
  //find the fundamental process, save B's and D's maybe
  std::vector<int> *_mcpartPDGID,*_mcpartParent1PDGID,*_mcpartParent1Stat,*_mcpartParent2PDGID,*_mcpartParent2Stat,*_mcpartDaughter1PDGID,*_mcpartDaughter1Stat,*_mcpartDaughter2PDGID,*_mcpartDaughter2Stat,*_mcpartDaughter3PDGID,*_mcpartDaughter3Stat;
    
  std::vector<float> *_elGENPx,*_elGENPy,*_elGENPz,*_elGENE,*_elGENPhi,*_elGENTheta,*_elGENIsoDR03,*_elGENIsoDR04,*_elGENIsoCosTheta995;
  std::vector<int> *_elGENPDGID,*_elGENParent1PDGID,*_elGENParent2PDGID,*_elGENGrandParent1PDGID;//grandparent1 is parent 1 of parent1 (if existing)

  std::vector<float> *_muGENPx,*_muGENPy,*_muGENPz,*_muGENE,*_muGENPhi,*_muGENTheta,*_muGENIsoDR03,*_muGENIsoDR04,*_muGENIsoCosTheta995;
  std::vector<int> *_muGENPDGID,*_muGENParent1PDGID,*_muGENParent2PDGID,*_muGENGrandParent1PDGID;//grandparent1 is parent 1 of parent1 (if existing)

  std::vector<float> *_phGENPx,*_phGENPy,*_phGENPz,*_phGENE,*_phGENPhi,*_phGENTheta,*_phGENIsoDR03,*_phGENIsoDR04,*_phGENIsoCosTheta995;
  std::vector<int> *_phGENParent1PDGID,*_phGENParent2PDGID,*_phGENGrandParent1PDGID;//grandparent1 is parent 1 of parent1 (if existing)

  std::vector<float> *_nuGENPx,*_nuGENPy,*_nuGENPz,*_nuGENE,*_nuGENPhi,*_nuGENTheta;
  std::vector<int> *_nuGENPDGID,*_nuGENParent1PDGID,*_nuGENParent2PDGID,*_nuGENGrandParent1PDGID;//grandparent1 is parent 1 of parent1 (if existing)

  //GenMET==>neutrino sum
  Float_t _MExGEN,_MEyGEN,_MEzGEN,_METGEN,_MEGEN,_SumEGEN,_SumPtGEN,_SumChargedHadronEGEN,_SumPhotonEGEN,_SumNeutralHadronEGEN,_SumChargedHadronPtGEN,_SumPhotonPtGEN,_SumNeutralHadronPtGEN,_SumMuonEGEN,_SumMuonPtGEN,_SumElectronEGEN,_SumElectronPtGEN;
  //VTX etc
  //int nV;
  //float PVx,PVy,PVz,PVchi2,PVprob;

  // ---- method that builds the tree -------------------------------
  virtual void buildTree();
  // ---- method that re-initializes the tree branches --------------
  virtual void clearTree();

} ;

#endif




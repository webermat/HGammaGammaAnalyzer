#include "MyDemoAnalyzer.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include "Objects/Helix.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDADemoAnalyzer.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif // MARLIN_USE_AIDA

 
using namespace lcio ;
using namespace marlin ;


MyDemoAnalyzer aMyDemoAnalyzer ;


MyDemoAnalyzer::MyDemoAnalyzer() : Processor("MyDemoAnalyzer") {

    // modify processor description
    _description = "MyDemoAnalyzer does whatever it does ..." ;


    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE,
            "MCCollectionName" , 
            "Name of the MCParticle collection"  ,
            _mcColName ,
            std::string("MCParticles")
    );
    registerInputCollection( LCIO::TRACK,
			     "TrackCollectionName",
			     "Name of the Track Collection",
			     _trackColName,
			     std::string("LDCTracks")
			     );
   registerInputCollection( LCIO::TRACK,
			     "SelTrackCollectionName",
			     "Name of the selected Track Collection",
			     _seltrackColName,
			     std::string("SelectedLDCTracks")
			     );


    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE,
            "MCCollectionNameReReco" , 
            "Name of the rerecoed skimmed MCParticle collection"  ,
            _mcColNameReReco ,
            std::string("MCParticles")
    );

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetCollectionName" , 
			     "Name of the RecoJet collection"  ,
			     _jetColName ,
			     std::string("JetOut")
			     );   
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "RecoParticles" , 
            "Name RECOParticle collection",
            _recoPartName ,
            std::string("SelectedPandoraPFANewPFO")
    );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "RecoJetCollectionNameReReco" , 
			     "Name of the ReReco RecoJet collection"  ,
			     _jetColNameReReco ,
			     std::string("JetOut_reprocess")
			     );
    
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "RecoParticlesReReco" , 
            "Name ReReco RECOParticle collection",
            _recoPartNameReReco ,
            std::string("SelectedPandoraPFANewPFO_reprocess")
    );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
            "RecoParticlesReRecoNoTrackSigmaPOverP" , 
            "Name ReReco RECOParticle without track sigma(p)/p qualitiy cuts collection",
			     _recoPartNameReRecoNoTrackSigmaPOverP,
            std::string("SelectedPandoraPFANewPFO_reprocessNoTrackSigmaPOverP")
    );


    registerInputCollection(LCIO::VERTEX,
            "V0Name" , 
            "Name of the PVtx Collection",
	    _v0ColName,          
	    std::string("V0Vertices")
    );



    // register steering parameters: name, description, class-variable, default value
    registerProcessorParameter(
	      "LepEMin" , 
	      "minimum Energy cut for RECO leptons",
	      _lepEMin ,
	      float(10.)
    );
    registerProcessorParameter(
	      "PhEMin" , 
	      "minimum Energy cut for RECO photons",
	      _phEMin ,
	      float(10.)
    );

    registerProcessorParameter(
	      "JetEMin" , 
	      "minimum Energy cut for RECO jets",
	      _jetEMin ,
	      float(20.)
    );
    registerProcessorParameter(
	      "MCPartEMin" , 
	      "minimum Energy cut for unstable MC particles",
	      _mcPartEMin,
	      float(1.)
    );
    registerProcessorParameter(
	      "LepGENEMin" , 
	      "minimum Energy cut for GEN leptons",
	      _lepGENEMin ,
	      float(8.)
    );
    registerProcessorParameter(
	      "PhGENEMin" , 
	      "minimum Energy cut for GEN photons",
	      _phGENEMin ,
	      float(8.)
    );

    registerProcessorParameter(
	      "JetGENEMin" , 
	      "minimum Energy cut for GenJets",
	      _jetGENEMin ,
	      float(17.)
    );
    registerProcessorParameter(
	      "NuGENEMin" , 
	      "minimum Energy cut for GEN neutrinos",
	      _nuGENEMin ,
	      float(6.)
    );

     registerProcessorParameter(
            "FileName" , 
            "Name of the outputfile",
            _fileName ,
            std::string("FileDefaultName.root")
				);
     registerProcessorParameter(
            "HasReRecoInfo" , 
            "rereco workflow present",
            _hasrereco_information ,
	    bool(false)
				);
     registerProcessorParameter(
            "HasReRecoNoTrkSigmaPOverPInfo" , 
            "rereco workflow with no cut on track sigma(p)/p present",				
	    _hasrerecoNoTrackSigmaPOverP_information,
	    bool(false)
				);
}



void MyDemoAnalyzer::init() { 

  //needed for track cluster distance
  m_bField=Global::GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z();
  m_eCalBarrelInnerSymmetry=marlin::Global::GEAR->getEcalBarrelParameters().getSymmetryOrder();
  m_eCalBarrelInnerPhi0=marlin::Global::GEAR->getEcalBarrelParameters().getPhi0();
  m_eCalBarrelInnerR=marlin::Global::GEAR->getEcalBarrelParameters().getExtent()[0];
  m_eCalEndCapInnerZ=marlin::Global::GEAR->getEcalEndcapParameters().getExtent()[2];
  maxSearchLayer=9;
  parallelDistanceCut=100;  
  minTrackClusterCosAngle=0;

  std::cout<<"nu gen e min stuff "<<_nuGENEMin<<std::endl;

    streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
    _outputFile= new TFile(_fileName.c_str(),"recreate");

    _processedEvents = new TTree("ProcessedEvents","ProcessedEvents");
    _processedEvents->Branch("runNum", &_numRun,"runNum/I");
    _processedEvents->Branch("eventNum", &_numEvt, "eventNum/I");
    _processedEvents->Branch("weight", &_weight, "weight/D");

    _tree = new TTree("events","events");
    buildTree();

    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}


void MyDemoAnalyzer::processRunHeader( LCRunHeader* run) { 

    _nRun++ ;
} 



void MyDemoAnalyzer::processEvent( LCEvent * evt ) { 

  clearTree();

  _numEvt=evt->getEventNumber();
  _numRun=evt->getRunNumber();
  _weight= evt->getWeight();

  _processedEvents->Fill();

  if(evt->getEventNumber()%100 ==0){
    std::cout<<"in run/event/weight "<<evt->getRunNumber ()<<"/"<<evt->getEventNumber()<<"/"<<evt->getWeight()<<"/"<<_lepEMin<<"/"<<"/"<<_phEMin<<"/"<<_jetEMin<<std::endl;
  }

#ifdef MARLIN_USE_AIDA

    // define a histogram pointer
    static AIDA::ICloud1D* hMCPEnergy ;    

    if( isFirstEvent() ) { 

        hMCPEnergy = AIDADemoAnalyzer::histogramFactory(this)->
            createCloud1D( "hMCPEnergy", "energy of the MCParticles", 100 ) ; 

    }
#endif // MARLIN_USE_AIDA

    // try to get lcio collection (exits if collection is not available)
    // NOTE: if the AIDADemoAnalyzer is activated in your steering file and Marlin is linked with
    //      RAIDA you may get the message: "*** Break *** segmentation violation" followed by a
    //      stack-trace (generated by ROOT) in case the collection is unavailable. This happens
    //      because ROOT is somehow catching the exit signal commonly used to exit a program
    //      intentionally. Sorry if this messsage may confuse you. Please ignore it!

    LCCollection* mccol = NULL;
    // Alternativelly if you do not want Marlin to exit in case of a non-existing collection
    // use the following (commented out) code:
    //run on H to gamma gamma -> in case there are no tracks around linker will fail to produce output collection
    try{
      mccol = evt->getCollection( _mcColNameReReco ) ;
    }
    catch( lcio::DataNotAvailableException e )
      {
	streamlog_out(WARNING) << _mcColNameReReco << " collection not available, use instead" << _mcColName << std::endl;
	mccol = evt->getCollection( _mcColName ) ;
      }

    // this will only be entered if the collection is available
    if( mccol != NULL ){
      for(int i=0; i< mccol->getNumberOfElements(); i++){
	MCParticle* mcpart = dynamic_cast<MCParticle*>( mccol->getElementAt( i ) ) ;
	//avoid later double filling
	if(mcpart->getGeneratorStatus()!=1){
	  //quarks, leptons, SUSY, gauge bosons/D and B mesons --> check on energy of photons and gluons to avoid being swamped by EM or QCD parton shower stuff
	  if( (abs(mcpart->getPDG())<10 || ( abs(mcpart->getPDG())>10 && abs(mcpart->getPDG())<20 ) || 
	     (( (mcpart->getPDG())==22|| (mcpart->getPDG())==21) && (mcpart->getEnergy())>_mcPartEMin ) || (abs(mcpart->getPDG())>22 && abs(mcpart->getPDG())<50) 
	     ||(abs(mcpart->getPDG())>1000000 && abs(mcpart->getPDG())<3000000)
	     || (abs(mcpart->getPDG())>410 && abs(mcpart->getPDG())<440) || (abs(mcpart->getPDG())>510 && abs(mcpart->getPDG())<550))  && (mcpart->getGeneratorStatus())==2 ){
	    bool veto_mcparticle=false;
	    //check if we have the following case -> we find a daughter, which is basically the particle itself (only one daughter, same momentum, i.e. baiscally identical then don't save!)--> often just a change in status
	    if(mcpart->getParents().size()>0 && mcpart->getParents()[0]->getGeneratorStatus()==1){
	      veto_mcparticle=true;
	    }
	    if(mcpart->getDaughters().size()==1 && (mcpart->getPDG()==mcpart->getDaughters()[0]->getPDG())){
	      TLorentzVector mctempP;
	      mctempP.SetPxPyPzE(mcpart->getMomentum()[0],mcpart->getMomentum()[1],mcpart->getMomentum()[2],mcpart->getEnergy());
	      TLorentzVector mctempD;
	      mctempD.SetPxPyPzE(mcpart->getDaughters()[0]->getMomentum()[0],mcpart->getDaughters()[0]->getMomentum()[1],mcpart->getDaughters()[0]->getMomentum()[2],mcpart->getDaughters()[0]->getEnergy());
	      double deltaRcheck=mctempP.DeltaR(mctempD);
	      double costhetacheck=((mctempP.Px()*mctempD.Px()+mctempP.Py()*mctempD.Py()+mctempP.Pz()*mctempD.Pz())/(mctempP.P()*mctempD.P()));
	      if(deltaRcheck < 10e-5 && fabs(costhetacheck-1.)<10e-5){
		veto_mcparticle=true;
	      }
	    }
	    if(!veto_mcparticle){
	      _mcpartStatus->push_back(mcpart->getGeneratorStatus());
	      _mcpartPDGID->push_back(mcpart->getPDG());
	      _mcpartPx->push_back(mcpart->getMomentum()[0]);
	      _mcpartPy->push_back(mcpart->getMomentum()[1]);
	      _mcpartPz->push_back(mcpart->getMomentum()[2]);
	      _mcpartE->push_back(mcpart->getEnergy());
	      _mcpartMass->push_back(mcpart->getMass());
	      TLorentzVector mctemp;
	      mctemp.SetPxPyPzE(mcpart->getMomentum()[0],mcpart->getMomentum()[1],mcpart->getMomentum()[2],mcpart->getEnergy());
	      _mcpartPhi->push_back(mctemp.Phi());
	      _mcpartTheta->push_back(mctemp.Theta());
	      _mcpartParent1PDGID->push_back(0);
	      _mcpartParent2PDGID->push_back(0);
	      if(mcpart->getParents().size()>0){
		(*_mcpartParent1PDGID)[_mcpartParent1PDGID->size()-1]=mcpart->getParents()[0]->getPDG();
		if(mcpart->getParents().size()>1){
		  (*_mcpartParent2PDGID)[_mcpartParent2PDGID->size()-1]=mcpart->getParents()[1]->getPDG();
		}
	      }
	      _mcpartDaughter1PDGID->push_back(0);
	      _mcpartDaughter2PDGID->push_back(0);
	      _mcpartDaughter3PDGID->push_back(0);
	      if(mcpart->getDaughters().size()>0){
		//if(mcpart->getPDG()==25){
		  //std::cout<<"daughter 0/status "<<mcpart->getDaughters()[0]->getPDG()<<"/"<<mcpart->getDaughters()[0]->getGeneratorStatus()<<std::endl;
		//}
		(*_mcpartDaughter1PDGID)[_mcpartDaughter1PDGID->size()-1]=mcpart->getDaughters()[0]->getPDG();
		if(mcpart->getDaughters().size()>1){
		  //if(mcpart->getPDG()==25){
		  //std::cout<<"daughter 1/status "<<mcpart->getDaughters()[1]->getPDG()<<"/"<<mcpart->getDaughters()[1]->getGeneratorStatus()<<std::endl;
		  //}
		  (*_mcpartDaughter2PDGID)[_mcpartDaughter2PDGID->size()-1]=mcpart->getDaughters()[1]->getPDG();
		  if(mcpart->getDaughters().size()>2){
		    (*_mcpartDaughter3PDGID)[_mcpartDaughter3PDGID->size()-1]=mcpart->getDaughters()[2]->getPDG();
		  }
		}
	      }
	    }
	  //std::vector<int> *_mcpartPDGID,*_mcpartParent1PDGID,*_mcpartParent1Stat,*_mcpartParent2PDGID,*_mcpartParent2Stat,*_mcpartDaughter1PDGID,*_mcpartDaughter1Stat,*_mcpartDaughter2PDGID,*_mcpartDaughter2Stat,*_mcpartDaughter3PDGID,*_mcpartDaughter3Stat;
	  }
	}
	if(mcpart->getGeneratorStatus()==1){
	  float sumGEN_DR03=0;
	  float sumGEN_DR04=0;
	  float sumGEN_cosTheta995=0;
	  TLorentzVector temp;
	  temp.SetPxPyPzE(mcpart->getMomentum()[0],mcpart->getMomentum()[1],mcpart->getMomentum()[2],mcpart->getEnergy());
	  if((abs(mcpart->getPDG())==13 && mcpart->getEnergy()>=_lepGENEMin) || (abs(mcpart->getPDG())==11 && mcpart->getEnergy()>=_lepGENEMin)|| (mcpart->getPDG()==22 && mcpart->getEnergy()>=_phGENEMin)){
	    for(int j=0; j< mccol->getNumberOfElements(); j++){
	      MCParticle* mcpart2 = dynamic_cast<MCParticle*>( mccol->getElementAt( j ) ) ;
	      if(i!=j){
		if(mcpart2->getGeneratorStatus()==1){//exclude neutrinos and neutralinos/gravitinos from isolation sum
		  if(abs(mcpart2->getPDG()!=12) && abs(mcpart2->getPDG()!=14) && abs(mcpart2->getPDG()!=16) && abs(mcpart2->getPDG())!=1000022 && abs(mcpart2->getPDG())!=1000039){
		    TLorentzVector temp2;
		    temp2.SetPxPyPzE(mcpart2->getMomentum()[0],mcpart2->getMomentum()[1],mcpart2->getMomentum()[2],mcpart2->getEnergy());
		    double deltaR=temp2.DeltaR(temp);
		    double cosTheta=(temp.Px()*temp2.Px()+temp.Py()*temp2.Py()+temp.Pz()*temp2.Pz())/(temp.P()*temp2.P());
		    if(deltaR<0.4){
		      sumGEN_DR04+=mcpart2->getEnergy();
		      if(deltaR<0.3){
			sumGEN_DR03+=mcpart2->getEnergy();
		      }
		    }
		    if(cosTheta>0.995){
		      sumGEN_cosTheta995+=mcpart2->getEnergy();
		    }
		  }
		}
	      }
	    }
	  }
	  //neutrinos, neutralino,gravitino
	  if((abs(mcpart->getPDG())==12)|| (abs(mcpart->getPDG())==14) || (abs(mcpart->getPDG())==16) || (abs(mcpart->getPDG())==1000022) || (abs(mcpart->getPDG())==1000039)){
	    _MExGEN+=mcpart->getMomentum()[0];
	    _MEyGEN+=mcpart->getMomentum()[1];
	    _MEzGEN+=mcpart->getMomentum()[2];
	    _MEGEN+=mcpart->getEnergy();
	    if(((abs(mcpart->getPDG())==12)||(abs(mcpart->getPDG())==14) || (abs(mcpart->getPDG())==16)) && (mcpart->getEnergy())>_nuGENEMin){
	      _nuGENPx->push_back(mcpart->getMomentum()[0]);
	      _nuGENPy->push_back(mcpart->getMomentum()[1]);
	      _nuGENPz->push_back(mcpart->getMomentum()[2]);
	      _nuGENE ->push_back(mcpart->getEnergy());
	      _nuGENPhi->push_back(temp.Phi());
	      _nuGENTheta->push_back(temp.Theta());
	      _nuGENPDGID->push_back(mcpart->getPDG());
	      _nuGENParent1PDGID->push_back(0);
	      _nuGENGrandParent1PDGID->push_back(0);
	      _nuGENParent2PDGID->push_back(0);
	      if(mcpart->getParents().size()>0){
		(*_nuGENParent1PDGID)[_nuGENParent1PDGID->size()-1]=mcpart->getParents()[0]->getPDG();
		if(mcpart->getParents()[0]->getParents().size()>0){
		  (*_nuGENGrandParent1PDGID)[_nuGENGrandParent1PDGID->size()-1]=mcpart->getParents()[0]->getParents()[0]->getPDG();
		}
		if(mcpart->getParents().size()>1){
		  (*_nuGENParent2PDGID)[_nuGENParent2PDGID->size()-1]=mcpart->getParents()[1]->getPDG();
		}
	      }
	    }
	  }//status 1 is set before
	  if(((abs(mcpart->getPDG())!=12)&&(abs(mcpart->getPDG())!=14) && (abs(mcpart->getPDG())!=16)) && (abs(mcpart->getPDG())!=1000022) && (abs(mcpart->getPDG())!=1000039)){
	    _SumEGEN+=mcpart->getEnergy();
	    _SumPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	    if(abs(mcpart->getPDG())==11){
	      _SumElectronEGEN+=mcpart->getEnergy();
	      _SumElectronPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	    }else if(abs(mcpart->getPDG())==13){
	      _SumMuonEGEN+=mcpart->getEnergy();
	      _SumMuonPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	    }else if(mcpart->getPDG()==22){
	      _SumPhotonEGEN+=mcpart->getEnergy();
	      _SumPhotonPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	    }else if(((abs(mcpart->getPDG())< 1000000) || (abs(mcpart->getPDG())> 4000000) ) && mcpart->getCharge()==0){//take SUSY particles out, check charge
	      _SumNeutralHadronEGEN+=mcpart->getEnergy();
	      _SumNeutralHadronPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	      //std::cout<<"neutral hadron in GEN particle "<<mcpart->getPDG()<<std::endl;
	    }else if(((abs(mcpart->getPDG())< 1000000) || (abs(mcpart->getPDG())> 4000000) ) && mcpart->getCharge()!=0){//take SUSY particles out, check charge
	      _SumChargedHadronEGEN+=mcpart->getEnergy();
	      _SumChargedHadronPtGEN+=sqrt(mcpart->getMomentum()[0]*mcpart->getMomentum()[0]+mcpart->getMomentum()[1]*mcpart->getMomentum()[1]);
	    }
	  }
	  //check electrons, but veto electrons from photon conversions (i.e. veto electons with one mother which is a photon and status 1)
	  if(abs(mcpart->getPDG())==11 && mcpart->getEnergy()>_lepGENEMin){
	    bool veto_photon_conversion=false;
	    if(mcpart->getParents().size()==1 && ((mcpart->getParents()[0]->getPDG())==22 && (mcpart->getParents()[0]->getGeneratorStatus())==1)){
	      veto_photon_conversion=true;
	      //std::cout<<"we have a photon conversion somewhere"<<std::endl;
	    }
	    if(!veto_photon_conversion){//electrons from conversions seem to have not status 1, most have no daughters, consider for now irrelevant
	      _elGENPx->push_back(mcpart->getMomentum()[0]);
	      _elGENPy->push_back(mcpart->getMomentum()[1]);
	      _elGENPz->push_back(mcpart->getMomentum()[2]);
	      _elGENE ->push_back(mcpart->getEnergy());
	      _elGENPhi->push_back(temp.Phi());
	      _elGENTheta->push_back(temp.Theta());
	      _elGENPDGID->push_back(mcpart->getPDG());
	      _elGENParent1PDGID->push_back(0);
	      _elGENGrandParent1PDGID->push_back(0);
	      _elGENParent2PDGID->push_back(0);
	      _elGENIsoDR03->push_back(sumGEN_DR03);
	      _elGENIsoDR04->push_back(sumGEN_DR04);
	      _elGENIsoCosTheta995->push_back(sumGEN_cosTheta995);
	      if(mcpart->getParents().size()>0){
		(*_elGENParent1PDGID)[_elGENParent1PDGID->size()-1]=mcpart->getParents()[0]->getPDG();
		if(mcpart->getParents()[0]->getParents().size()>0){
		  (*_elGENGrandParent1PDGID)[_elGENGrandParent1PDGID->size()-1]=mcpart->getParents()[0]->getParents()[0]->getPDG();
		}
		if(mcpart->getParents().size()>1){
		  (*_elGENParent2PDGID)[_elGENParent2PDGID->size()-1]=mcpart->getParents()[1]->getPDG();
		}
	      }
	    }
	  }
	  if(abs(mcpart->getPDG())==13 && (mcpart->getEnergy())>_lepGENEMin){
	    _muGENPx->push_back(mcpart->getMomentum()[0]);
	    _muGENPy->push_back(mcpart->getMomentum()[1]);
	    _muGENPz->push_back(mcpart->getMomentum()[2]);
	    _muGENE ->push_back(mcpart->getEnergy());
	    _muGENPhi->push_back(temp.Phi());
	    _muGENTheta->push_back(temp.Theta());
	    _muGENPDGID->push_back(mcpart->getPDG());
	    _muGENParent1PDGID->push_back(0);
	    _muGENParent2PDGID->push_back(0);
	    _muGENGrandParent1PDGID->push_back(0);
	    _muGENIsoDR03->push_back(sumGEN_DR03);
	    _muGENIsoDR04->push_back(sumGEN_DR04);
	    _muGENIsoCosTheta995->push_back(sumGEN_cosTheta995);
	    if(mcpart->getParents().size()>0){
	      (*_muGENParent1PDGID)[_muGENParent1PDGID->size()-1]=mcpart->getParents()[0]->getPDG();
	      if(mcpart->getParents()[0]->getParents().size()>0){
		(*_muGENGrandParent1PDGID)[_muGENGrandParent1PDGID->size()-1]=mcpart->getParents()[0]->getParents()[0]->getPDG();
	      }
	      if(mcpart->getParents().size()>1){
		(*_muGENParent2PDGID)[_muGENParent2PDGID->size()-1]=mcpart->getParents()[1]->getPDG();
	      }
	    }
	  }
	  if(mcpart->getPDG()==22 && (mcpart->getEnergy())>_phGENEMin){
	    _phGENPx->push_back(mcpart->getMomentum()[0]);
	    _phGENPy->push_back(mcpart->getMomentum()[1]);
	    _phGENPz->push_back(mcpart->getMomentum()[2]);
	    _phGENE ->push_back(mcpart->getEnergy());
	    _phGENPhi->push_back(temp.Phi());
	    _phGENTheta->push_back(temp.Theta());
	    _phGENParent1PDGID->push_back(0);
	    _phGENGrandParent1PDGID->push_back(0);
	    _phGENParent2PDGID->push_back(0);
	    _phGENIsoDR03->push_back(sumGEN_DR03);
	    _phGENIsoDR04->push_back(sumGEN_DR04);
	    _phGENIsoCosTheta995->push_back(sumGEN_cosTheta995);
	    if(mcpart->getParents().size()>0){
	      (*_phGENParent1PDGID)[_phGENParent1PDGID->size()-1]=mcpart->getParents()[0]->getPDG();
	      if(mcpart->getParents()[0]->getParents().size()>0){
		(*_phGENGrandParent1PDGID)[_phGENGrandParent1PDGID->size()-1]=mcpart->getParents()[0]->getParents()[0]->getPDG();
	      }
	      if(mcpart->getParents().size()>1){
		(*_phGENParent2PDGID)[_phGENParent2PDGID->size()-1]=mcpart->getParents()[1]->getPDG();
	      }
	    }
	  }	
	}
	//std::vector<float> *_mcpartStatus,*_mcpartPx,*_mcpartPy,*_mcpartPz,*_mcpartE,*_mcpartPhi,*_mcpartTheta,*_mcpartMass,*_mcpartIsoDR03,*_mcpartIsoDR04,*_mcpartIsoCosTheta995;
	//std::vector<int> *_mcpartPDGID,*_mcpartParent1PDGID,*_mcpartParent2PDGID,*_mcpartDaughter1PDGID,*_mcpartDaughter2PDGID,*_mcpartDaughter3PDGID,

      }
      _METGEN=sqrt(pow(_MExGEN,2)+pow(_MEyGEN,2));
    }

    LCCollection* recoPartCol = evt->getCollection( _recoPartName ) ;
    LCCollection* trackCol = evt->getCollection( _trackColName  ) ;
    LCCollection* seltrackCol = evt->getCollection( _seltrackColName  ) ;
    
    _MEx=0;_MEy=0;_MEz=0;_MET=0;_SumE=0;_SumPt=0;_SumChargedHadronE=0;_SumPhotonE=0;_SumNeutralHadronE=0;_SumChargedHadronPt=0;_SumPhotonPt=0;_SumNeutralHadronPt=0;_SumMuonE=0;_SumMuonPt=0;_SumElectronE=0;_SumElectronPt=0;
    if(recoPartCol!=NULL){
      //PandoraCandidate loop
     for(int i=0;i<recoPartCol->getNumberOfElements();i++){
       //std::cout<<"in reco particle stuff"<<std::endl;
       ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoPartCol->getElementAt(i));
       //if(pandorapart->getType()==130 || pandorapart->getType()==310){
       //std::cout<<"in reco check type 130 pr 310 "<<pandorapart->getType()<<"/"<<pandorapart->getMass()<<std::endl;
       //}
       _MEx-=pandorapart->getMomentum()[0];
       _MEy-=pandorapart->getMomentum()[1];
       _MEz-=pandorapart->getMomentum()[2];
       _SumE+=pandorapart->getEnergy();
       _SumPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       TLorentzVector temp;
       temp.SetPxPyPzE(pandorapart->getMomentum()[0],pandorapart->getMomentum()[1],pandorapart->getMomentum()[2],pandorapart->getEnergy());
       if(abs(pandorapart->getType())==211){
	 _SumChargedHadronE+=pandorapart->getEnergy();
	 _SumChargedHadronPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       }else if(pandorapart->getType()==22){
	 _SumPhotonE+=pandorapart->getEnergy();
	 _SumPhotonPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       }else if(abs(pandorapart->getType())==11){
	 _SumElectronE+=pandorapart->getEnergy();
	 _SumElectronPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       }else if(abs(pandorapart->getType())==13){
	 _SumMuonE+=pandorapart->getEnergy();
	 _SumMuonPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       }else{
	 //if(pandorapart->getType()==130){
	 //std::cout<<"in reco neutral hadrons 130 "<<pandorapart->getType()<<"/"<<pandorapart->getMass()<<std::endl;
	 //}else{
	 //std::cout<<"in reco neutral hadrons "<<pandorapart->getType()<<"/"<<pandorapart->getMass()<<std::endl;
	 //}
	 _SumNeutralHadronE+=pandorapart->getEnergy();
	 _SumNeutralHadronPt+=sqrt(pow(pandorapart->getMomentum()[1],2)+pow(pandorapart->getMomentum()[2],2));
       }
       if(abs(pandorapart->getType())==11 && pandorapart->getEnergy()>=_lepEMin){
	 //std::cout<<"in reco particle stuff e"<<std::endl;
	 _elID->push_back(pandorapart->getType());
	 _elE->push_back(pandorapart->getEnergy());
	 _elPx->push_back(pandorapart->getMomentum()[0]);
	 _elPy->push_back(pandorapart->getMomentum()[1]);
	 _elPz->push_back(pandorapart->getMomentum()[2]);
	 _elPhi->push_back(temp.Phi());
	 _elTheta->push_back(temp.Theta());
	 _elMass->push_back(pandorapart->getMass());
	 //start here electoton specific cluster/composition/shower shape etc etc
	 float En_EM=0;
	 float En_HAD=0;
	 float En_MUON=0;	  
	 float En_ECAL_Barrel=0;
	 float En_ECAL_Endcap=0;
	 float En_ECAL_else=0;
	 float En_HCAL_Barrel=0;
	 float En_HCAL_Endcap=0;
	 float En_HCAL_else=0;
	 float En_ELSE_Barrel=0;
	 float En_ELSE_Endcap=0;
	 float En_ELSE_else=0;
	 //*_elLongShowerProfile,*_elTransShowerProfile,*_elTrackBasedTransProfile,*_elPeakEnergy,*_elNRadiationLength;
	 int elFirstLayerHCAL=-1;
	 int elFirstLayerECAL=-1;
	 int elLastLayerHCAL=-1;
	 int elLastLayerECAL=-1;
	 int hitsECAL=0;
	 int hitsHCAL=0;
	 //40 layers decided at the moment, 20+10 in samples
	 std::vector<int>hits_per_layer_ecal(50,0);
	 std::vector<float>energy_per_layer_ecal(50,0);
	 //75 layers at the moment
	 std::vector<int>hits_per_layer_hcal(100,0);
	 std::vector<float>energy_per_layer_hcal(100,0);
	 for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	   for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
	     TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
	     tempP.SetXYZT(pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
	     int type=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
	     static const int fCaloType =     1 ;
	     static const int fCaloID   =    10 ;
	     static const int fLayout   =  1000 ;
	     static const int fLayer    = 10000 ;
	     int caloLayer=type/fLayer;
	     int caloLayout=(type-caloLayer*fLayer)/fLayout;
	     int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
	     int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
	     //std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
	     if(caloType==0){//0 em, 1 had, 2 mu
	       En_EM+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloType==1){     
	       En_HAD+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==2){
	       En_MUON+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	     if(caloID==1){//ecal
	       if(elFirstLayerECAL<0){
		 elFirstLayerECAL=caloLayer;
	       }else if (caloLayer<elFirstLayerECAL){
		 elFirstLayerECAL=caloLayer;
	       }
	       if(elLastLayerECAL<0){
		 elLastLayerECAL=caloLayer;
	       }else if (caloLayer>elLastLayerECAL){
		 elLastLayerECAL=caloLayer;
	       }
	       hitsECAL+=1;
	       hits_per_layer_ecal[caloLayer]+=1;
	       energy_per_layer_ecal[caloLayer]+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if (caloID==2){
	       if(elFirstLayerHCAL<0){
		 elFirstLayerHCAL=caloLayer;
	       }else if (caloLayer<elFirstLayerHCAL){
		 elFirstLayerHCAL=caloLayer;
	       }
	       if(elLastLayerHCAL<0){
		 elLastLayerHCAL=caloLayer;
	       }else if (caloLayer>elLastLayerHCAL){
		 elLastLayerHCAL=caloLayer;
	       }
	       hitsHCAL+=1;
	       hits_per_layer_hcal[caloLayer]+=1;
	       energy_per_layer_hcal[caloLayer]+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	     if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
	       En_ECAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==1 && caloLayout==2){
	       En_ECAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==1){
	       En_ECAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==2 && caloLayout==1){
	       En_HCAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==2 && caloLayout==2){
	       En_HCAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==2){
	       En_HCAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if((caloID!=1 && caloID!=2) && caloLayout==1){
	       En_ELSE_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if((caloID!=1 && caloID!=2) && caloLayout==2){
	       En_ELSE_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID!=1 && caloID!=2){
	       En_ELSE_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	   }
	 }//loop over all clusters and hits for energy in different detectors
	 //*_elLongShowerProfile,*_elTransShowerProfile,*_elTrackBasedTransProfile,*_elPeakEnergy,*_elNRadiationLength;
	 int maxHitsPerLayer=0;
	 float maxEnergyPerLayer=0;
	 _elHitsECAL->push_back(hitsECAL);
	 _elFirstLayerECAL->push_back(elFirstLayerECAL);
	 _elLastLayerECAL->push_back(elLastLayerECAL);
	 if(hitsECAL>0){
	   _elLayerMaxHitsECAL->push_back(-1);
	   _elLayerMaxEECAL->push_back(-1);
	   for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
	     if(hits_per_layer_ecal[j]>maxHitsPerLayer){
	       maxHitsPerLayer=hits_per_layer_ecal[j];
	       (*_elLayerMaxHitsECAL)[_elLayerMaxHitsECAL->size()-1]=j;
	     }
	     if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
	       maxEnergyPerLayer=energy_per_layer_ecal[j];
	       (*_elLayerMaxEECAL)[_elLayerMaxEECAL->size()-1]=j;
	     }
	   }
	   _elMaxHitsPerLayerECAL->push_back(maxHitsPerLayer);
	   _elMaxEInECALLayer->push_back(maxEnergyPerLayer);
	 }else{//no ECAL hits
	   _elLayerMaxHitsECAL->push_back(-1);
	   _elLayerMaxEECAL->push_back(-1);
	   _elMaxHitsPerLayerECAL->push_back(0);
	   _elMaxEInECALLayer->push_back(0);
	 }
	 _elHitsHCAL->push_back(hitsHCAL);
	 _elFirstLayerHCAL->push_back(elFirstLayerHCAL);
	 _elLastLayerHCAL->push_back(elLastLayerHCAL);
	 //initialize to 0 again for HCAL
	 maxHitsPerLayer=0;
	 maxEnergyPerLayer=0;
	 if(hitsHCAL>0){
	   for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
	     if(hits_per_layer_ecal[j]>maxHitsPerLayer){
	       maxHitsPerLayer=hits_per_layer_ecal[j];
	     }
	     if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
	       maxEnergyPerLayer=energy_per_layer_ecal[j];
	     }
	   }
	   _elMaxHitsPerLayerHCAL->push_back(maxHitsPerLayer);
	   _elMaxEInHCALLayer->push_back(maxEnergyPerLayer);
	 }else{//no HCAL hits
	   _elMaxHitsPerLayerHCAL->push_back(0);
	   _elMaxEInHCALLayer->push_back(0);
	 }
	 _elEnergyEB->push_back(En_ECAL_Barrel);
	 _elEnergyEE->push_back(En_ECAL_Endcap);
	 _elEnergyEElse->push_back(En_ECAL_else);
	 _elEnergyHE->push_back(En_HCAL_Barrel);
	 _elEnergyHB->push_back(En_HCAL_Endcap);
	 _elEnergyHElse->push_back( En_HCAL_else);
	 _elEnergyElseB->push_back(En_ELSE_Barrel);
	 _elEnergyElseE->push_back(En_ELSE_Endcap);
	 _elEnergyElseElse->push_back(En_ELSE_else);
	 _elEnergyEM->push_back(En_EM);
	 _elEnergyHAD->push_back(En_HAD);
	 _elEnergyMuon->push_back(En_MUON);
	 if(pandorapart->getClusters().size()==1){
	   _elClusterE->push_back(pandorapart->getClusters()[0]->getEnergy());
	   //sudetector 0:ECAL
	   _elECAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]);
	   //sudetector 1:HCAL
	   _elHCAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]);
	   _elEoverTot->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapart->getClusters()[0]->getEnergy());
	   Track* track=pandorapart->getTracks()[0] ;		
	   PandoraApi::Track::Parameters trackParametersAll;	
	   //std::cout<<"track state electron "<<std::endl;	
	   this->GetTrackStatesOld(track, trackParametersAll);
	   const pandora::CartesianVector &trackPosition=trackParametersAll.m_trackStateAtCalorimeter.Get().GetPosition();
	   //in pandora cluster initialdirection is the unit-vector of the sum of associated calohit position unit vectors
	   //define CartesianVector of the initial direction as untit vector of the  cluster position
	   const float *pClusterPosition(pandorapart->getClusters()[0]->getPosition());
	   //position is NOT normalized but doesn't matter for the opening angle definition
	   const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
	   _elTrackClusterOpeningAngle->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
	   int m_maxSearchLayer=9;
	   float m_parallelDistanceCut=100;
	   float m_minTrackClusterCosAngle=0;
	   float trackClusterDistanceTest=-1;
	   this->GetTrackClusterDistance(trackParametersAll, pandorapart->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest);
	   if(trackClusterDistanceTest==(-1)){
	     std::cout<<"no track cluster distance -> el should HAVE this track as clostest"<<std::endl;
	   }
	   _elMinTrackClustDiff->push_back(trackClusterDistanceTest);
	   _elSigmaPOverPTrack->push_back(std::sqrt(track->getCovMatrix()[5])/std::fabs(track->getOmega()));
	   _elTrackPt->push_back(sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ()));
	   _elTrackP->push_back(trackParametersAll.m_momentumAtDca.Get().GetMagnitude());
	 }else if(pandorapart->getClusters().size()!=1){
	   //std::cout<<"el has not exactly one cluster"<<std::endl;
	   //if(fabs(temp.Theta())>2.4 && temp.Energy()>10){
	   //for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	   //for(unsigned int k=0;k<pandorapart->getClusters()[j]->getSubdetectorEnergies().size();k++){
	   //  std::cout<<" el in cluster["<<j<<" ] of "<<temp.Theta()<<" "<<pandorapart->getClusters().size()<<" E tot "<<pandorapart->getClusters()[j]->getEnergy()<<" E sub "<<k <<" "<<pandorapart->getClusters()[j]->getSubdetectorEnergies()[k]<<" subcluster size? "<<pandorapart->getClusters()[j]->getClusters().size()<<std::endl;
	   //}
	   //}
	 }
	 if(pandorapart->getTracks().size()>0){
	  _el0D0->push_back(pandorapart->getTracks()[0]->getD0());
	  _el0Z0->push_back(pandorapart->getTracks()[0]->getZ0());
	  _el0Phi->push_back(pandorapart->getTracks()[0]->getPhi());
	  _el0Chi2_NDOF->push_back(pandorapart->getTracks()[0]->getChi2()/(float)pandorapart->getTracks()[0]->getNdf());
	  _elOmega->push_back(pandorapart->getTracks()[0]->getOmega());
	  _el1D0->push_back(-1);
	  _el1Z0->push_back(-1);
	  _el1Phi->push_back(-10);
	  _el1Chi2_NDOF->push_back(-1);
	  _el1NHits->push_back(-1);
	  _el1NHitsVeto->push_back(-1);
	  if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()==12){
	    //at the moment 0-5 are NHits used in the fit,6-11 are all hits including those not used
	    //order VTX, FTD, SIT, TPC, SET, ETD
	    _elNHitsVTX->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]);
	    _elNHitsFTD->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]);
	    _elNHitsSIT->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]);
	    _elNHitsTPC->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]);
	    _elNHitsSET->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]);
	    _elNHitsETD->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	    _elNHitsVTXVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[6]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]);
	    _elNHitsFTDVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[7]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]);
	    _elNHitsSITVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[8]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]);
	    _elNHitsTPCVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[9]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]);
	    _elNHitsSETVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[10]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]);
	    _elNHitsETDVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[11]+
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	    //ECAL is 0, HCAL is 1 in calo subdetector stuff
	    if(pandorapart->getTracks().size()>1){
	      (*_el1D0)[_el1D0->size()-1]=pandorapart->getTracks()[1]->getD0();
	      (*_el1Z0)[_el1Z0->size()-1]=pandorapart->getTracks()[1]->getZ0();
	      (*_el1Phi)[_el1Phi->size()-1]=pandorapart->getTracks()[1]->getPhi();
	      (*_el1Chi2_NDOF)[_el1Chi2_NDOF->size()-1]=pandorapart->getTracks()[1]->getChi2()/(float)pandorapart->getTracks()[1]->getNdf();
	      (*_el1NHits)[_el1NHits->size()-1]=pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[0]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[1]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[2]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[3]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[4]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[5];
	      (*_el1NHitsVeto)[_el1NHitsVeto->size()-1]=pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[6]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[7]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[8]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[9]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[10]
		+pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[11]
		-(*_el1NHits)[_el1NHits->size()-1];
	    }
	  }else if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<12){
	    std::cout<<"size of sudetector tracker hits changed, check for new model and implementation"<<pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	  }else if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()>12){
	    std::cout<<"size of sudetector tracker hits changed, check for new model and implementation, maybe 9 very soon ? "<<pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	  }	
	}
       }else if(abs(pandorapart->getType())==13 && pandorapart->getEnergy()>=_lepEMin){
	 //std::cout<<"reco mu"<<" in run/event/weight "<<evt->getRunNumber ()<<"/"<<evt->getEventNumber()<<"/"<<evt->getWeight()<<"/"<<_lepEMin<<"/"<<"/"<<_phEMin<<"/"<<_jetEMin<<std::endl;
	 _muID->push_back(pandorapart->getType());
	 _muE->push_back(pandorapart->getEnergy());
	 _muPx->push_back(pandorapart->getMomentum()[0]);
	 _muPy->push_back(pandorapart->getMomentum()[1]);
	 _muPz->push_back(pandorapart->getMomentum()[2]);
	 _muPhi->push_back(temp.Phi());
	 _muTheta->push_back(temp.Theta());
	 _muMass->push_back(pandorapart->getMass());
	 if(pandorapart->getClusters().size()==1){
	   _muClusterE->push_back(pandorapart->getClusters()[0]->getEnergy());
	   //sudetector 0:ECAL
	   _muECAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]);
	   //sudetector 1:HCAL
	   _muHCAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]);
	   _muEoverTot->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapart->getClusters()[0]->getEnergy());
	 }else if(pandorapart->getClusters().size()!=1){
	   std::cout<<"mu has not exactly one cluster"<<std::endl;
	   for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	     for(unsigned int k=0;k<pandorapart->getClusters()[j]->getSubdetectorEnergies().size();k++){
	       std::cout<<" mu in cluster["<<j<<" ] of "<<temp.Theta()<<" "<<pandorapart->getClusters().size()<<" E tot "<<pandorapart->getClusters()[j]->getEnergy()<<" E sub "<<k <<" "<<pandorapart->getClusters()[j]->getSubdetectorEnergies()[k]<<" subcluster size? "<<pandorapart->getClusters()[j]->getClusters().size()<<std::endl;
	     }
	   }
	 }
	 if(pandorapart->getTracks().size()>0){
	   Track* mutrack=pandorapart->getTracks()[0] ;		
	   PandoraApi::Track::Parameters trackParametersAll;
	   //std::cout<<"track state MU "<<std::endl;
	   this->GetTrackStatesOld(mutrack, trackParametersAll);
	   _muSigmaPOverPTrack->push_back(std::sqrt(mutrack->getCovMatrix()[5])/std::fabs(mutrack->getOmega()));
	   _muTrackPt->push_back(sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ()));
	   _muTrackP->push_back(trackParametersAll.m_momentumAtDca.Get().GetMagnitude());
	   _mu0D0->push_back(pandorapart->getTracks()[0]->getD0());
	   _mu0Z0->push_back(pandorapart->getTracks()[0]->getZ0());
	   _mu0Phi->push_back(pandorapart->getTracks()[0]->getPhi());
	   _mu0Chi2_NDOF->push_back(pandorapart->getTracks()[0]->getChi2()/(float)pandorapart->getTracks()[0]->getNdf());
	   _muOmega->push_back(pandorapart->getTracks()[0]->getOmega());
	   _mu1D0->push_back(-1);
	   _mu1Z0->push_back(-1);
	   _mu1Phi->push_back(-10);
	   _mu1Chi2_NDOF->push_back(-1);
	   _mu1NHits->push_back(-1);
	   _mu1NHitsVeto->push_back(-1);
	   if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()==12){
	     //at the moment 0-5 are NHits used in the fit,6-12 are all hits including those not used
	     //order VTX, FTD, SIT, TPC, SET, ETD
	     _muNHitsVTX->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]);
	     _muNHitsFTD->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]);
	     _muNHitsSIT->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]);
	     _muNHitsTPC->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]);
	     _muNHitsSET->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]);
	     _muNHitsETD->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	     _muNHitsVTXVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[6]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]);
	     _muNHitsFTDVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[7]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]);
	     _muNHitsSITVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[8]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]);
	     _muNHitsTPCVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[9]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]);
	     _muNHitsSETVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[10]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]);
	     _muNHitsETDVeto->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[11]+
					      -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	     //ECAL is 0, HCAL is 1 in calo subdetector stuff
	     if(pandorapart->getTracks().size()>1){
	       (*_mu1D0)[_mu1D0->size()-1]=pandorapart->getTracks()[1]->getD0();
	       (*_mu1Z0)[_mu1Z0->size()-1]=pandorapart->getTracks()[1]->getZ0();
	       (*_mu1Phi)[_mu1Phi->size()-1]=pandorapart->getTracks()[1]->getPhi();
	       (*_mu1Chi2_NDOF)[_mu1Chi2_NDOF->size()-1]=pandorapart->getTracks()[1]->getChi2()/(float)pandorapart->getTracks()[1]->getNdf();
	       (*_mu1NHits)[_mu1NHits->size()-1]=pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[0]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[1]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[2]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[3]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[4]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[5];
	       (*_mu1NHitsVeto)[_mu1NHitsVeto->size()-1]=pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[6]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[7]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[8]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[9]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[10]
		 +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[11]
		 -(*_mu1NHits)[_mu1NHits->size()-1];
	     }
	   }else if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<12){
	     std::cout<<"size of sudetector tracker hits changed, check for new model and implementation "<<pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	   }else if(pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()>12){
	     std::cout<<"size of sudetector tracker hits changed, check for new model and implementation "<<pandorapart->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	   }
	 }
       }else if(pandorapart->getType()==22 && pandorapart->getEnergy()>=_phEMin){
	 //std::cout<<"in reco particle stuff ph"<<std::endl;
	 _phE->push_back(pandorapart->getEnergy());
	 _phPx->push_back(pandorapart->getMomentum()[0]);
	 _phPy->push_back(pandorapart->getMomentum()[1]);
	 _phPz->push_back(pandorapart->getMomentum()[2]);
	 _phPhi->push_back(temp.Phi());
	 _phTheta->push_back(temp.Theta());
	 _phMass->push_back(pandorapart->getMass());
	 float En_EM=0;
	 float En_HAD=0;
	 float En_MUON=0;	  
	 float En_ECAL_Barrel=0;
	 float En_ECAL_Endcap=0;
	 float En_ECAL_else=0;
	 float En_HCAL_Barrel=0;
	 float En_HCAL_Endcap=0;
	 float En_HCAL_else=0;
	 float En_ELSE_Barrel=0;
	 float En_ELSE_Endcap=0;
	 float En_ELSE_else=0;
	 //*_phLongShowerProfile,*_phTransShowerProfile,*_phTrackBasedTransProfile,*_phPeakEnergy,*_phNRadiationLength;
	 int phFirstLayerHCAL=-1;
	 int phFirstLayerECAL=-1;
	 int phLastLayerHCAL=-1;
	 int phLastLayerECAL=-1;
	 int hitsECAL=0;
	 int hitsHCAL=0;
	 //40 layers decided at the moment, 20+10 in samples
	 std::vector<int>hits_per_layer_ecal(50,0);
	 std::vector<float>energy_per_layer_ecal(50,0);
	 //75 layers at the moment
	 std::vector<int>hits_per_layer_hcal(100,0);
	 std::vector<float>energy_per_layer_hcal(100,0);
	 //std::cout<<"in reco particle stuff ph cluster"<<std::endl;
	 for(unsigned int j=0;j<pandorapart->getClusters().size();j++){
	   for(unsigned int l=0; l<pandorapart->getClusters()[j]->getCalorimeterHits().size();l++){
	     TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
	     tempP.SetXYZT(pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
	     int type=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getType();
	     static const int fCaloType =     1 ;
	     static const int fCaloID   =    10 ;
	     static const int fLayout   =  1000 ;
	     static const int fLayer    = 10000 ;
	     int caloLayer=type/fLayer;
	     int caloLayout=(type-caloLayer*fLayer)/fLayout;
	     int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
	     int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
	     //std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
	     if(caloType==0){//0 em, 1 had, 2 mu
	       En_EM+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloType==1){     
	       En_HAD+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloType==2){
	       En_MUON+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	     if(caloID==1){//ecal
	       if(phFirstLayerECAL<0){
		 phFirstLayerECAL=caloLayer;
	       }else if (caloLayer<phFirstLayerECAL){
		 phFirstLayerECAL=caloLayer;
	       }
	       if(phLastLayerECAL<0){
		 phLastLayerECAL=caloLayer;
	       }else if (caloLayer>phLastLayerECAL){
		 phLastLayerECAL=caloLayer;
	       }
	       hitsECAL+=1;
	       hits_per_layer_ecal[caloLayer]+=1;
	       energy_per_layer_ecal[caloLayer]+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if (caloID==2){
	       if(phFirstLayerHCAL<0){
		 phFirstLayerHCAL=caloLayer;
	       }else if (caloLayer<phFirstLayerHCAL){
		 phFirstLayerHCAL=caloLayer;
	       }
	       if(phLastLayerHCAL<0){
		 phLastLayerHCAL=caloLayer;
	       }else if (caloLayer>phLastLayerHCAL){
		 phLastLayerHCAL=caloLayer;
	       }
	       hitsHCAL+=1;
	       hits_per_layer_hcal[caloLayer]+=1;
	       energy_per_layer_hcal[caloLayer]+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	     if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
	       En_ECAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==1 && caloLayout==2){
	       En_ECAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else  if(caloID==1){
	       En_ECAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==2 && caloLayout==1){
	       En_HCAL_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if(caloID==2 && caloLayout==2){
	       En_HCAL_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else	if(caloID==2){
	       En_HCAL_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if((caloID!=1 && caloID!=2) && caloLayout==1){
	       En_ELSE_Barrel+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else if((caloID!=1 && caloID!=2) && caloLayout==2){
	       En_ELSE_Endcap+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }else	if(caloID!=1 && caloID!=2){
	       En_ELSE_else+=pandorapart->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
	     }
	   }
	 }//loop over all clusters and hits for energy in different detectors
	 //*_phLongShowerProfile,*_phTransShowerProfile,*_phTrackBasedTransProfile,*_phPeakEnergy,*_phNRadiationLength;
	 _phHitsECAL->push_back(hitsECAL);
	 _phFirstLayerECAL->push_back(phFirstLayerECAL);
	 _phLastLayerECAL->push_back(phLastLayerECAL);
	 int maxHitsPerLayer=0;
	 float maxEnergyPerLayer=0;
	 if(hitsECAL>0){
	   _phLayerMaxHitsECAL->push_back(-1);
	   _phLayerMaxEECAL->push_back(-1);
	   for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
	     if(hits_per_layer_ecal[j]>maxHitsPerLayer){
	       maxHitsPerLayer=hits_per_layer_ecal[j];
	       (*_phLayerMaxHitsECAL)[_phLayerMaxHitsECAL->size()-1]=j;
	     }
	     if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
	       maxEnergyPerLayer=energy_per_layer_ecal[j];
	       (*_phLayerMaxEECAL)[_phLayerMaxEECAL->size()-1]=j;
	     }
	   }
	   _phMaxHitsPerLayerECAL->push_back(maxHitsPerLayer);
	   _phMaxEInECALLayer->push_back(maxEnergyPerLayer);
	 }else{//no ECAL hits
	   _phLayerMaxHitsECAL->push_back(-1);
	   _phLayerMaxEECAL->push_back(-1);
	   _phMaxHitsPerLayerECAL->push_back(0);
	   _phMaxEInECALLayer->push_back(0);
	 }
	 _phHitsHCAL->push_back(hitsHCAL);
	 _phFirstLayerHCAL->push_back(phFirstLayerHCAL);
	 _phLastLayerHCAL->push_back(phLastLayerHCAL);
	 //initialize to 0 again for HCAL
	 maxHitsPerLayer=0;
	 maxEnergyPerLayer=0;
	 if(hitsHCAL>0){
	   for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
	     if(hits_per_layer_ecal[j]>maxHitsPerLayer){
	       maxHitsPerLayer=hits_per_layer_ecal[j];
	     }
	     if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
	       maxEnergyPerLayer=energy_per_layer_ecal[j];
	     }
	   }
	   _phMaxHitsPerLayerHCAL->push_back(maxHitsPerLayer);
	   _phMaxEInHCALLayer->push_back(maxEnergyPerLayer);
	 }else{//no HCAL hits
	   _phMaxHitsPerLayerHCAL->push_back(0);
	   _phMaxEInHCALLayer->push_back(0);
	 }
	 _phEnergyEB->push_back(En_ECAL_Barrel);
	 _phEnergyEE->push_back(En_ECAL_Endcap);
	 _phEnergyEElse->push_back(En_ECAL_else);
	 _phEnergyHE->push_back(En_HCAL_Barrel);
	 _phEnergyHB->push_back(En_HCAL_Endcap);
	 _phEnergyHElse->push_back( En_HCAL_else);
	 _phEnergyElseB->push_back(En_ELSE_Barrel);
	 _phEnergyElseE->push_back(En_ELSE_Endcap);
	 _phEnergyElseElse->push_back(En_ELSE_else);
	 _phEnergyEM->push_back(En_EM);
	 _phEnergyHAD->push_back(En_HAD);
	 _phEnergyMuon->push_back(En_MUON);
	 //until here for photonID/fake/cluster shape etc etc
	 if(pandorapart->getClusters().size()==1){
	   _phClusterE->push_back(pandorapart->getClusters()[0]->getEnergy());
	   //sudetector 0:ECAL
	   _phECAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]);
	   //sudetector 1:HCAL
	   _phHCAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]);
	   _phEoverTot->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapart->getClusters()[0]->getEnergy());
	   //fraction of the first subcluster -> 1 cluster, obviously 100 %
	   _phf1ClusterE->push_back(1.);
	   _phf1ECAL->push_back(1.);
	   _phf1HCAL->push_back(1.);
	   //fraction of the first subcluster -> 1 cluster, obviously 100 %
	   _phf1ClusterE->push_back(1.);
	   _phf1ECAL->push_back(1.);
	   _phf1HCAL->push_back(1.);
	   //now for unconverted photon check for shower track cluster difference (else we KNOW about the track anyway
	   //only consider the track here, to check for track veto etc
	   float minDistance = std::numeric_limits<float>::max();
	   float minEnergyDifference(std::numeric_limits<float>::max());
	   float sigmaPOverPClosestTrack=-1;
	   float PClosestTrack=-1;
	   float PtClosestTrack=-1;
	   int ind_cTrk=-1;
	   //std::cout<<"reco ph tracks"<<std::endl;
	   //std::cout<<"in reco particle stuff ph track"<<std::endl;
	   for(int t=0;t<trackCol->getNumberOfElements();t++){
	     //std::cout<<"inside ph tracks"<<std::endl;
	     Track* track=dynamic_cast<Track*>(trackCol->getElementAt(t));		
	     PandoraApi::Track::Parameters trackParametersAll;		
	     bool trackstatefails=false;
	     try{
	       this->GetTrackStatesOld(track, trackParametersAll);
	     }
	     catch (pandora::StatusCodeException &statusCodeException)
	       {
		 trackstatefails=true;
		 streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
	       }
	     int m_maxSearchLayer=9;
	     float m_parallelDistanceCut=100;
	     float m_minTrackClusterCosAngle=0;
	     float trackClusterDistanceTest=-1;
	     if(!trackstatefails){
	       this->GetTrackClusterDistance(trackParametersAll, pandorapart->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest);
	       if(trackClusterDistanceTest!=(-1)){
		 float trackparticleMass=trackParametersAll.m_mass.Get();
		 float trackenergy=std::sqrt(trackparticleMass*trackparticleMass + trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared());
		 const float energyDifference=std::fabs(En_HAD - trackenergy);
		 //std::cout<<"energy difference"<<std::endl;
		 if ((trackClusterDistanceTest < minDistance) || ((trackClusterDistanceTest == minDistance) && (energyDifference < minEnergyDifference))){
		   minDistance = trackClusterDistanceTest;
		   minEnergyDifference = energyDifference;
		   //std::cout<<"cov matrix"<<std::endl;
		   sigmaPOverPClosestTrack=std::sqrt(track->getCovMatrix()[5]) / std::fabs(track->getOmega());
		   PClosestTrack=sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ());
		   PtClosestTrack=trackParametersAll.m_momentumAtDca.Get().GetMagnitude();
		   //std::cout<<"DCAstuff"<<std::endl;
		   ind_cTrk=t;
		 }
	       }
	     }
	   }
	   if(PClosestTrack<0){
	     //i.e. no track found
	     minDistance=-1;
	   }
	   //std::cout<<"filling now"<<std::endl;
	   _phSigmaPOverPClosestTrack->push_back(sigmaPOverPClosestTrack);
	   _phClosestTrackP->push_back(PClosestTrack);
	   _phClosestTrackPt->push_back(PtClosestTrack);
	   _phMinTrackClustDiff->push_back(minDistance);
	   //std::cout<<"reco ph sel tracks"<<std::endl;
	   if(ind_cTrk>=0){//closest track was found
	     Track* cltrack=dynamic_cast<Track*>(trackCol->getElementAt(ind_cTrk)) ;	
	     PandoraApi::Track::Parameters trackParametersCT;	
	     //std::cout<<"track staTE PH CLOSEST AGAIN "<<std::endl;
	     this->GetTrackStatesOld(cltrack, trackParametersCT);
	     const pandora::CartesianVector &trackPosition=trackParametersCT.m_trackStateAtCalorimeter.Get().GetPosition();
	     const float *pClusterPosition(pandorapart->getClusters()[0]->getPosition());
	     const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
	     _phClosestTrackClusterOpeningAngle->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
	     //at the moment 0-5 are NHits used in the fit,6-11 are all hits including those not used
	     //order VTX, FTD, SIT, TPC, SET, ETD
	     _phClosestTrackD0->push_back(cltrack->getD0());
	     _phClosestTrackZ0->push_back(cltrack->getZ0());
	     _phClosestTrackPhi->push_back(cltrack->getPhi());
	     _phClosestTrackChi2_NDOF->push_back(cltrack->getChi2()/(float)cltrack->getNdf());
	     _phClosestTrackOmega->push_back(cltrack->getOmega());
	     _phClosestTrackNHitsVTX->push_back(cltrack->getSubdetectorHitNumbers()[0]);
	     _phClosestTrackNHitsFTD->push_back(cltrack->getSubdetectorHitNumbers()[1]);
	     _phClosestTrackNHitsSIT->push_back(cltrack->getSubdetectorHitNumbers()[2]);
	     _phClosestTrackNHitsTPC->push_back(cltrack->getSubdetectorHitNumbers()[3]);
	     _phClosestTrackNHitsSET->push_back(cltrack->getSubdetectorHitNumbers()[4]);
	     _phClosestTrackNHitsETD->push_back(cltrack->getSubdetectorHitNumbers()[5]);
	     _phClosestTrackNHitsVTXVeto->push_back(cltrack->getSubdetectorHitNumbers()[6] - cltrack->getSubdetectorHitNumbers()[0]);
	     _phClosestTrackNHitsFTDVeto->push_back(cltrack->getSubdetectorHitNumbers()[7] - cltrack->getSubdetectorHitNumbers()[1]);
	     _phClosestTrackNHitsSITVeto->push_back(cltrack->getSubdetectorHitNumbers()[8] - cltrack->getSubdetectorHitNumbers()[2]);
	     _phClosestTrackNHitsTPCVeto->push_back(cltrack->getSubdetectorHitNumbers()[9] - cltrack->getSubdetectorHitNumbers()[3]);
	     _phClosestTrackNHitsSETVeto->push_back(cltrack->getSubdetectorHitNumbers()[10] - cltrack->getSubdetectorHitNumbers()[4]);
	     _phClosestTrackNHitsETDVeto->push_back(cltrack->getSubdetectorHitNumbers()[11] - cltrack->getSubdetectorHitNumbers()[5]);
	   }else{
	     _phClosestTrackClusterOpeningAngle->push_back(-1);
	     _phClosestTrackD0->push_back(-1);
	     _phClosestTrackZ0->push_back(-1);
	     _phClosestTrackPhi->push_back(-1);
	     _phClosestTrackChi2_NDOF->push_back(-1);
	     _phClosestTrackOmega->push_back(-1);
	     _phClosestTrackNHitsVTX->push_back(-1);
	     _phClosestTrackNHitsFTD->push_back(-1);
	     _phClosestTrackNHitsSIT->push_back(-1);
	     _phClosestTrackNHitsTPC->push_back(-1);
	     _phClosestTrackNHitsSET->push_back(-1);
	     _phClosestTrackNHitsETD->push_back(-1);
	     _phClosestTrackNHitsVTXVeto->push_back(-1);
	     _phClosestTrackNHitsFTDVeto->push_back(-1);
	     _phClosestTrackNHitsSITVeto->push_back(-1);
	     _phClosestTrackNHitsTPCVeto->push_back(-1);
	     _phClosestTrackNHitsSETVeto->push_back(-1);
	     _phClosestTrackNHitsETDVeto->push_back(-1);
	   }
	   float minDistanceSelTrack = std::numeric_limits<float>::max();
	   float minEnergyDifferenceSelTrack(std::numeric_limits<float>::max());
	   float sigmaPOverPClosestSelTrack=-1;
	   float PClosestSelTrack=-1;
	   float PtClosestSelTrack=-1;
	   int ind_selcTrk =-1;
	   //std::cout<<"reco ph seltracks"<<std::endl;
	   for(int t=0;t<seltrackCol->getNumberOfElements();t++){
	     Track* seltrack=dynamic_cast<Track*>(seltrackCol->getElementAt(t)) ;		
	     PandoraApi::Track::Parameters seltrackParametersAll;	
	     bool seltrackstatefails=false;
	     try{
	       this->GetTrackStatesOld(seltrack, seltrackParametersAll);
	     }
	     catch (pandora::StatusCodeException &statusCodeException)
	       {
		 seltrackstatefails=true;
		 streamlog_out(ERROR) << "Failed to extract a seltrack: " << statusCodeException.ToString() << std::endl;
	       }
	     if(!seltrackstatefails){
	       int m_maxSearchLayerSelTrack=9;
	       float m_parallelDistanceCutSelTrack=100;
	       float m_minSelTrackClusterCosAngle=0;
	       float seltrackClusterDistanceTest=-1;
	       this->GetTrackClusterDistance(seltrackParametersAll, pandorapart->getClusters()[0], m_maxSearchLayerSelTrack, m_parallelDistanceCutSelTrack,m_minSelTrackClusterCosAngle, seltrackClusterDistanceTest);
	       if(seltrackClusterDistanceTest!=(-1)){
		 float seltrackparticleMass=seltrackParametersAll.m_mass.Get();
		 float seltrackenergy=std::sqrt(seltrackparticleMass*seltrackparticleMass + seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared());
		 const float selenergyDifference=std::fabs(En_HAD - seltrackenergy);
		 if ((seltrackClusterDistanceTest < minDistance) || ((seltrackClusterDistanceTest == minDistance) && (selenergyDifference < minEnergyDifferenceSelTrack))){
		   minDistanceSelTrack = seltrackClusterDistanceTest;
		   minEnergyDifferenceSelTrack = selenergyDifference;
		   sigmaPOverPClosestSelTrack=std::sqrt(seltrack->getCovMatrix()[5]) / std::fabs(seltrack->getOmega());
		   PClosestSelTrack=sqrt(seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-seltrackParametersAll.m_momentumAtDca.Get().GetZ()*seltrackParametersAll.m_momentumAtDca.Get().GetZ());
		   PtClosestSelTrack=seltrackParametersAll.m_momentumAtDca.Get().GetMagnitude();
		   ind_selcTrk=t;
		 }
	       }
	     }
	   }
	   if(ind_selcTrk>=0){//closest track was found
	     Track* clseltrack=dynamic_cast<Track*>(seltrackCol->getElementAt(ind_selcTrk)) ;			
	     PandoraApi::Track::Parameters trackParametersCT;	
	     //std::cout<<"track state ph sel track cloests rewind"<<std::endl;
	     this->GetTrackStatesOld(clseltrack, trackParametersCT);
	     const pandora::CartesianVector &trackPosition=trackParametersCT.m_trackStateAtCalorimeter.Get().GetPosition();
	     const float *pClusterPosition(pandorapart->getClusters()[0]->getPosition());
	     const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
	     _phClosestSelTrackClusterOpeningAngle->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
	   }else{
	     _phClosestSelTrackClusterOpeningAngle->push_back(-1);
	   }
	   if(PClosestSelTrack<0){
	     //i.e. no track found
	     minDistanceSelTrack=-1;
	   }
	   _phSigmaPOverPClosestSelTrack->push_back(sigmaPOverPClosestSelTrack);
	   _phClosestSelTrackP->push_back(PClosestSelTrack);
	   _phClosestSelTrackPt->push_back(PtClosestSelTrack);
	   _phMinSelTrackClustDiff->push_back(minDistanceSelTrack);
	 }
	 if(pandorapart->getTracks().size()==0){
	   _phT1Chi2_NDOF  ->push_back(-1.);
	   _phT1Phi        ->push_back(-10.);
	   _phT1NHitsFitted->push_back(0);
	   _phT1NHitsVeto  ->push_back(0);
	   _phT2Chi2_NDOF  ->push_back(-1.);
	   _phT2Phi        ->push_back(-10.);
	   _phT2NHitsFitted->push_back(0);
	   _phT2NHitsVeto  ->push_back(0);
	 }else if(pandorapart->getTracks().size()==2){
	   //converted photon: track values play a role now
	   _phT1Chi2_NDOF  ->push_back(pandorapart->getTracks()[0]->getChi2()/(float)pandorapart->getTracks()[0]->getNdf());
	   _phT1Phi        ->push_back(pandorapart->getTracks()[0]->getPhi());
	   _phT1NHitsFitted->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]
				       +pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]
				       +pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]
				       +pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]
				       +pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]
				       +pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	   _phT1NHitsVeto  ->push_back(pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[6]+
				       pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[7]+
				       pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[8]+
				       pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[9]+
				       pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[10]+
				       pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[11]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[0]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[1]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[2]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[3]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[4]
				       -pandorapart->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	   _phT2Chi2_NDOF  ->push_back(pandorapart->getTracks()[1]->getChi2()/(float)pandorapart->getTracks()[1]->getNdf());
	   _phT2Phi        ->push_back(pandorapart->getTracks()[1]->getPhi());
	   _phT2NHitsFitted->push_back(pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[0]
				       +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[1]
				       +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[2]
				       +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[3]
				       +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[4]
				       +pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	   _phT2NHitsVeto  ->push_back(pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[6]+
				       pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[7]+
				       pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[8]+
				       pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[9]+
				       pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[10]+
				       pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[11]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[0]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[1]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[2]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[3]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[4]
				       -pandorapart->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	   if(pandorapart->getClusters().size()==0){
	     _phClusterE->push_back(0.);
	     //sudetector 0:ECAL
	     _phECAL->push_back(0.);
	     //sudetector 1:HCAL
	     _phHCAL->push_back(0.);
	     _phEoverTot->push_back(-1.);
	     //fraction of the first subcluster -> 0 cluster, obviously 0
	     _phf1ClusterE->push_back(0.);
	     _phf1ECAL->push_back(0.);
	     _phf1HCAL->push_back(0.);
	   }else if(pandorapart->getClusters().size()==2){
	     _phClusterE->push_back(pandorapart->getClusters()[0]->getEnergy()+pandorapart->getClusters()[1]->getEnergy());
	     //sudetector 0:ECAL
	     _phECAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapart->getClusters()[1]->getSubdetectorEnergies()[0]);
	     //sudetector 1:HCAL
	     _phHCAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapart->getClusters()[1]->getSubdetectorEnergies()[1]);
	     _phEoverTot->push_back((pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapart->getClusters()[1]->getSubdetectorEnergies()[0])/(pandorapart->getClusters()[0]->getEnergy()+pandorapart->getClusters()[1]->getEnergy()));
	     //fraction of the first subcluster -> 2 clusters, obviously NOT 100 % 
	     _phf1ClusterE->push_back(pandorapart->getClusters()[0]->getEnergy()/(pandorapart->getClusters()[0]->getEnergy()+pandorapart->getClusters()[1]->getEnergy()));
	     _phf1ECAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]/(pandorapart->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapart->getClusters()[1]->getSubdetectorEnergies()[0]));
	     _phf1HCAL->push_back(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]/(pandorapart->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapart->getClusters()[1]->getSubdetectorEnergies()[1]));
	   }else if (pandorapart->getClusters().size()!=1){
	     std::cout<<"converted photon, but unexpected size of clusters (neither 0,1,2) "<<pandorapart->getClusters().size()<<std::endl;
	   }
	 }else{
	   std::cout<<"photon unusual track number: neither 0 (unconverted) nor 2 (conversion) "<<pandorapart->getTracks().size()<<std::endl;
	 }
       }
       //for isolation of interesting particles (i.e. photons ==22, electrons | |==11 and muons | |==13
       if((abs(pandorapart->getType())==13 && pandorapart->getEnergy()>=_lepEMin) || (abs(pandorapart->getType())==11 && pandorapart->getEnergy()>=_lepEMin)|| (pandorapart->getType()==22 && pandorapart->getEnergy()>=_phEMin)){
	 float sum_DR03_chargedHadrons=0;
	 float sum_DR03_neutralHadrons=0;
	 float sum_DR03_photons=0;
	 float sum_DR03_muons=0;
	 float sum_DR03_electrons=0;
	 float sum_DR04_chargedHadrons=0;
	 float sum_DR04_neutralHadrons=0;
	 float sum_DR04_photons=0;
	 float sum_DR04_muons=0;
	 float sum_DR04_electrons=0;
	 float sum_cosTheta995_chargedHadrons=0;
	 float sum_cosTheta995_neutralHadrons=0;
	 float sum_cosTheta995_photons=0;
	 float sum_cosTheta995_muons=0;
	 float sum_cosTheta995_electrons=0;
	 for(int j=0;j<recoPartCol->getNumberOfElements();j++){
	   if(i!=j){//exclude particle from its own sum
	     ReconstructedParticle* pandorapart2 = dynamic_cast<ReconstructedParticle*>(recoPartCol->getElementAt(j));
	     TLorentzVector temp2;
	     temp2.SetPxPyPzE(pandorapart2->getMomentum()[0],pandorapart2->getMomentum()[1],pandorapart2->getMomentum()[2],pandorapart2->getEnergy());
	     double deltaR=temp2.DeltaR(temp);
	     double cosTheta=(temp.Px()*temp2.Px()+temp.Py()*temp2.Py()+temp.Pz()*temp2.Pz())/(temp.P()*temp2.P());
	     if(cosTheta>0.995){
	       if(abs(pandorapart2->getType())==211){
		 sum_cosTheta995_chargedHadrons+=pandorapart2->getEnergy();
	       }else if(pandorapart2->getType()==22){
		 sum_cosTheta995_photons+=pandorapart2->getEnergy();
	       }else if(abs(pandorapart2->getType())==11){
		 sum_cosTheta995_electrons+=pandorapart2->getEnergy();
	       }else if(abs(pandorapart2->getType())==13){
		 sum_cosTheta995_muons+=pandorapart2->getEnergy();
	       }else if(abs(pandorapart2->getType())>300){
		 sum_cosTheta995_neutralHadrons+=pandorapart2->getEnergy();
	       }
	     }
	     if(deltaR<0.4){
	       if(abs(pandorapart2->getType())==211){
		 sum_DR04_chargedHadrons+=pandorapart2->getEnergy();
		 if(deltaR<0.3){
		   sum_DR03_chargedHadrons+=pandorapart2->getEnergy();
		 }	      
	       }else if(pandorapart2->getType()==22){
		 sum_DR04_photons+=pandorapart2->getEnergy();
		 if(deltaR<0.3){
		   sum_DR03_photons+=pandorapart2->getEnergy();
		 }	      
	       }else if(abs(pandorapart2->getType()==11)){
		 sum_DR04_electrons+=pandorapart2->getEnergy();
		 if(deltaR<0.3){
		   sum_DR03_electrons+=pandorapart2->getEnergy();
		 }
	       }	      
	       else if(abs(pandorapart2->getType()==13)){
		 sum_DR04_muons+=pandorapart2->getEnergy();
		 if(deltaR<0.3){
		   sum_DR03_muons+=pandorapart2->getEnergy();
		 }	      
	       }else if(abs(pandorapart2->getType())>300){
		 //so far neutral I found are 310 for K_short, 2112 for neutrons and 3122 for Lambdas --> should include all neutrals
		 //obviously larger than leptons, charged pions and photons
		 sum_DR04_neutralHadrons+=pandorapart2->getEnergy();
		 if(deltaR<0.3){
		   sum_DR03_neutralHadrons+=pandorapart2->getEnergy();
		 }
	       }
	     }
	   }
	 }//loop over second recopart collection for isolation calculation
	 if(abs(pandorapart->getType())==11){
	   _elCHIsoDR03->push_back(sum_DR03_chargedHadrons);
	   _elCHIsoDR04->push_back(sum_DR04_chargedHadrons);
	   _elNHIsoDR03->push_back(sum_DR03_neutralHadrons);
	   _elNHIsoDR04->push_back(sum_DR04_neutralHadrons);
	   _elPhIsoDR03->push_back(sum_DR03_photons);
	   _elPhIsoDR04->push_back(sum_DR04_photons);
	   _elMuIsoDR03->push_back(sum_DR03_muons);
	   _elMuIsoDR04->push_back(sum_DR04_muons);
	   _elElIsoDR03->push_back(sum_DR03_electrons);
	   _elElIsoDR04->push_back(sum_DR04_electrons);
	   _elCHIsoCosTheta995->push_back(sum_cosTheta995_chargedHadrons);
	   _elNHIsoCosTheta995->push_back(sum_cosTheta995_neutralHadrons);
	   _elPhIsoCosTheta995->push_back(sum_cosTheta995_photons);
	   _elMuIsoCosTheta995->push_back(sum_cosTheta995_muons);
	   _elElIsoCosTheta995->push_back(sum_cosTheta995_electrons);
	 }else if(abs(pandorapart->getType())==13){
	   _muCHIsoDR03->push_back(sum_DR03_chargedHadrons);
	   _muCHIsoDR04->push_back(sum_DR04_chargedHadrons);
	   _muNHIsoDR03->push_back(sum_DR03_neutralHadrons);
	   _muNHIsoDR04->push_back(sum_DR04_neutralHadrons);
	   _muPhIsoDR03->push_back(sum_DR03_photons);
	   _muPhIsoDR04->push_back(sum_DR04_photons);
	   _muMuIsoDR03->push_back(sum_DR03_muons);
	   _muMuIsoDR04->push_back(sum_DR04_muons);
	   _muElIsoDR03->push_back(sum_DR03_electrons);
	   _muElIsoDR04->push_back(sum_DR04_electrons);
	   _muCHIsoCosTheta995->push_back(sum_cosTheta995_chargedHadrons);
	   _muNHIsoCosTheta995->push_back(sum_cosTheta995_neutralHadrons);
	   _muPhIsoCosTheta995->push_back(sum_cosTheta995_photons);
	   _muMuIsoCosTheta995->push_back(sum_cosTheta995_muons);
	   _muElIsoCosTheta995->push_back(sum_cosTheta995_electrons);
	 }else if(pandorapart->getType()==22){
	   _phCHIsoDR03->push_back(sum_DR03_chargedHadrons);
	   _phCHIsoDR04->push_back(sum_DR04_chargedHadrons);
	   _phNHIsoDR03->push_back(sum_DR03_neutralHadrons);
	   _phNHIsoDR04->push_back(sum_DR04_neutralHadrons);
	   _phPhIsoDR03->push_back(sum_DR03_photons);
	   _phPhIsoDR04->push_back(sum_DR04_photons);
	   _phMuIsoDR03->push_back(sum_DR03_muons);
	   _phMuIsoDR04->push_back(sum_DR04_muons);
	   _phElIsoDR03->push_back(sum_DR03_electrons);
	   _phElIsoDR04->push_back(sum_DR04_electrons);
	   _phCHIsoCosTheta995->push_back(sum_cosTheta995_chargedHadrons);
	   _phNHIsoCosTheta995->push_back(sum_cosTheta995_neutralHadrons);
	   _phPhIsoCosTheta995->push_back(sum_cosTheta995_photons);
	   _phMuIsoCosTheta995->push_back(sum_cosTheta995_muons);
	   _phElIsoCosTheta995->push_back(sum_cosTheta995_electrons);
	 }
       }//recopart is either muon, electron or photon
     }
    }

    LCCollection* jetCol = NULL;
    try{
      evt->getCollection( _jetColName ) ;
    }catch( lcio::DataNotAvailableException e ){
      std::cout<<"in run/event/weight no jets "<<evt->getRunNumber ()<<"/"<<evt->getEventNumber()<<"/"<<evt->getWeight()<<"/"<<_lepEMin<<"/"<<"/"<<_phEMin<<"/"<<_jetEMin<<std::endl;
    }
    if(jetCol!=NULL){
      //std::cout<<"reco jet col"<<std::endl;
      for(int i=0;i<jetCol->getNumberOfElements();i++){
	ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(jetCol->getElementAt(i));
	if(recojet->getEnergy()>_jetEMin){
	  _jetE->push_back(recojet->getEnergy());
	  _jetPx->push_back(recojet->getMomentum()[0]);
	  _jetPy->push_back(recojet->getMomentum()[1]);
	  _jetPz->push_back(recojet->getMomentum()[2]);
	  TLorentzVector temp;
	  temp.SetPxPyPzE(recojet->getMomentum()[0],recojet->getMomentum()[1],recojet->getMomentum()[2],recojet->getEnergy());
	  _jetPhi->push_back(temp.Phi());
	  _jetTheta->push_back(temp.Theta());
	  _jetMass->push_back(recojet->getMass());
	  _jetNF->push_back(0);
	  _jetLF->push_back(0);
	  _jetKF ->push_back(0);
	  _jetCHF->push_back(0);
	  _jetPhF->push_back(0);
	  _jetElF->push_back(0);
	  _jetMuF->push_back(0);
	  _jetNMult->push_back(0);
	  _jetLMult->push_back(0);
	  _jetKMult->push_back(0);
	  _jetCHMult->push_back(0);
	  _jetPhMult->push_back(0);
	  _jetElMult->push_back(0);
	  _jetMuMult->push_back(0);
	  unsigned int jet_index=_jetPhi->size()-1;
	  for(unsigned int j=0;j<recojet->getParticles().size();j++){
	    if(abs(recojet->getParticles()[j]->getType())==211){ 
	      (*_jetCHF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetCHMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==22){
	      (*_jetPhF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetPhMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==11){
	      (*_jetElF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetElMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==13){
	      (*_jetMuF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetMuMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==2112){
	      (*_jetNF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetNMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==3122){
	      (*_jetLF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	      (*_jetLMult)[jet_index]+=1;
	    }else if(abs(recojet->getParticles()[j]->getType())==310){
	      (*_jetKF)[jet_index]+=recojet->getParticles()[j]->getEnergy();
	    (*_jetKMult)[jet_index]+=1;
	    }else{
	      std::cout<<"don't take into account particle type "<<recojet->getParticles()[j]->getType()<<std::endl;
	    }
	  }
	  float tot_energy=(*_jetCHF)[jet_index]+(*_jetPhF)[jet_index]+(*_jetMuF)[jet_index]+(*_jetElF)[jet_index]+(*_jetNF)[jet_index]+(*_jetLF)[jet_index]+(*_jetKF)[jet_index];
	  (*_jetCHF)[jet_index]=(*_jetCHF)[jet_index]/recojet->getEnergy();
	  (*_jetPhF)[jet_index]=(*_jetPhF)[jet_index]/recojet->getEnergy();
	  (*_jetMuF)[jet_index]=(*_jetMuF)[jet_index]/recojet->getEnergy();
	  (*_jetElF)[jet_index]=(*_jetElF)[jet_index]/recojet->getEnergy();
	  (*_jetNF)[jet_index]=(*_jetNF)[jet_index]/recojet->getEnergy();
	  (*_jetLF)[jet_index]=(*_jetLF)[jet_index]/recojet->getEnergy();
	  (*_jetKF)[jet_index]=(*_jetKF)[jet_index]/recojet->getEnergy();
	  float tot_fraction=(*_jetCHF)[jet_index]+(*_jetPhF)[jet_index]+(*_jetMuF)[jet_index]+(*_jetElF)[jet_index]+(*_jetNF)[jet_index]+(*_jetLF)[jet_index]+(*_jetKF)[jet_index];
	  if(fabs(1.-(*_jetCHF)[jet_index]-(*_jetPhF)[jet_index]-(*_jetMuF)[jet_index]-(*_jetElF)[jet_index]-(*_jetNF)[jet_index]-(*_jetLF)[jet_index]-(*_jetKF)[jet_index])>1.e-4){
	    std::cout<<"big difference in fractions "<<fabs(1.-(*_jetCHF)[jet_index]-(*_jetPhF)[jet_index]-(*_jetMuF)[jet_index]-(*_jetElF)[jet_index]-(*_jetNF)[jet_index]-(*_jetLF)[jet_index]-(*_jetKF)[jet_index])<<"/"<<tot_energy<<"/"<<recojet->getEnergy()<<"/"<<tot_energy/recojet->getEnergy()<<"/"<<tot_fraction<<std::endl;
	  }
	}
      }
    }


    _MET=sqrt(pow(_MEx,2)+pow(_MEy,2));

    
    if(_hasrereco_information){
      //std::cout<<"rereco"<<std::endl;
      LCCollection* recoPartColReReco = evt->getCollection( _recoPartNameReReco ) ;
      _MExReReco=0;_MEyReReco=0;_MEzReReco=0;_METReReco=0;_SumEReReco=0;_SumPtReReco=0;_SumChargedHadronEReReco=0;_SumPhotonEReReco=0;_SumNeutralHadronEReReco=0;_SumChargedHadronPtReReco=0;_SumPhotonPtReReco=0;_SumNeutralHadronPtReReco=0;_SumMuonEReReco=0;_SumMuonPtReReco=0;_SumElectronEReReco=0;_SumElectronPtReReco=0;
      if(recoPartColReReco!=NULL){
	//PandoraCandidate loop
	for(int i=0;i<recoPartColReReco->getNumberOfElements();i++){
	  ReconstructedParticle* pandorapartReReco = dynamic_cast<ReconstructedParticle*>(recoPartColReReco->getElementAt(i));
	  _MExReReco-=pandorapartReReco->getMomentum()[0];
	  _MEyReReco-=pandorapartReReco->getMomentum()[1];
	  _MEzReReco-=pandorapartReReco->getMomentum()[2];
	  _SumEReReco+=pandorapartReReco->getEnergy();
	  _SumPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  TLorentzVector temp;
	  temp.SetPxPyPzE(pandorapartReReco->getMomentum()[0],pandorapartReReco->getMomentum()[1],pandorapartReReco->getMomentum()[2],pandorapartReReco->getEnergy());
	  if(abs(pandorapartReReco->getType())==211){
	    _SumChargedHadronEReReco+=pandorapartReReco->getEnergy();
	    _SumChargedHadronPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  }else if(pandorapartReReco->getType()==22){
	    _SumPhotonEReReco+=pandorapartReReco->getEnergy();
	    _SumPhotonPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  }else if(abs(pandorapartReReco->getType())==11){
	    _SumElectronEReReco+=pandorapartReReco->getEnergy();
	    _SumElectronPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  }else if(abs(pandorapartReReco->getType())==13){
	    _SumMuonEReReco+=pandorapartReReco->getEnergy();
	    _SumMuonPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  }else{
	    _SumNeutralHadronEReReco+=pandorapartReReco->getEnergy();
	    _SumNeutralHadronPtReReco+=sqrt(pow(pandorapartReReco->getMomentum()[1],2)+pow(pandorapartReReco->getMomentum()[2],2));
	  }
	  if(abs(pandorapartReReco->getType())==11 && pandorapartReReco->getEnergy()>=_lepEMin){
	    //std::cout<<"rereco e"<<std::endl;
	    _elIDReReco->push_back(pandorapartReReco->getType());
	    _elEReReco->push_back(pandorapartReReco->getEnergy());
	    _elPxReReco->push_back(pandorapartReReco->getMomentum()[0]);
	    _elPyReReco->push_back(pandorapartReReco->getMomentum()[1]);
	    _elPzReReco->push_back(pandorapartReReco->getMomentum()[2]);
	    _elPhiReReco->push_back(temp.Phi());
	    _elThetaReReco->push_back(temp.Theta());
	    _elMassReReco->push_back(pandorapartReReco->getMass());
	    //start here eloton specific cluster/composition/shower shape etc etc
	    float En_EM=0;
	    float En_HAD=0;
	    float En_MUON=0;	  
	    float En_ECAL_Barrel=0;
	    float En_ECAL_Endcap=0;
	    float En_ECAL_else=0;
	    float En_HCAL_Barrel=0;
	    float En_HCAL_Endcap=0;
	    float En_HCAL_else=0;
	    float En_ELSE_Barrel=0;
	    float En_ELSE_Endcap=0;
	    float En_ELSE_else=0;
	    //*_elLongShowerProfileReReco,*_elTransShowerProfileReReco,*_elTrackBasedTransProfileReReco,*_elPeakEnergyReReco,*_elNRadiationLengthReReco;
	    int elFirstLayerHCAL=-1;
	    int elFirstLayerECAL=-1;
	    int elLastLayerHCAL=-1;
	    int elLastLayerECAL=-1;
	    int hitsECAL=0;
	    int hitsHCAL=0;
	    //40 layers decided at the moment, 20+10 in samples
	    std::vector<int>hits_per_layer_ecal(50,0);
	    std::vector<float>energy_per_layer_ecal(50,0);
	    //75 layers at the moment
	    std::vector<int>hits_per_layer_hcal(100,0);
	    std::vector<float>energy_per_layer_hcal(100,0);
	    for(unsigned int j=0;j<pandorapartReReco->getClusters().size();j++){
	      for(unsigned int l=0; l<pandorapartReReco->getClusters()[j]->getCalorimeterHits().size();l++){
		TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
		tempP.SetXYZT(pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
		int type=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getType();
		static const int fCaloType =     1 ;
		static const int fCaloID   =    10 ;
		static const int fLayout   =  1000 ;
		static const int fLayer    = 10000 ;
		int caloLayer=type/fLayer;
		int caloLayout=(type-caloLayer*fLayer)/fLayout;
		int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
		int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
		//std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
		if(caloType==0){//0 em, 1 had, 2 mu
		  En_EM+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==1){     
		  En_HAD+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==2){
		  En_MUON+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1){//ecal
		  if(elFirstLayerECAL<0){
		    elFirstLayerECAL=caloLayer;
		  }else if (caloLayer<elFirstLayerECAL){
		    elFirstLayerECAL=caloLayer;
		  }
		  if(elLastLayerECAL<0){
		    elLastLayerECAL=caloLayer;
		  }else if (caloLayer>elLastLayerECAL){
		    elLastLayerECAL=caloLayer;
		  }
		  hitsECAL+=1;
		  hits_per_layer_ecal[caloLayer]+=1;
		  energy_per_layer_ecal[caloLayer]+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if (caloID==2){
		  if(elFirstLayerHCAL<0){
		    elFirstLayerHCAL=caloLayer;
		  }else if (caloLayer<elFirstLayerHCAL){
		    elFirstLayerHCAL=caloLayer;
		  }
		  if(elLastLayerHCAL<0){
		    elLastLayerHCAL=caloLayer;
		  }else if (caloLayer>elLastLayerHCAL){
		    elLastLayerHCAL=caloLayer;
		  }
		  hitsHCAL+=1;
		  hits_per_layer_hcal[caloLayer]+=1;
		  energy_per_layer_hcal[caloLayer]+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_ECAL_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==1 && caloLayout==2){
		  En_ECAL_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==1){
		  En_ECAL_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==1){
		  En_HCAL_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==2){
		  En_HCAL_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==2){
		  En_HCAL_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==1){
		  En_ELSE_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==2){
		  En_ELSE_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID!=1 && caloID!=2){
		  En_ELSE_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
	      }
	    }//loop over all clusters and hits for energy in different detectors
	    //*_elLongShowerProfileReReco,*_elTransShowerProfileReReco,*_elTrackBasedTransProfileReReco,*_elPeakEnergyReReco,*_elNRadiationLengthReReco;
	    _elHitsECALReReco->push_back(hitsECAL);
	    _elFirstLayerECALReReco->push_back(elFirstLayerECAL);
	    _elLastLayerECALReReco->push_back(elLastLayerECAL);
	    int maxHitsPerLayer=0;
	    float maxEnergyPerLayer=0;
	    if(hitsECAL>0){
	      _elLayerMaxHitsECALReReco->push_back(-1);
	      _elLayerMaxEECALReReco->push_back(-1);
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		  (*_elLayerMaxHitsECALReReco)[_elLayerMaxHitsECALReReco->size()-1]=j;
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		  (*_elLayerMaxEECALReReco)[_elLayerMaxEECALReReco->size()-1]=j;
		}
	      }
	      _elMaxHitsPerLayerECALReReco->push_back(maxHitsPerLayer);
	      _elMaxEInECALLayerReReco->push_back(maxEnergyPerLayer);
	    }else{//no ECAL hits
	      _elLayerMaxHitsECALReReco->push_back(-1);
	      _elLayerMaxEECALReReco->push_back(-1);
	      _elMaxHitsPerLayerECALReReco->push_back(0);
	      _elMaxEInECALLayerReReco->push_back(0);
	    }
	    _elHitsHCALReReco->push_back(hitsHCAL);
	    _elFirstLayerHCALReReco->push_back(elFirstLayerHCAL);
	    _elLastLayerHCALReReco->push_back(elLastLayerHCAL);
	    //initialize to 0 again for HCAL
	    maxHitsPerLayer=0;
	    maxEnergyPerLayer=0;
	    if(hitsHCAL>0){
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		}
	      }
	      _elMaxHitsPerLayerHCALReReco->push_back(maxHitsPerLayer);
	      _elMaxEInHCALLayerReReco->push_back(maxEnergyPerLayer);
	    }else{//no HCAL hits
	      _elLayerMaxHitsECALReReco->push_back(-1);
	      _elMaxHitsPerLayerHCALReReco->push_back(0);
	      _elMaxEInHCALLayerReReco->push_back(0);
	    }
	    _elEnergyEBReReco->push_back(En_ECAL_Barrel);
	    _elEnergyEEReReco->push_back(En_ECAL_Endcap);
	    _elEnergyEElseReReco->push_back(En_ECAL_else);
	    _elEnergyHEReReco->push_back(En_HCAL_Barrel);
	    _elEnergyHBReReco->push_back(En_HCAL_Endcap);
	    _elEnergyHElseReReco->push_back( En_HCAL_else);
	    _elEnergyElseBReReco->push_back(En_ELSE_Barrel);
	    _elEnergyElseEReReco->push_back(En_ELSE_Endcap);
	    _elEnergyElseElseReReco->push_back(En_ELSE_else);
	    _elEnergyEMReReco->push_back(En_EM);
	    _elEnergyHADReReco->push_back(En_HAD);
	    _elEnergyMuonReReco->push_back(En_MUON);
	    Track* eltrack=pandorapartReReco->getTracks()[0];
	    PandoraApi::Track::Parameters trackParametersAll;	
	    //std::cout<<"track state e rereo "<<std::endl;
	    this->GetTrackStatesOld(eltrack, trackParametersAll);
	    _elSigmaPOverPTrackReReco->push_back(std::sqrt(eltrack->getCovMatrix()[5])/std::fabs(eltrack->getOmega()));
	    _elTrackPtReReco->push_back(sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ()));
	    _elTrackPReReco->push_back(trackParametersAll.m_momentumAtDca.Get().GetMagnitude());
	    if(pandorapartReReco->getClusters().size()==1){
	      const pandora::CartesianVector &trackPosition=trackParametersAll.m_trackStateAtCalorimeter.Get().GetPosition();
	      const float *pClusterPosition(pandorapartReReco->getClusters()[0]->getPosition());
	      const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
	      _elTrackClusterOpeningAngleReReco->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
	      _elClusterEReReco->push_back(pandorapartReReco->getClusters()[0]->getEnergy());
	      //sudetector 0:ECAL
	      _elECALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]);
	      //sudetector 1:HCAL
	      _elHCALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]);
	      _elEoverTotReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapartReReco->getClusters()[0]->getEnergy());
	      //unlike for photons for electrons we KNOW about the track anyway
	      //only consider this track here, to check for track veto etc		
	      int m_maxSearchLayer=9;
	      float m_parallelDistanceCut=100;
	      float m_minTrackClusterCosAngle=0;
	      float trackClusterDistanceTest=-1;
	      this->GetTrackClusterDistance(trackParametersAll, pandorapartReReco->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest);
	      if(trackClusterDistanceTest==(-1)){
		std::cout<<"no track cluster distance -> el rereco should HAVE this track as closest"<<std::endl;
	      }
	      _elMinTrackClustDiffReReco->push_back(trackClusterDistanceTest);
	    }else if(pandorapartReReco->getClusters().size()!=1){
	      std::cout<<"el rereco has not exactly one cluster"<<std::endl;
	      //if(fabs(temp.Theta())>2.4 && temp.Energy()>10){
	      for(unsigned int j=0;j<pandorapartReReco->getClusters().size();j++){
		for(unsigned int k=0;k<pandorapartReReco->getClusters()[j]->getSubdetectorEnergies().size();k++){
		  std::cout<<" el rereco in cluster["<<j<<" ] of "<<temp.Theta()<<" "<<pandorapartReReco->getClusters().size()<<" E tot "<<pandorapartReReco->getClusters()[j]->getEnergy()<<" E sub "<<k <<" "<<pandorapartReReco->getClusters()[j]->getSubdetectorEnergies()[k]<<" subcluster size? "<<pandorapartReReco->getClusters()[j]->getClusters().size()<<std::endl;
		}
	      }
	    }
	    if(pandorapartReReco->getTracks().size()>0){
	      _el0D0ReReco->push_back(pandorapartReReco->getTracks()[0]->getD0());
	      _el0Z0ReReco->push_back(pandorapartReReco->getTracks()[0]->getZ0());
	      _el0PhiReReco->push_back(pandorapartReReco->getTracks()[0]->getPhi());
	      _el0Chi2_NDOFReReco->push_back(pandorapartReReco->getTracks()[0]->getChi2()/(float)pandorapartReReco->getTracks()[0]->getNdf());
	      _elOmegaReReco->push_back(pandorapartReReco->getTracks()[0]->getOmega());
	      _el1D0ReReco->push_back(-1);
	      _el1Z0ReReco->push_back(-1);
	      _el1PhiReReco->push_back(-10);
	      _el1Chi2_NDOFReReco->push_back(-1);
	      _el1NHitsReReco->push_back(-1);
	      _el1NHitsVetoReReco->push_back(-1);
	      if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()==12){
		//at the moment 0-5 are NHits used in the fit,6-11 are all hits including those not used
		//order VTX, FTD, SIT, TPC, SET, ETD
		_elNHitsVTXReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_elNHitsFTDReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_elNHitsSITReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_elNHitsTPCReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_elNHitsSETReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_elNHitsETDReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		_elNHitsVTXVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[6]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_elNHitsFTDVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[7]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_elNHitsSITVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[8]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_elNHitsTPCVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[9]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_elNHitsSETVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[10]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_elNHitsETDVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[11]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		//ECAL is 0, HCAL is 1 in calo subdetector stuff
		if(pandorapartReReco->getTracks().size()>1){
		  (*_el1D0ReReco)[_el1D0ReReco->size()-1]=pandorapartReReco->getTracks()[1]->getD0();
		  (*_el1Z0ReReco)[_el1Z0ReReco->size()-1]=pandorapartReReco->getTracks()[1]->getZ0();
		  (*_el1PhiReReco)[_el1PhiReReco->size()-1]=pandorapartReReco->getTracks()[1]->getPhi();
		  (*_el1Chi2_NDOFReReco)[_el1Chi2_NDOFReReco->size()-1]=pandorapartReReco->getTracks()[1]->getChi2()/(float)pandorapartReReco->getTracks()[1]->getNdf();
		  (*_el1NHitsReReco)[_el1NHitsReReco->size()-1]=pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[0]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[1]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[2]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[3]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[4]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[5];
		  (*_el1NHitsVetoReReco)[_el1NHitsVetoReReco->size()-1]=pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[6]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[7]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[8]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[9]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[10]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[11]
		    -(*_el1NHitsReReco)[_el1NHitsReReco->size()-1];
		}
	      }else if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation"<<pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }else if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()>12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation, maybe 9 very soon ? "<<pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }	
	    }
	  }else if(abs(pandorapartReReco->getType())==13 && pandorapartReReco->getEnergy()>=_lepEMin){
	    _muIDReReco->push_back(pandorapartReReco->getType());
	    _muEReReco->push_back(pandorapartReReco->getEnergy());
	    _muPxReReco->push_back(pandorapartReReco->getMomentum()[0]);
	    _muPyReReco->push_back(pandorapartReReco->getMomentum()[1]);
	    _muPzReReco->push_back(pandorapartReReco->getMomentum()[2]);
	    _muPhiReReco->push_back(temp.Phi());
	    _muThetaReReco->push_back(temp.Theta());
	    _muMassReReco->push_back(pandorapartReReco->getMass());
	    if(pandorapartReReco->getClusters().size()==1){
	      _muClusterEReReco->push_back(pandorapartReReco->getClusters()[0]->getEnergy());
	      //sudetector 0:ECAL
	      _muECALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]);
	      //sudetector 1:HCAL
	      _muHCALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]);
	      _muEoverTotReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapartReReco->getClusters()[0]->getEnergy());
	    }else if(pandorapartReReco->getClusters().size()!=1){
	      //if(fabs(temp.Theta())>2.4 && temp.Energy()>10){
	      for(unsigned int j=0;j<pandorapartReReco->getClusters().size();j++){
		for(unsigned int k=0;k<pandorapartReReco->getClusters()[j]->getSubdetectorEnergies().size();k++){
		  std::cout<<" mu in cluster["<<j<<" ] of "<<temp.Theta()<<" "<<pandorapartReReco->getClusters().size()<<" E tot "<<pandorapartReReco->getClusters()[j]->getEnergy()<<" E sub "<<k <<" "<<pandorapartReReco->getClusters()[j]->getSubdetectorEnergies()[k]<<" subcluster size? "<<pandorapartReReco->getClusters()[j]->getClusters().size()<<std::endl;
		}
	      }
	    }
	    if(pandorapartReReco->getTracks().size()>0){
	      Track* mutrack=pandorapartReReco->getTracks()[0];
	      PandoraApi::Track::Parameters trackParametersAll;	
	      //std::cout<<"track state MU rereco "<<std::endl;
	      this->GetTrackStatesOld(mutrack, trackParametersAll);
	      _muSigmaPOverPTrackReReco->push_back(std::sqrt(mutrack->getCovMatrix()[5])/std::fabs(mutrack->getOmega()));
	      _muTrackPtReReco->push_back(sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ()));
	      _muTrackPReReco->push_back(trackParametersAll.m_momentumAtDca.Get().GetMagnitude());
	      _mu0D0ReReco->push_back(pandorapartReReco->getTracks()[0]->getD0());
	      _mu0Z0ReReco->push_back(pandorapartReReco->getTracks()[0]->getZ0());
	      _mu0PhiReReco->push_back(pandorapartReReco->getTracks()[0]->getPhi());
	      _mu0Chi2_NDOFReReco->push_back(pandorapartReReco->getTracks()[0]->getChi2()/(float)pandorapartReReco->getTracks()[0]->getNdf());
	      _muOmegaReReco->push_back(pandorapartReReco->getTracks()[0]->getOmega());
	      _mu1D0ReReco->push_back(-1);
	      _mu1Z0ReReco->push_back(-1);
	      _mu1PhiReReco->push_back(-10);
	      _mu1Chi2_NDOFReReco->push_back(-1);
	      _mu1NHitsReReco->push_back(-1);
	      _mu1NHitsVetoReReco->push_back(-1);
	      if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()==12){
		//at the moment 0-5 are NHits used in the fit,6-12 are all hits including those not used
		//order VTX, FTD, SIT, TPC, SET, ETD
		_muNHitsVTXReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_muNHitsFTDReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_muNHitsSITReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_muNHitsTPCReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_muNHitsSETReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_muNHitsETDReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		_muNHitsVTXVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[6]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_muNHitsFTDVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[7]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_muNHitsSITVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[8]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_muNHitsTPCVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[9]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_muNHitsSETVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[10]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_muNHitsETDVetoReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[11]+
						 -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		//ECAL is 0, HCAL is 1 in calo subdetector stuff
		if(pandorapartReReco->getTracks().size()>1){
		  (*_mu1D0ReReco)[_mu1D0ReReco->size()-1]=pandorapartReReco->getTracks()[1]->getD0();
		  (*_mu1Z0ReReco)[_mu1Z0ReReco->size()-1]=pandorapartReReco->getTracks()[1]->getZ0();
		  (*_mu1PhiReReco)[_mu1PhiReReco->size()-1]=pandorapartReReco->getTracks()[1]->getPhi();
		  (*_mu1Chi2_NDOFReReco)[_mu1Chi2_NDOFReReco->size()-1]=pandorapartReReco->getTracks()[1]->getChi2()/(float)pandorapartReReco->getTracks()[1]->getNdf();
		  (*_mu1NHitsReReco)[_mu1NHitsReReco->size()-1]=pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[0]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[1]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[2]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[3]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[4]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[5];
		  (*_mu1NHitsVetoReReco)[_mu1NHitsVetoReReco->size()-1]=pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[6]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[7]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[8]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[9]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[10]
		    +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[11]
		    -(*_mu1NHitsReReco)[_mu1NHitsReReco->size()-1];
		}
	      }else if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation "<<pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }else if(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()>12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation "<<pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }
	    }
	  }else if(pandorapartReReco->getType()==22 && pandorapartReReco->getEnergy()>=_phEMin){
	    _phEReReco->push_back(pandorapartReReco->getEnergy());
	    _phPxReReco->push_back(pandorapartReReco->getMomentum()[0]);
	    _phPyReReco->push_back(pandorapartReReco->getMomentum()[1]);
	    _phPzReReco->push_back(pandorapartReReco->getMomentum()[2]);
	    _phPhiReReco->push_back(temp.Phi());
	    _phThetaReReco->push_back(temp.Theta());
	    _phMassReReco->push_back(pandorapartReReco->getMass());
	    //start here photon specific cluster/composition/shower shape etc etc
	    float En_EM=0;
	    float En_HAD=0;
	    float En_MUON=0;	  
	    float En_ECAL_Barrel=0;
	    float En_ECAL_Endcap=0;
	    float En_ECAL_else=0;
	    float En_HCAL_Barrel=0;
	    float En_HCAL_Endcap=0;
	    float En_HCAL_else=0;
	    float En_ELSE_Barrel=0;
	    float En_ELSE_Endcap=0;
	    float En_ELSE_else=0;
	    //*_phLongShowerProfileReReco,*_phTransShowerProfileReReco,*_phTrackBasedTransProfileReReco,*_phPeakEnergyReReco,*_phNRadiationLengthReReco;
	    int phFirstLayerHCAL=-1;
	    int phFirstLayerECAL=-1;
	    int phLastLayerHCAL=-1;
	    int phLastLayerECAL=-1;
	    int hitsECAL=0;
	    int hitsHCAL=0;
	    //40 layers decided at the moment, 20+10 in samples
	    std::vector<int>hits_per_layer_ecal(50,0);
	    std::vector<float>energy_per_layer_ecal(50,0);
	    //75 layers at the moment
	    std::vector<int>hits_per_layer_hcal(100,0);
	    std::vector<float>energy_per_layer_hcal(100,0);
	    for(unsigned int j=0;j<pandorapartReReco->getClusters().size();j++){
	      for(unsigned int l=0; l<pandorapartReReco->getClusters()[j]->getCalorimeterHits().size();l++){
		TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
		tempP.SetXYZT(pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
		int type=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getType();
		static const int fCaloType =     1 ;
		static const int fCaloID   =    10 ;
		static const int fLayout   =  1000 ;
		static const int fLayer    = 10000 ;
		int caloLayer=type/fLayer;
		int caloLayout=(type-caloLayer*fLayer)/fLayout;
		int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
		int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
		//std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
		if(caloType==0){//0 em, 1 had, 2 mu
		  En_EM+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==1){     
		  En_HAD+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==2){
		  En_MUON+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1){//ecal
		  if(phFirstLayerECAL<0){
		    phFirstLayerECAL=caloLayer;
		  }else if (caloLayer<phFirstLayerECAL){
		    phFirstLayerECAL=caloLayer;
		  }
		  if(phLastLayerECAL<0){
		    phLastLayerECAL=caloLayer;
		  }else if (caloLayer>phLastLayerECAL){
		    phLastLayerECAL=caloLayer;
		  }
		  hitsECAL+=1;
		  hits_per_layer_ecal[caloLayer]+=1;
		  energy_per_layer_ecal[caloLayer]+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if (caloID==2){
		  if(phFirstLayerHCAL<0){
		    phFirstLayerHCAL=caloLayer;
		  }else if (caloLayer<phFirstLayerHCAL){
		    phFirstLayerHCAL=caloLayer;
		  }
		  if(phLastLayerHCAL<0){
		    phLastLayerHCAL=caloLayer;
		  }else if (caloLayer>phLastLayerHCAL){
		    phLastLayerHCAL=caloLayer;
		  }
		  hitsHCAL+=1;
		  hits_per_layer_hcal[caloLayer]+=1;
		  energy_per_layer_hcal[caloLayer]+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_ECAL_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==1 && caloLayout==2){
		  En_ECAL_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==1){
		  En_ECAL_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==1){
		  En_HCAL_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==2){
		  En_HCAL_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==2){
		  En_HCAL_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==1){
		  En_ELSE_Barrel+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==2){
		  En_ELSE_Endcap+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID!=1 && caloID!=2){
		  En_ELSE_else+=pandorapartReReco->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
	      }
	    }//loop over all clusters and hits for energy in different detectors
	    //*_phLongShowerProfileReReco,*_phTransShowerProfileReReco,*_phTrackBasedTransProfileReReco,*_phPeakEnergyReReco,*_phNRadiationLengthReReco;
	    _phHitsECALReReco->push_back(hitsECAL);
	    _phFirstLayerECALReReco->push_back(phFirstLayerECAL);
	    _phLastLayerECALReReco->push_back(phLastLayerECAL);
	    int maxHitsPerLayer=0;
	    float maxEnergyPerLayer=0;
	    if(hitsECAL>0){
	      _phLayerMaxHitsECALReReco->push_back(-1);
	      _phLayerMaxEECALReReco->push_back(-1);
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		  (*_phLayerMaxHitsECALReReco)[_phLayerMaxHitsECALReReco->size()-1]=j;
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		  (*_phLayerMaxEECALReReco)[_phLayerMaxEECALReReco->size()-1]=j;
		}
	      }
	      _phMaxHitsPerLayerECALReReco->push_back(maxHitsPerLayer);
	      _phMaxEInECALLayerReReco->push_back(maxEnergyPerLayer);
	    }else{//no ECAL hits
	      _phLayerMaxHitsECALReReco->push_back(-1);
	      _phLayerMaxEECALReReco->push_back(-1);
	      _phMaxHitsPerLayerECALReReco->push_back(0);
	      _phMaxEInECALLayerReReco->push_back(0);
	    }
	    _phHitsHCALReReco->push_back(hitsHCAL);
	    _phFirstLayerHCALReReco->push_back(phFirstLayerHCAL);
	    _phLastLayerHCALReReco->push_back(phLastLayerHCAL);
	    //initialize to 0 again for HCAL
	    maxHitsPerLayer=0;
	    maxEnergyPerLayer=0;
	    if(hitsHCAL>0){
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		}
	      }
	      _phMaxHitsPerLayerHCALReReco->push_back(maxHitsPerLayer);
	      _phMaxEInHCALLayerReReco->push_back(maxEnergyPerLayer);
	    }else{//no HCAL hits
	      _phLayerMaxHitsECALReReco->push_back(-1);
	      _phMaxHitsPerLayerHCALReReco->push_back(0);
	      _phMaxEInHCALLayerReReco->push_back(0);
	    }
	    _phEnergyEBReReco->push_back(En_ECAL_Barrel);
	    _phEnergyEEReReco->push_back(En_ECAL_Endcap);
	    _phEnergyEElseReReco->push_back(En_ECAL_else);
	    _phEnergyHEReReco->push_back(En_HCAL_Barrel);
	    _phEnergyHBReReco->push_back(En_HCAL_Endcap);
	    _phEnergyHElseReReco->push_back( En_HCAL_else);
	    _phEnergyElseBReReco->push_back(En_ELSE_Barrel);
	    _phEnergyElseEReReco->push_back(En_ELSE_Endcap);
	    _phEnergyElseElseReReco->push_back(En_ELSE_else);
	    _phEnergyEMReReco->push_back(En_EM);
	    _phEnergyHADReReco->push_back(En_HAD);
	    _phEnergyMuonReReco->push_back(En_MUON);
	    //until here for photonID/fake/cluster shape etc etc
	    if(pandorapartReReco->getClusters().size()==1){
	      _phClusterEReReco->push_back(pandorapartReReco->getClusters()[0]->getEnergy());
	      //sudetector 0:ECAL
	      _phECALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]);
	      //sudetector 1:HCAL
	      _phHCALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]);
	      _phEoverTotReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapartReReco->getClusters()[0]->getEnergy());
	      //fraction of the first subcluster -> 1 cluster, obviously 100 %
	      _phf1ClusterEReReco->push_back(1.);
	      _phf1ECALReReco->push_back(1.);
	      _phf1HCALReReco->push_back(1.);
	      //now for unconverted photon check for shower track cluster difference (else we KNOW about the track anyway
	      //only consider the track here, to check for track veto etc
	      float minDistance = std::numeric_limits<float>::max();
	      float minEnergyDifference(std::numeric_limits<float>::max());
	      float sigmaPOverPClosestTrack=-1;
	      float PClosestTrackReReco=-1;
	      float PtClosestTrackReReco=-1;
	      int ind_cTrk=-1;
	      for(int t=0;t<trackCol->getNumberOfElements();t++){
		Track* track=dynamic_cast<Track*>(trackCol->getElementAt(t)) ;		
		PandoraApi::Track::Parameters trackParametersAllPhReReco;
		bool trackstatefails=false;
		try{
		  this->GetTrackStatesOld(track, trackParametersAllPhReReco);
		}catch (pandora::StatusCodeException &statusCodeException)
		  {
		    trackstatefails=true;
		    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
		  }
		if(!trackstatefails){
		  int m_maxSearchLayer=9;
		  float m_parallelDistanceCut=100;
		  float m_minTrackClusterCosAngle=0;
		  float trackClusterDistanceTest=-1;
		  this->GetTrackClusterDistance(trackParametersAllPhReReco, pandorapartReReco->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest); 	
		  if(trackClusterDistanceTest!=(-1)){
		    float trackparticleMass=trackParametersAllPhReReco.m_mass.Get();
		    float trackenergy=std::sqrt(trackparticleMass*trackparticleMass + trackParametersAllPhReReco.m_momentumAtDca.Get().GetMagnitudeSquared());
		    const float energyDifference=std::fabs(En_HAD - trackenergy);
		    if ((trackClusterDistanceTest < minDistance) || ((trackClusterDistanceTest == minDistance) && (energyDifference < minEnergyDifference))){
		      minDistance = trackClusterDistanceTest;
		      minEnergyDifference = energyDifference;
		      sigmaPOverPClosestTrack=std::sqrt(track->getCovMatrix()[5]) / std::fabs(track->getOmega());
		      PClosestTrackReReco=sqrt(trackParametersAllPhReReco.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAllPhReReco.m_momentumAtDca.Get().GetZ()*trackParametersAllPhReReco.m_momentumAtDca.Get().GetZ());
		      PtClosestTrackReReco=trackParametersAllPhReReco.m_momentumAtDca.Get().GetMagnitude();
		      ind_cTrk=t;
		    }
		  }
		}
	      }
	      if(PClosestTrackReReco<0){
		//i.e. no track found
		minDistance=-1;
	      }
	      _phSigmaPOverPClosestTrackReReco->push_back(sigmaPOverPClosestTrack);
	      _phMinTrackClustDiffReReco->push_back(minDistance);
	      _phClosestTrackPtReReco->push_back(PtClosestTrackReReco);
	      _phClosestTrackPReReco->push_back(PClosestTrackReReco);
	      if(ind_cTrk>=0){//closest track was found
		Track* cltrack=dynamic_cast<Track*>(trackCol->getElementAt(ind_cTrk)) ;			
		PandoraApi::Track::Parameters trackParametersCT;
		//std::cout<<"track state ph rereco rewind "<<std::endl;
		this->GetTrackStatesOld(cltrack, trackParametersCT);
		const pandora::CartesianVector &trackPosition=trackParametersCT.m_trackStateAtCalorimeter.Get().GetPosition();
		const float *pClusterPosition(pandorapartReReco->getClusters()[0]->getPosition());
		const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
		_phClosestTrackClusterOpeningAngleReReco->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
		//order VTX, FTD, SIT, TPC, SET, ETD
		_phClosestTrackD0ReReco->push_back(cltrack->getD0());
		_phClosestTrackZ0ReReco->push_back(cltrack->getZ0());
		_phClosestTrackPhiReReco->push_back(cltrack->getPhi());
		_phClosestTrackChi2_NDOFReReco->push_back(cltrack->getChi2()/(float)cltrack->getNdf());
		_phClosestTrackOmegaReReco->push_back(cltrack->getOmega());
		_phClosestTrackNHitsVTXReReco->push_back(cltrack->getSubdetectorHitNumbers()[0]);
		_phClosestTrackNHitsFTDReReco->push_back(cltrack->getSubdetectorHitNumbers()[1]);
		_phClosestTrackNHitsSITReReco->push_back(cltrack->getSubdetectorHitNumbers()[2]);
		_phClosestTrackNHitsTPCReReco->push_back(cltrack->getSubdetectorHitNumbers()[3]);
		_phClosestTrackNHitsSETReReco->push_back(cltrack->getSubdetectorHitNumbers()[4]);
		_phClosestTrackNHitsETDReReco->push_back(cltrack->getSubdetectorHitNumbers()[5]);
		_phClosestTrackNHitsVTXVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[6] - cltrack->getSubdetectorHitNumbers()[0]);
		_phClosestTrackNHitsFTDVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[7] - cltrack->getSubdetectorHitNumbers()[1]);
		_phClosestTrackNHitsSITVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[8] - cltrack->getSubdetectorHitNumbers()[2]);
		_phClosestTrackNHitsTPCVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[9] - cltrack->getSubdetectorHitNumbers()[3]);
		_phClosestTrackNHitsSETVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[10] - cltrack->getSubdetectorHitNumbers()[4]);
		_phClosestTrackNHitsETDVetoReReco->push_back(cltrack->getSubdetectorHitNumbers()[11] - cltrack->getSubdetectorHitNumbers()[5]);
	      }else{
		_phClosestTrackClusterOpeningAngleReReco->push_back(-1);
		_phClosestTrackD0ReReco->push_back(-1);
		_phClosestTrackZ0ReReco->push_back(-1);
		_phClosestTrackPhiReReco->push_back(-1);
		_phClosestTrackChi2_NDOFReReco->push_back(-1);
		_phClosestTrackOmegaReReco->push_back(-1);
		_phClosestTrackNHitsVTXReReco->push_back(-1);
		_phClosestTrackNHitsFTDReReco->push_back(-1);
		_phClosestTrackNHitsSITReReco->push_back(-1);
		_phClosestTrackNHitsTPCReReco->push_back(-1);
		_phClosestTrackNHitsSETReReco->push_back(-1);
		_phClosestTrackNHitsETDReReco->push_back(-1);
		_phClosestTrackNHitsVTXVetoReReco->push_back(-1);
		_phClosestTrackNHitsFTDVetoReReco->push_back(-1);
		_phClosestTrackNHitsSITVetoReReco->push_back(-1);
		_phClosestTrackNHitsTPCVetoReReco->push_back(-1);
		_phClosestTrackNHitsSETVetoReReco->push_back(-1);
		_phClosestTrackNHitsETDVetoReReco->push_back(-1);
	      }
	      float minDistanceSelTrack = std::numeric_limits<float>::max();
	      float minEnergyDifferenceSelTrack(std::numeric_limits<float>::max());
	      float sigmaPOverPClosestSelTrack=-1;
	      float PClosestSelTrack=-1;
	      float PtClosestSelTrack=-1;
	      int ind_selcTrk=-1;
	      //std::cout<<"before rereco seltrack stuff"<<std::endl;
	      for(int t=0;t<seltrackCol->getNumberOfElements();t++){
		Track* seltrack=dynamic_cast<Track*>(seltrackCol->getElementAt(t)) ;		
		PandoraApi::Track::Parameters seltrackParametersAll;	
		bool trackstatefails=false;
		try{
		  this->GetTrackStatesOld(seltrack, seltrackParametersAll);
		  
		}catch (pandora::StatusCodeException &statusCodeException)
		  {
		    trackstatefails=true;
		    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
		  }
		if(!trackstatefails){
		  int m_maxSearchLayerSelTrack=9;
		  float m_parallelDistanceCutSelTrack=100;
		  float m_minSelTrackClusterCosAngle=0;
		  float seltrackClusterDistanceTest=-1;
		  this->GetTrackClusterDistance(seltrackParametersAll, pandorapartReReco->getClusters()[0], m_maxSearchLayerSelTrack, m_parallelDistanceCutSelTrack,m_minSelTrackClusterCosAngle, seltrackClusterDistanceTest);
		  if(seltrackClusterDistanceTest!=(-1)){
		    float seltrackparticleMass=seltrackParametersAll.m_mass.Get();
		    float seltrackenergy=std::sqrt(seltrackparticleMass*seltrackparticleMass + seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared());
		    const float selenergyDifference=std::fabs(En_HAD - seltrackenergy);
		    if ((seltrackClusterDistanceTest < minDistance) || ((seltrackClusterDistanceTest == minDistance) && (selenergyDifference < minEnergyDifferenceSelTrack))){
		      minDistanceSelTrack = seltrackClusterDistanceTest;
		      minEnergyDifferenceSelTrack = selenergyDifference;
		      sigmaPOverPClosestSelTrack=std::sqrt(seltrack->getCovMatrix()[5]) / std::fabs(seltrack->getOmega());
		      PClosestSelTrack=sqrt(seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-seltrackParametersAll.m_momentumAtDca.Get().GetZ()*seltrackParametersAll.m_momentumAtDca.Get().GetZ());
		      PtClosestSelTrack=seltrackParametersAll.m_momentumAtDca.Get().GetMagnitude();
		      ind_selcTrk=t;
		    }
		  }
		}
	      }
	      if(ind_selcTrk>=0){//closest track was found
		Track* clseltrack=dynamic_cast<Track*>(seltrackCol->getElementAt(ind_selcTrk)) ;			
		PandoraApi::Track::Parameters trackParametersCT;	
		//std::cout<<"track state ph rereco rewind "<<std::endl;
		this->GetTrackStatesOld(clseltrack, trackParametersCT);
		const pandora::CartesianVector &trackPosition=trackParametersCT.m_trackStateAtCalorimeter.Get().GetPosition();
		const float *pClusterPosition(pandorapartReReco->getClusters()[0]->getPosition());
		const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0], pClusterPosition[1], pClusterPosition[2]);
		_phClosestSelTrackClusterOpeningAngleReReco->push_back(clusterInitialDirection.GetOpeningAngle(trackPosition));
	      }else{
		_phClosestSelTrackClusterOpeningAngleReReco->push_back(-1);
	      }
	      if(PClosestSelTrack<0){
		//i.e. no track found
		minDistanceSelTrack=-1;
	      }
	      _phSigmaPOverPClosestSelTrackReReco->push_back(sigmaPOverPClosestSelTrack);
	      _phClosestSelTrackPReReco->push_back(PClosestSelTrack);
	      _phClosestSelTrackPtReReco->push_back(PtClosestSelTrack);
	      _phMinSelTrackClustDiffReReco->push_back(minDistanceSelTrack);
	    }
	    if(pandorapartReReco->getTracks().size()==0){
	      _phT1Chi2_NDOFReReco  ->push_back(-1.);
	      _phT1PhiReReco        ->push_back(-10.);
	      _phT1NHitsFittedReReco->push_back(0);
	      _phT1NHitsVetoReReco  ->push_back(0);
	      _phT2Chi2_NDOFReReco  ->push_back(-1.);
	      _phT2PhiReReco        ->push_back(-10.);
	      _phT2NHitsFittedReReco->push_back(0);
	      _phT2NHitsVetoReReco  ->push_back(0);
	    }else if(pandorapartReReco->getTracks().size()==2){
	      //converted photon: track values play a role now
	      _phT1Chi2_NDOFReReco  ->push_back(pandorapartReReco->getTracks()[0]->getChi2()/(float)pandorapartReReco->getTracks()[0]->getNdf());
	      _phT1PhiReReco        ->push_back(pandorapartReReco->getTracks()[0]->getPhi());
	      _phT1NHitsFittedReReco->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]
					  +pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]
					  +pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]
					  +pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]
					  +pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]
					  +pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	      _phT1NHitsVetoReReco  ->push_back(pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[6]+
					  pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[7]+
					  pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[8]+
					  pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[9]+
					  pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[10]+
					  pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[11]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[0]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[1]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[2]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[3]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[4]
					  -pandorapartReReco->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	      _phT2Chi2_NDOFReReco  ->push_back(pandorapartReReco->getTracks()[1]->getChi2()/(float)pandorapartReReco->getTracks()[1]->getNdf());
	      _phT2PhiReReco        ->push_back(pandorapartReReco->getTracks()[1]->getPhi());
	      _phT2NHitsFittedReReco->push_back(pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[0]
					  +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[1]
					  +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[2]
					  +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[3]
					  +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[4]
					  +pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	      _phT2NHitsVetoReReco  ->push_back(pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[6]+
					  pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[7]+
					  pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[8]+
					  pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[9]+
					  pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[10]+
					  pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[11]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[0]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[1]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[2]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[3]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[4]
					  -pandorapartReReco->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	      if(pandorapartReReco->getClusters().size()==0){
		_phClusterEReReco->push_back(0.);
		//sudetector 0:ECAL
		_phECALReReco->push_back(0.);
		//sudetector 1:HCAL
		_phHCALReReco->push_back(0.);
		_phEoverTotReReco->push_back(-1.);
		//fraction of the first subcluster -> 0 cluster, obviously 0
		_phf1ClusterEReReco->push_back(0.);
		_phf1ECALReReco->push_back(0.);
		_phf1HCALReReco->push_back(0.);
	      }else if(pandorapartReReco->getClusters().size()==2){
		//sudetector 0:ECAL
		_phECALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReReco->getClusters()[1]->getSubdetectorEnergies()[0]);
		//sudetector 1:HCAL
		_phHCALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapartReReco->getClusters()[1]->getSubdetectorEnergies()[1]);
		_phEoverTotReReco->push_back((pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReReco->getClusters()[1]->getSubdetectorEnergies()[0])/(pandorapartReReco->getClusters()[0]->getEnergy()+pandorapartReReco->getClusters()[1]->getEnergy()));
		//fraction of the first subcluster -> 2 clusters, obviously NOT 100 % 
		_phf1ClusterEReReco->push_back(pandorapartReReco->getClusters()[0]->getEnergy()/(pandorapartReReco->getClusters()[0]->getEnergy()+pandorapartReReco->getClusters()[1]->getEnergy()));
		_phf1ECALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]/(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReReco->getClusters()[1]->getSubdetectorEnergies()[0]));
		_phf1HCALReReco->push_back(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]/(pandorapartReReco->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapartReReco->getClusters()[1]->getSubdetectorEnergies()[1]));
	      }else if (pandorapartReReco->getClusters().size()!=1){
		std::cout<<"converted photon, but unexpected size of clusters (neither 0,1,2) "<<pandorapartReReco->getClusters().size()<<std::endl;
	      }
	    }else{
	      std::cout<<"photon unusual track number: neither 0 (unconverted) nor 2 (conversion) "<<pandorapartReReco->getTracks().size()<<std::endl;
	    }
	    //fill energy contribution


	  }
	  //for isolation of interesting particles (i.e. photons ==22, electrons | |==11 and muons | |==13
	  if((abs(pandorapartReReco->getType())==13 && pandorapartReReco->getEnergy()>=_lepEMin) || (abs(pandorapartReReco->getType())==11 && pandorapartReReco->getEnergy()>=_lepEMin)|| (pandorapartReReco->getType()==22 && pandorapartReReco->getEnergy()>=_phEMin)){
	    //std::cout<<"i get here now el "<<pandorapartReReco->getEnergy()<<" iso sum "<<pandorapartReReco->getType()<<std::endl;
	    float sum_rereco_DR03_chargedHadrons=0;
	    float sum_rereco_DR03_neutralHadrons=0;
	    float sum_rereco_DR03_photons=0;
	    float sum_rereco_DR03_muons=0;
	    float sum_rereco_DR03_electrons=0;
	    float sum_rereco_DR04_chargedHadrons=0;
	    float sum_rereco_DR04_neutralHadrons=0;
	    float sum_rereco_DR04_photons=0;
	    float sum_rereco_DR04_muons=0;
	    float sum_rereco_DR04_electrons=0;
	    float sum_rereco_cosTheta995_chargedHadrons=0;
	    float sum_rereco_cosTheta995_neutralHadrons=0;
	    float sum_rereco_cosTheta995_photons=0;
	    float sum_rereco_cosTheta995_muons=0;
	    float sum_rereco_cosTheta995_electrons=0;
	    for(int j=0;j<recoPartColReReco->getNumberOfElements();j++){
	      if(i!=j){//exclude particle from its own sum
		ReconstructedParticle* pandorapartReReco2 = dynamic_cast<ReconstructedParticle*>(recoPartColReReco->getElementAt(j));
		TLorentzVector temp2;
		temp2.SetPxPyPzE(pandorapartReReco2->getMomentum()[0],pandorapartReReco2->getMomentum()[1],pandorapartReReco2->getMomentum()[2],pandorapartReReco2->getEnergy());
		double deltaR=temp2.DeltaR(temp);
		double cosTheta=(temp.Px()*temp2.Px()+temp.Py()*temp2.Py()+temp.Pz()*temp2.Pz())/(temp.P()*temp2.P());
		if(cosTheta>0.995){
		  if(abs(pandorapartReReco2->getType())==211){
		    sum_rereco_cosTheta995_chargedHadrons+=pandorapartReReco2->getEnergy();
		  }else if(pandorapartReReco2->getType()==22){
		    sum_rereco_cosTheta995_photons+=pandorapartReReco2->getEnergy();
		  }else if(abs(pandorapartReReco2->getType())==13){
		    sum_rereco_cosTheta995_muons+=pandorapartReReco2->getEnergy();
		  }else if(abs(pandorapartReReco2->getType())==11){
		    sum_rereco_cosTheta995_electrons+=pandorapartReReco2->getEnergy();
		  }else if(abs(pandorapartReReco2->getType())>300){
		    sum_rereco_cosTheta995_neutralHadrons+=pandorapartReReco2->getEnergy();
		  }
		}
		if(deltaR<0.4){
		  if(abs(pandorapartReReco2->getType())==211){
		    sum_rereco_DR04_chargedHadrons+=pandorapartReReco2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_chargedHadrons+=pandorapartReReco2->getEnergy();
		    }	      
		  }else if(pandorapartReReco2->getType()==22){
		    sum_rereco_DR04_photons+=pandorapartReReco2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_photons+=pandorapartReReco2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReReco2->getType())==13){
		    sum_rereco_DR04_muons+=pandorapartReReco2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_muons+=pandorapartReReco2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReReco2->getType())==11){
		    sum_rereco_DR04_electrons+=pandorapartReReco2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_electrons+=pandorapartReReco2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReReco2->getType())>300){
		    //so far neutral I found are 310 for K_short, 2112 for neutrons and 3122 for Lambdas --> should include all neutrals
		    //obviously larger than leptons, charged pions and photons
		    sum_rereco_DR04_neutralHadrons+=pandorapartReReco2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_neutralHadrons+=pandorapartReReco2->getEnergy();
		    }
		  }
		}
	      }
	    }//loop over second recopart collection for isolation calculation
	    if(abs(pandorapartReReco->getType())==11){
	      _elCHIsoDR03ReReco->push_back(sum_rereco_DR03_chargedHadrons);
	      _elCHIsoDR04ReReco->push_back(sum_rereco_DR04_chargedHadrons);
	      _elNHIsoDR03ReReco->push_back(sum_rereco_DR03_neutralHadrons);
	      _elNHIsoDR04ReReco->push_back(sum_rereco_DR04_neutralHadrons);
	      _elPhIsoDR03ReReco->push_back(sum_rereco_DR03_photons);
	      _elPhIsoDR04ReReco->push_back(sum_rereco_DR04_photons);
	      _elMuIsoDR03ReReco->push_back(sum_rereco_DR03_muons);
	      _elMuIsoDR04ReReco->push_back(sum_rereco_DR04_muons);
	      _elElIsoDR03ReReco->push_back(sum_rereco_DR03_electrons);
	      _elElIsoDR04ReReco->push_back(sum_rereco_DR04_electrons);
	      _elCHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_chargedHadrons);
	      _elNHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_neutralHadrons);
	      _elPhIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_photons);
	      _elMuIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_muons);
	      _elElIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_electrons);
	    }else if(abs(pandorapartReReco->getType())==13){
	      _muCHIsoDR03ReReco->push_back(sum_rereco_DR03_chargedHadrons);
	      _muCHIsoDR04ReReco->push_back(sum_rereco_DR04_chargedHadrons);
	      _muNHIsoDR03ReReco->push_back(sum_rereco_DR03_neutralHadrons);
	      _muNHIsoDR04ReReco->push_back(sum_rereco_DR04_neutralHadrons);
	      _muPhIsoDR03ReReco->push_back(sum_rereco_DR03_photons);
	      _muPhIsoDR04ReReco->push_back(sum_rereco_DR04_photons);
	      _muMuIsoDR03ReReco->push_back(sum_rereco_DR03_muons);
	      _muMuIsoDR04ReReco->push_back(sum_rereco_DR04_muons);
	      _muElIsoDR03ReReco->push_back(sum_rereco_DR03_electrons);
	      _muElIsoDR04ReReco->push_back(sum_rereco_DR04_electrons);
	      _muCHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_chargedHadrons);
	      _muNHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_neutralHadrons);
	      _muPhIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_photons);
	      _muMuIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_muons);
	      _muElIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_electrons);
	    }else if(pandorapartReReco->getType()==22){
	      _phCHIsoDR03ReReco->push_back(sum_rereco_DR03_chargedHadrons);
	      _phCHIsoDR04ReReco->push_back(sum_rereco_DR04_chargedHadrons);
	      _phNHIsoDR03ReReco->push_back(sum_rereco_DR03_neutralHadrons);
	      _phNHIsoDR04ReReco->push_back(sum_rereco_DR04_neutralHadrons);
	      _phPhIsoDR03ReReco->push_back(sum_rereco_DR03_photons);
	      _phPhIsoDR04ReReco->push_back(sum_rereco_DR04_photons);
	      _phMuIsoDR03ReReco->push_back(sum_rereco_DR03_muons);
	      _phMuIsoDR04ReReco->push_back(sum_rereco_DR04_muons);
	      _phElIsoDR03ReReco->push_back(sum_rereco_DR03_electrons);
	      _phElIsoDR04ReReco->push_back(sum_rereco_DR04_electrons);
	      _phCHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_chargedHadrons);
	      _phNHIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_neutralHadrons);
	      _phPhIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_photons);
	      _phMuIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_muons);
	      _phElIsoCosTheta995ReReco->push_back(sum_rereco_cosTheta995_electrons);
	    }
	  }//recopart is either muon, electron or photon
	}
      }

      LCCollection* jetColReReco = NULL;
      try{
	evt->getCollection( _jetColNameReReco ) ;
      }catch( lcio::DataNotAvailableException e ){
	std::cout<<"in run/event/weight no rereco jets "<<evt->getRunNumber ()<<"/"<<evt->getEventNumber()<<"/"<<evt->getWeight()<<"/"<<_lepEMin<<"/"<<"/"<<_phEMin<<"/"<<_jetEMin<<std::endl;
      }
      if(jetColReReco!=NULL){
	//std::cout<<"rereco jet"<<std::endl;
	for(int i=0;i<jetColReReco->getNumberOfElements();i++){
	  ReconstructedParticle* recojetReReco = dynamic_cast<ReconstructedParticle*>(jetColReReco->getElementAt(i));
	  if(recojetReReco->getEnergy()>_jetEMin){
	    _jetEReReco->push_back(recojetReReco->getEnergy());
	    _jetPxReReco->push_back(recojetReReco->getMomentum()[0]);
	    _jetPyReReco->push_back(recojetReReco->getMomentum()[1]);
	    _jetPzReReco->push_back(recojetReReco->getMomentum()[2]);
	    TLorentzVector temp;
	    temp.SetPxPyPzE(recojetReReco->getMomentum()[0],recojetReReco->getMomentum()[1],recojetReReco->getMomentum()[2],recojetReReco->getEnergy());
	    _jetPhiReReco->push_back(temp.Phi());
	    _jetThetaReReco->push_back(temp.Theta());
	    _jetMassReReco->push_back(recojetReReco->getMass());
	    _jetNFReReco->push_back(0);
	    _jetLFReReco->push_back(0);
	    _jetKFReReco->push_back(0);
	    _jetCHFReReco->push_back(0);
	    _jetPhFReReco->push_back(0);
	    _jetElFReReco->push_back(0);
	    _jetMuFReReco->push_back(0);
	    _jetNMultReReco->push_back(0);
	    _jetLMultReReco->push_back(0);
	    _jetKMultReReco->push_back(0);
	    _jetCHMultReReco->push_back(0);
	    _jetPhMultReReco->push_back(0);
	    _jetElMultReReco->push_back(0);
	    _jetMuMultReReco->push_back(0);
	    unsigned int jet_index=_jetPhiReReco->size()-1;
	    for(unsigned int j=0;j<recojetReReco->getParticles().size();j++){
	      if(abs(recojetReReco->getParticles()[j]->getType())==211){ 
		(*_jetCHFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetCHMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==22){
		(*_jetPhFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetPhMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==11){
		(*_jetElFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetElMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==13){
		(*_jetMuFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetMuMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==2112){
		(*_jetNFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetNMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==3122){
		(*_jetLFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetLMultReReco)[jet_index]+=1;
	      }else if(abs(recojetReReco->getParticles()[j]->getType())==310 || abs(recojetReReco->getParticles()[j]->getType())==130 || abs(recojetReReco->getParticles()[j]->getType())==311){
		std::cout<<"somewhere in k0 "<<recojetReReco->getParticles()[j]->getType()<<" mass "<<recojetReReco->getParticles()[j]->getMass()<<std::endl; 
		(*_jetKFReReco)[jet_index]+=recojetReReco->getParticles()[j]->getEnergy();
		(*_jetKMultReReco)[jet_index]+=1;
	      }else{
		std::cout<<"don't take into account particle type "<<recojetReReco->getParticles()[j]->getType()<<std::endl;
	      }
	      //&& recojetReReco->getParticles()[j]->getType()!=3122 && recojetReReco->getParticles()[j]->getType()!=310){
	      //std::cout<<"part "<<j<<" ID,mass,E/px/py/pz"<<recojetReReco->getParticles()[j]->getType()<<"/"<<std::endl;
	      //}
	    }
	    float tot_energy=(*_jetCHFReReco)[jet_index]+(*_jetPhFReReco)[jet_index]+(*_jetMuFReReco)[jet_index]+(*_jetElFReReco)[jet_index]+(*_jetNFReReco)[jet_index]+(*_jetLFReReco)[jet_index]+(*_jetKFReReco)[jet_index];
	    (*_jetCHFReReco)[jet_index]=(*_jetCHFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetPhFReReco)[jet_index]=(*_jetPhFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetMuFReReco)[jet_index]=(*_jetMuFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetElFReReco)[jet_index]=(*_jetElFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetNFReReco)[jet_index]=(*_jetNFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetLFReReco)[jet_index]=(*_jetLFReReco)[jet_index]/recojetReReco->getEnergy();
	    (*_jetKFReReco)[jet_index]=(*_jetKFReReco)[jet_index]/recojetReReco->getEnergy();
	    float tot_fraction=(*_jetCHFReReco)[jet_index]+(*_jetPhFReReco)[jet_index]+(*_jetMuFReReco)[jet_index]+(*_jetElFReReco)[jet_index]+(*_jetNFReReco)[jet_index]+(*_jetLFReReco)[jet_index]+(*_jetKFReReco)[jet_index];
	    if(fabs(1.-(*_jetCHFReReco)[jet_index]-(*_jetPhFReReco)[jet_index]-(*_jetMuFReReco)[jet_index]-(*_jetElFReReco)[jet_index]-(*_jetNFReReco)[jet_index]-(*_jetLFReReco)[jet_index]-(*_jetKFReReco)[jet_index])>1.e-4){
	      std::cout<<"big difference in fractions jetReREco "<<fabs(1.-(*_jetCHFReReco)[jet_index]-(*_jetPhFReReco)[jet_index]-(*_jetMuFReReco)[jet_index]-(*_jetElFReReco)[jet_index]-(*_jetNFReReco)[jet_index]-(*_jetLFReReco)[jet_index]-(*_jetKFReReco)[jet_index])<<"/"<<tot_energy<<"/"<<recojetReReco->getEnergy()<<"/"<<tot_energy/recojetReReco->getEnergy()<<"/"<<tot_fraction<<std::endl;
	    }
	  }
	}
      }
      _METReReco=sqrt(pow(_MExReReco,2)+pow(_MEyReReco,2));
    }
 
     //std::cout<<"i get here now el "<<pandorapartReRecoNoTrkSigmaPOverP->getEnergy()<<" iso sum "<<pandorapartReRecoNoTrkSigmaPOverP->getType()<<std::endl;


    ////////no rereco stuff without tracks

    if(_hasrerecoNoTrackSigmaPOverP_information){
      //std::cout<<"rereco sigma "<<std::endl;
      LCCollection* recoPartColReRecoNoTrkSigmaPOverP = evt->getCollection( _recoPartNameReReco ) ;
      if(recoPartColReRecoNoTrkSigmaPOverP!=NULL){
	//PandoraCandidate loop
	for(int i=0;i<recoPartColReRecoNoTrkSigmaPOverP->getNumberOfElements();i++){
	  ReconstructedParticle* pandorapartReRecoNoTrkSigmaPOverP = dynamic_cast<ReconstructedParticle*>(recoPartColReRecoNoTrkSigmaPOverP->getElementAt(i));
	  TLorentzVector temp;
	  temp.SetPxPyPzE(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[0],pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[1],pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[2],pandorapartReRecoNoTrkSigmaPOverP->getEnergy());

	  if(abs(pandorapartReRecoNoTrkSigmaPOverP->getType())==11 && pandorapartReRecoNoTrkSigmaPOverP->getEnergy()>=_lepEMin){
	    _elIDReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getType());
	    _elEReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getEnergy());
	    _elPxReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[0]);
	    _elPyReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[1]);
	    _elPzReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[2]);
	    _elPhiReRecoNoTrkSigmaPOverP->push_back(temp.Phi());
	    _elThetaReRecoNoTrkSigmaPOverP->push_back(temp.Theta());
	    _elMassReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMass());
	    //start here eloton specific cluster/composition/shower shape etc etc
	    float En_EM=0;
	    float En_HAD=0;
	    float En_MUON=0;	  
	    float En_ECAL_Barrel=0;
	    float En_ECAL_Endcap=0;
	    float En_ECAL_else=0;
	    float En_HCAL_Barrel=0;
	    float En_HCAL_Endcap=0;
	    float En_HCAL_else=0;
	    float En_ELSE_Barrel=0;
	    float En_ELSE_Endcap=0;
	    float En_ELSE_else=0;
	    //*_elLongShowerProfileReRecoNoTrkSigmaPOverP,*_elTransShowerProfileReRecoNoTrkSigmaPOverP,*_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*_elPeakEnergyReRecoNoTrkSigmaPOverP,*_elNRadiationLengthReRecoNoTrkSigmaPOverP;
	    int elFirstLayerHCAL=-1;
	    int elFirstLayerECAL=-1;
	    int elLastLayerHCAL=-1;
	    int elLastLayerECAL=-1;
	    int hitsECAL=0;
	    int hitsHCAL=0;
	    //40 layers decided at the moment, 20+10 in samples
	    std::vector<int>hits_per_layer_ecal(50,0);
	    std::vector<float>energy_per_layer_ecal(50,0);
	    //75 layers at the moment
	    std::vector<int>hits_per_layer_hcal(100,0);
	    std::vector<float>energy_per_layer_hcal(100,0);
	    for(unsigned int j=0;j<pandorapartReRecoNoTrkSigmaPOverP->getClusters().size();j++){
	      for(unsigned int l=0; l<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits().size();l++){
		TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
		tempP.SetXYZT(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
		int type=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getType();
		static const int fCaloType =     1 ;
		static const int fCaloID   =    10 ;
		static const int fLayout   =  1000 ;
		static const int fLayer    = 10000 ;
		int caloLayer=type/fLayer;
		int caloLayout=(type-caloLayer*fLayer)/fLayout;
		int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
		int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
		//std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
		if(caloType==0){//0 em, 1 had, 2 mu
		  En_EM+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==1){     
		  En_HAD+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==2){
		  En_MUON+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1){//ecal
		  if(elFirstLayerECAL<0){
		    elFirstLayerECAL=caloLayer;
		  }else if (caloLayer<elFirstLayerECAL){
		    elFirstLayerECAL=caloLayer;
		  }
		  if(elLastLayerECAL<0){
		    elLastLayerECAL=caloLayer;
		  }else if (caloLayer>elLastLayerECAL){
		    elLastLayerECAL=caloLayer;
		  }
		  hitsECAL+=1;
		  hits_per_layer_ecal[caloLayer]+=1;
		  energy_per_layer_ecal[caloLayer]+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if (caloID==2){
		  if(elFirstLayerHCAL<0){
		    elFirstLayerHCAL=caloLayer;
		  }else if (caloLayer<elFirstLayerHCAL){
		    elFirstLayerHCAL=caloLayer;
		  }
		  if(elLastLayerHCAL<0){
		    elLastLayerHCAL=caloLayer;
		  }else if (caloLayer>elLastLayerHCAL){
		    elLastLayerHCAL=caloLayer;
		  }
		  hitsHCAL+=1;
		  hits_per_layer_hcal[caloLayer]+=1;
		  energy_per_layer_hcal[caloLayer]+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_ECAL_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==1 && caloLayout==2){
		  En_ECAL_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==1){
		  En_ECAL_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==1){
		  En_HCAL_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==2){
		  En_HCAL_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==2){
		  En_HCAL_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==1){
		  En_ELSE_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==2){
		  En_ELSE_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID!=1 && caloID!=2){
		  En_ELSE_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
	      }
	    }//loop over all clusters and hits for energy in different detectors
	    //*_elLongShowerProfileReRecoNoTrkSigmaPOverP,*_elTransShowerProfileReRecoNoTrkSigmaPOverP,*_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*_elPeakEnergyReRecoNoTrkSigmaPOverP,*_elNRadiationLengthReRecoNoTrkSigmaPOverP;
	    _elHitsECALReRecoNoTrkSigmaPOverP->push_back(hitsECAL);
	    _elFirstLayerECALReRecoNoTrkSigmaPOverP->push_back(elFirstLayerECAL);
	    _elLastLayerECALReRecoNoTrkSigmaPOverP->push_back(elLastLayerECAL);
	    int maxHitsPerLayer=0;
	    float maxEnergyPerLayer=0;
	    if(hitsECAL>0){
	      _elLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _elLayerMaxEECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		  (*_elLayerMaxHitsECALReRecoNoTrkSigmaPOverP)[_elLayerMaxHitsECALReRecoNoTrkSigmaPOverP->size()-1]=j;
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		  (*_elLayerMaxEECALReRecoNoTrkSigmaPOverP)[_elLayerMaxEECALReRecoNoTrkSigmaPOverP->size()-1]=j;
		}
	      }
	      _elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP->push_back(maxHitsPerLayer);
	      _elMaxEInECALLayerReRecoNoTrkSigmaPOverP->push_back(maxEnergyPerLayer);
	    }else{//no ECAL hits
	      _elLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _elLayerMaxEECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP->push_back(0);
	      _elMaxEInECALLayerReRecoNoTrkSigmaPOverP->push_back(0);
	    }
	    _elHitsHCALReRecoNoTrkSigmaPOverP->push_back(hitsHCAL);
	    _elFirstLayerHCALReRecoNoTrkSigmaPOverP->push_back(elFirstLayerHCAL);
	    _elLastLayerHCALReRecoNoTrkSigmaPOverP->push_back(elLastLayerHCAL);
	    //initialize to 0 again for HCAL
	    maxHitsPerLayer=0;
	    maxEnergyPerLayer=0;
	    if(hitsHCAL>0){
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		}
	      }
	      _elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP->push_back(maxHitsPerLayer);
	      _elMaxEInHCALLayerReRecoNoTrkSigmaPOverP->push_back(maxEnergyPerLayer);
	    }else{//no HCAL hits
	      _elLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP->push_back(0);
	      _elMaxEInHCALLayerReRecoNoTrkSigmaPOverP->push_back(0);
	    }
	    _elEnergyEBReRecoNoTrkSigmaPOverP->push_back(En_ECAL_Barrel);
	    _elEnergyEEReRecoNoTrkSigmaPOverP->push_back(En_ECAL_Endcap);
	    _elEnergyEElseReRecoNoTrkSigmaPOverP->push_back(En_ECAL_else);
	    _elEnergyHEReRecoNoTrkSigmaPOverP->push_back(En_HCAL_Barrel);
	    _elEnergyHBReRecoNoTrkSigmaPOverP->push_back(En_HCAL_Endcap);
	    _elEnergyHElseReRecoNoTrkSigmaPOverP->push_back( En_HCAL_else);
	    _elEnergyElseBReRecoNoTrkSigmaPOverP->push_back(En_ELSE_Barrel);
	    _elEnergyElseEReRecoNoTrkSigmaPOverP->push_back(En_ELSE_Endcap);
	    _elEnergyElseElseReRecoNoTrkSigmaPOverP->push_back(En_ELSE_else);
	    _elEnergyEMReRecoNoTrkSigmaPOverP->push_back(En_EM);
	    _elEnergyHADReRecoNoTrkSigmaPOverP->push_back(En_HAD);
	    _elEnergyMuonReRecoNoTrkSigmaPOverP->push_back(En_MUON);
	    Track* eltrack=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0];
	    PandoraApi::Track::Parameters trackParametersAll;	
	    //std::cout<<"track state e no sigm "<<std::endl;
	    this->GetTrackStatesOld(eltrack, trackParametersAll);
	    _elSigmaPOverPTrackReRecoNoTrkSigmaPOverP->push_back(std::sqrt(eltrack->getCovMatrix()[5])/std::fabs(eltrack->getOmega()));
	    _elTrackPtReRecoNoTrkSigmaPOverP->push_back(sqrt(trackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAll.m_momentumAtDca.Get().GetZ()*trackParametersAll.m_momentumAtDca.Get().GetZ()));
	    _elTrackPReRecoNoTrkSigmaPOverP->push_back(trackParametersAll.m_momentumAtDca.Get().GetMagnitude());
	    if(pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()==1){
	      _elClusterEReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy());
	      //sudetector 0:ECAL
	      _elECALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]);
	      //sudetector 1:HCAL
	      _elHCALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[1]);
	      _elEoverTotReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy());
	      //unlike for photons for electrons we KNOW about the track anyway
	      //only consider this track here, to check for track veto etc		
	      int m_maxSearchLayer=9;
	      float m_parallelDistanceCut=100;
	      float m_minTrackClusterCosAngle=0;
	      float trackClusterDistanceTest=-1;
	      this->GetTrackClusterDistance(trackParametersAll, pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest);
	      if(trackClusterDistanceTest==(-1)){
		std::cout<<"no track cluster distance -> el rereco should HAVE this track as closest"<<std::endl;
	      }
	      _elMinTrackClustDiffReRecoNoTrkSigmaPOverP->push_back(trackClusterDistanceTest);
	    }else if(pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()!=1){
	      std::cout<<"el rereco has not exactly one cluster"<<std::endl;
	      //if(fabs(temp.Theta())>2.4 && temp.Energy()>10){
	      for(unsigned int j=0;j<pandorapartReRecoNoTrkSigmaPOverP->getClusters().size();j++){
		for(unsigned int k=0;k<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getSubdetectorEnergies().size();k++){
		  std::cout<<" el rereco in cluster["<<j<<" ] of "<<temp.Theta()<<" "<<pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()<<" E tot "<<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getEnergy()<<" E sub "<<k <<" "<<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getSubdetectorEnergies()[k]<<" subcluster size? "<<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getClusters().size()<<std::endl;
		}
	      }
	    }
	    if(pandorapartReRecoNoTrkSigmaPOverP->getTracks().size()>0){
	      _el0D0ReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getD0());
	      _el0Z0ReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getZ0());
	      _el0PhiReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getPhi());
	      _el0Chi2_NDOFReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getChi2()/(float)pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getNdf());
	      _elOmegaReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getOmega());
	      _el1D0ReRecoNoTrkSigmaPOverP->push_back(-1);
	      _el1Z0ReRecoNoTrkSigmaPOverP->push_back(-1);
	      _el1PhiReRecoNoTrkSigmaPOverP->push_back(-10);
	      _el1Chi2_NDOFReRecoNoTrkSigmaPOverP->push_back(-1);
	      _el1NHitsReRecoNoTrkSigmaPOverP->push_back(-1);
	      _el1NHitsVetoReRecoNoTrkSigmaPOverP->push_back(-1);
	      if(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers().size()==12){
		//at the moment 0-5 are NHits used in the fit,6-11 are all hits including those not used
		//order VTX, FTD, SIT, TPC, SET, ETD
		_elNHitsVTXReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_elNHitsFTDReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_elNHitsSITReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_elNHitsTPCReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_elNHitsSETReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_elNHitsETDReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		_elNHitsVTXVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[6]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[0]);
		_elNHitsFTDVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[7]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[1]);
		_elNHitsSITVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[8]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[2]);
		_elNHitsTPCVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[9]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[3]);
		_elNHitsSETVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[10]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[4]);
		_elNHitsETDVetoReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[11]+
						 -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[5]);
		//ECAL is 0, HCAL is 1 in calo subdetector stuff
		if(pandorapartReRecoNoTrkSigmaPOverP->getTracks().size()>1){
		  (*_el1D0ReRecoNoTrkSigmaPOverP)[_el1D0ReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getD0();
		  (*_el1Z0ReRecoNoTrkSigmaPOverP)[_el1Z0ReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getZ0();
		  (*_el1PhiReRecoNoTrkSigmaPOverP)[_el1PhiReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getPhi();
		  (*_el1Chi2_NDOFReRecoNoTrkSigmaPOverP)[_el1Chi2_NDOFReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getChi2()/(float)pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getNdf();
		  (*_el1NHitsReRecoNoTrkSigmaPOverP)[_el1NHitsReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[0]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[1]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[2]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[3]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[4]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[5];
		  (*_el1NHitsVetoReRecoNoTrkSigmaPOverP)[_el1NHitsVetoReRecoNoTrkSigmaPOverP->size()-1]=pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[6]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[7]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[8]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[9]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[10]
		    +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[11]
		    -(*_el1NHitsReRecoNoTrkSigmaPOverP)[_el1NHitsReRecoNoTrkSigmaPOverP->size()-1];
		}
	      }else if(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers().size()<12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation"<<pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }else if(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers().size()>12){
		std::cout<<"size of sudetector tracker hits changed, check for new model and implementation, maybe 9 very soon ? "<<pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers().size()<<std::endl;
	      }	
	    }
	  }else if(pandorapartReRecoNoTrkSigmaPOverP->getType()==22 && pandorapartReRecoNoTrkSigmaPOverP->getEnergy()>=_phEMin){
	    //std::cout<<"rereco sigma ph"<<std::endl;
	    _phEReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getEnergy());
	    _phPxReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[0]);
	    _phPyReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[1]);
	    _phPzReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMomentum()[2]);
	    _phPhiReRecoNoTrkSigmaPOverP->push_back(temp.Phi());
	    _phThetaReRecoNoTrkSigmaPOverP->push_back(temp.Theta());
	    _phMassReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getMass());
	    //start here photon specific cluster/composition/shower shape etc etc
	    float En_EM=0;
	    float En_HAD=0;
	    float En_MUON=0;	  
	    float En_ECAL_Barrel=0;
	    float En_ECAL_Endcap=0;
	    float En_ECAL_else=0;
	    float En_HCAL_Barrel=0;
	    float En_HCAL_Endcap=0;
	    float En_HCAL_else=0;
	    float En_ELSE_Barrel=0;
	    float En_ELSE_Endcap=0;
	    float En_ELSE_else=0;
	    //*_phLongShowerProfileReRecoNoTrkSigmaPOverP,*_phTransShowerProfileReRecoNoTrkSigmaPOverP,*_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*_phPeakEnergyReRecoNoTrkSigmaPOverP,*_phNRadiationLengthReRecoNoTrkSigmaPOverP;
	    int phFirstLayerHCAL=-1;
	    int phFirstLayerECAL=-1;
	    int phLastLayerHCAL=-1;
	    int phLastLayerECAL=-1;
	    int hitsECAL=0;
	    int hitsHCAL=0;
	    //40 layers decided at the moment, 20+10 in samples
	    std::vector<int>hits_per_layer_ecal(50,0);
	    std::vector<float>energy_per_layer_ecal(50,0);
	    //75 layers at the moment
	    std::vector<int>hits_per_layer_hcal(100,0);
	    std::vector<float>energy_per_layer_hcal(100,0);
	    for(unsigned int j=0;j<pandorapartReRecoNoTrkSigmaPOverP->getClusters().size();j++){
	      for(unsigned int l=0; l<pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits().size();l++){
		TLorentzVector tempP; //works get correct values for phi/eta/energy (calculate R yourself)
		tempP.SetXYZT(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[0],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[1],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getPosition()[2],pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getTime());		
		int type=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getType();
		static const int fCaloType =     1 ;
		static const int fCaloID   =    10 ;
		static const int fLayout   =  1000 ;
		static const int fLayer    = 10000 ;
		int caloLayer=type/fLayer;
		int caloLayout=(type-caloLayer*fLayer)/fLayout;
		int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
		int caloType=(type-caloLayer*fLayer-caloLayout*fLayout-caloID*fCaloID)/fCaloType;
		//std::cout<<"cluster "<<j<<" hit "<<l<<" type/caloType/caloID/calolayer "<<type<<"/"<<caloType<<"/"<<caloID<<"/"<<caloLayout<<"/"<<caloLayer<<std::endl;
		if(caloType==0){//0 em, 1 had, 2 mu
		  En_EM+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==1){     
		  En_HAD+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloType==2){
		  En_MUON+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1){//ecal
		  if(phFirstLayerECAL<0){
		    phFirstLayerECAL=caloLayer;
		  }else if (caloLayer<phFirstLayerECAL){
		    phFirstLayerECAL=caloLayer;
		  }
		  if(phLastLayerECAL<0){
		    phLastLayerECAL=caloLayer;
		  }else if (caloLayer>phLastLayerECAL){
		    phLastLayerECAL=caloLayer;
		  }
		  hitsECAL+=1;
		  hits_per_layer_ecal[caloLayer]+=1;
		  energy_per_layer_ecal[caloLayer]+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if (caloID==2){
		  if(phFirstLayerHCAL<0){
		    phFirstLayerHCAL=caloLayer;
		  }else if (caloLayer<phFirstLayerHCAL){
		    phFirstLayerHCAL=caloLayer;
		  }
		  if(phLastLayerHCAL<0){
		    phLastLayerHCAL=caloLayer;
		  }else if (caloLayer>phLastLayerHCAL){
		    phLastLayerHCAL=caloLayer;
		  }
		  hitsHCAL+=1;
		  hits_per_layer_hcal[caloLayer]+=1;
		  energy_per_layer_hcal[caloLayer]+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
		if(caloID==1 && caloLayout==1){//id 1=ecal, 2 hcal, layout 1 barrel, 2 endcap
		  En_ECAL_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==1 && caloLayout==2){
		  En_ECAL_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==1){
		  En_ECAL_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==1){
		  En_HCAL_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if(caloID==2 && caloLayout==2){
		  En_HCAL_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID==2){
		  En_HCAL_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==1){
		  En_ELSE_Barrel+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else if((caloID!=1 && caloID!=2) && caloLayout==2){
		  En_ELSE_Endcap+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}else	if(caloID!=1 && caloID!=2){
		  En_ELSE_else+=pandorapartReRecoNoTrkSigmaPOverP->getClusters()[j]->getCalorimeterHits()[l]->getEnergy();
		}
	      }
	    }//loop over all clusters and hits for energy in different detectors
	    //*_phLongShowerProfileReRecoNoTrkSigmaPOverP,*_phTransShowerProfileReRecoNoTrkSigmaPOverP,*_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP,*_phPeakEnergyReRecoNoTrkSigmaPOverP,*_phNRadiationLengthReRecoNoTrkSigmaPOverP;
	    _phHitsECALReRecoNoTrkSigmaPOverP->push_back(hitsECAL);
	    _phFirstLayerECALReRecoNoTrkSigmaPOverP->push_back(phFirstLayerECAL);
	    _phLastLayerECALReRecoNoTrkSigmaPOverP->push_back(phLastLayerECAL);
	    int maxHitsPerLayer=0;
	    float maxEnergyPerLayer=0;
	    if(hitsECAL>0){
	      _phLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _phLayerMaxEECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		  (*_phLayerMaxHitsECALReRecoNoTrkSigmaPOverP)[_phLayerMaxHitsECALReRecoNoTrkSigmaPOverP->size()-1]=j;
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		  (*_phLayerMaxEECALReRecoNoTrkSigmaPOverP)[_phLayerMaxEECALReRecoNoTrkSigmaPOverP->size()-1]=j;
		}
	      }
	      _phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP->push_back(maxHitsPerLayer);
	      _phMaxEInECALLayerReRecoNoTrkSigmaPOverP->push_back(maxEnergyPerLayer);
	    }else{//no ECAL hits
	      _phLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _phLayerMaxEECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP->push_back(0);
	      _phMaxEInECALLayerReRecoNoTrkSigmaPOverP->push_back(0);
	    }
	    _phHitsHCALReRecoNoTrkSigmaPOverP->push_back(hitsHCAL);
	    _phFirstLayerHCALReRecoNoTrkSigmaPOverP->push_back(phFirstLayerHCAL);
	    _phLastLayerHCALReRecoNoTrkSigmaPOverP->push_back(phLastLayerHCAL);
	    //initialize to 0 again for HCAL
	    maxHitsPerLayer=0;
	    maxEnergyPerLayer=0;
	    if(hitsHCAL>0){
	      for(unsigned int j=0;j<hits_per_layer_ecal.size();j++){
		if(hits_per_layer_ecal[j]>maxHitsPerLayer){
		  maxHitsPerLayer=hits_per_layer_ecal[j];
		}
		if(energy_per_layer_ecal[j]>maxEnergyPerLayer){
		  maxEnergyPerLayer=energy_per_layer_ecal[j];
		}
	      }
	      _phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP->push_back(maxHitsPerLayer);
	      _phMaxEInHCALLayerReRecoNoTrkSigmaPOverP->push_back(maxEnergyPerLayer);
	    }else{//no HCAL hits
	      _phLayerMaxHitsECALReRecoNoTrkSigmaPOverP->push_back(-1);
	      _phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP->push_back(0);
	      _phMaxEInHCALLayerReRecoNoTrkSigmaPOverP->push_back(0);
	    }
	    _phEnergyEBReRecoNoTrkSigmaPOverP->push_back(En_ECAL_Barrel);
	    _phEnergyEEReRecoNoTrkSigmaPOverP->push_back(En_ECAL_Endcap);
	    _phEnergyEElseReRecoNoTrkSigmaPOverP->push_back(En_ECAL_else);
	    _phEnergyHEReRecoNoTrkSigmaPOverP->push_back(En_HCAL_Barrel);
	    _phEnergyHBReRecoNoTrkSigmaPOverP->push_back(En_HCAL_Endcap);
	    _phEnergyHElseReRecoNoTrkSigmaPOverP->push_back( En_HCAL_else);
	    _phEnergyElseBReRecoNoTrkSigmaPOverP->push_back(En_ELSE_Barrel);
	    _phEnergyElseEReRecoNoTrkSigmaPOverP->push_back(En_ELSE_Endcap);
	    _phEnergyElseElseReRecoNoTrkSigmaPOverP->push_back(En_ELSE_else);
	    _phEnergyEMReRecoNoTrkSigmaPOverP->push_back(En_EM);
	    _phEnergyHADReRecoNoTrkSigmaPOverP->push_back(En_HAD);
	    _phEnergyMuonReRecoNoTrkSigmaPOverP->push_back(En_MUON);
	    //until here for photonID/fake/cluster shape etc etc
	    if(pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()==1){
	      _phClusterEReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy());
	      //sudetector 0:ECAL
	      _phECALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]);
	      //sudetector 1:HCAL
	      _phHCALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[1]);
	      _phEoverTotReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]/pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy());
	      //fraction of the first subcluster -> 1 cluster, obviously 100 %
	      _phf1ClusterEReRecoNoTrkSigmaPOverP->push_back(1.);
	      _phf1ECALReRecoNoTrkSigmaPOverP->push_back(1.);
	      _phf1HCALReRecoNoTrkSigmaPOverP->push_back(1.);
	      //now for unconverted photon check for shower track cluster difference (else we KNOW about the track anyway
	      //only consider the track here, to check for track veto etc
	      float minDistance = std::numeric_limits<float>::max();
	      float minEnergyDifference(std::numeric_limits<float>::max());
	      float sigmaPOverPClosestTrack=-1;
	      float PClosestTrackReRecoNoTrkSigmaPOverP=-1;
	      float PtClosestTrackReRecoNoTrkSigmaPOverP=-1;
	      for(int t=0;t<trackCol->getNumberOfElements();t++){
		Track* track=dynamic_cast<Track*>(trackCol->getElementAt(t)) ;		
		PandoraApi::Track::Parameters trackParametersAllPhReRecoNoTrkSigmaPOverP;
		bool trackstatefails=false;
		try{
		  this->GetTrackStatesOld(track, trackParametersAllPhReRecoNoTrkSigmaPOverP);
		  
		}catch (pandora::StatusCodeException &statusCodeException)
		  {
		    trackstatefails=true;
		    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
		  }
		if(!trackstatefails){
		  int m_maxSearchLayer=9;
		  float m_parallelDistanceCut=100;
		  float m_minTrackClusterCosAngle=0;
		  float trackClusterDistanceTest=-1;
		  this->GetTrackClusterDistance(trackParametersAllPhReRecoNoTrkSigmaPOverP, pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0], m_maxSearchLayer, m_parallelDistanceCut,m_minTrackClusterCosAngle, trackClusterDistanceTest);
		  if(trackClusterDistanceTest!=(-1)){
		    float trackparticleMass=trackParametersAllPhReRecoNoTrkSigmaPOverP.m_mass.Get();
		    float trackenergy=std::sqrt(trackparticleMass*trackparticleMass + trackParametersAllPhReRecoNoTrkSigmaPOverP.m_momentumAtDca.Get().GetMagnitudeSquared());
		    const float energyDifference=std::fabs(En_HAD - trackenergy);
		    if ((trackClusterDistanceTest < minDistance) || ((trackClusterDistanceTest == minDistance) && (energyDifference < minEnergyDifference))){
		      minDistance = trackClusterDistanceTest;
		      minEnergyDifference = energyDifference;
		      sigmaPOverPClosestTrack=std::sqrt(track->getCovMatrix()[5]) / std::fabs(track->getOmega());
		      PClosestTrackReRecoNoTrkSigmaPOverP=sqrt(trackParametersAllPhReRecoNoTrkSigmaPOverP.m_momentumAtDca.Get().GetMagnitudeSquared()-trackParametersAllPhReRecoNoTrkSigmaPOverP.m_momentumAtDca.Get().GetZ()*trackParametersAllPhReRecoNoTrkSigmaPOverP.m_momentumAtDca.Get().GetZ());
		      PtClosestTrackReRecoNoTrkSigmaPOverP=trackParametersAllPhReRecoNoTrkSigmaPOverP.m_momentumAtDca.Get().GetMagnitude();
		    }
		  }
		}
	    }
	      if(PClosestTrackReRecoNoTrkSigmaPOverP<0){
		//i.e. no track found
		minDistance=-1;
	      }
	      _phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP->push_back(sigmaPOverPClosestTrack);
	      _phMinTrackClustDiffReRecoNoTrkSigmaPOverP->push_back(minDistance);
	      _phClosestTrackPtReRecoNoTrkSigmaPOverP->push_back(PtClosestTrackReRecoNoTrkSigmaPOverP);
	      _phClosestTrackPReRecoNoTrkSigmaPOverP->push_back(PClosestTrackReRecoNoTrkSigmaPOverP);
	      float minDistanceSelTrack = std::numeric_limits<float>::max();
	      float minEnergyDifferenceSelTrack(std::numeric_limits<float>::max());
	      float sigmaPOverPClosestSelTrack=-1;
	      float PClosestSelTrack=-1;
	      float PtClosestSelTrack=-1;
	      for(int t=0;t<seltrackCol->getNumberOfElements();t++){
		Track* seltrack=dynamic_cast<Track*>(seltrackCol->getElementAt(t)) ;		
		PandoraApi::Track::Parameters seltrackParametersAll;
		bool trackstatefails=false;
		try{
		  this->GetTrackStatesOld(seltrack, seltrackParametersAll);	
		}catch (pandora::StatusCodeException &statusCodeException)
		  {
		    trackstatefails=true;
		    streamlog_out(ERROR) << "Failed to extract a track: " << statusCodeException.ToString() << std::endl;
		  }
		if(!trackstatefails){
		  int m_maxSearchLayerSelTrack=9;
		  float m_parallelDistanceCutSelTrack=100;
		  float m_minSelTrackClusterCosAngle=0;
		  float seltrackClusterDistanceTest=-1;
		  this->GetTrackClusterDistance(seltrackParametersAll, pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0], m_maxSearchLayerSelTrack, m_parallelDistanceCutSelTrack,m_minSelTrackClusterCosAngle, seltrackClusterDistanceTest);
		  if(seltrackClusterDistanceTest!=(-1)){
		    float seltrackparticleMass=seltrackParametersAll.m_mass.Get();
		    float seltrackenergy=std::sqrt(seltrackparticleMass*seltrackparticleMass + seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared());
		    const float selenergyDifference=std::fabs(En_HAD - seltrackenergy);
		    if ((seltrackClusterDistanceTest < minDistance) || ((seltrackClusterDistanceTest == minDistance) && (selenergyDifference < minEnergyDifferenceSelTrack))){
		      minDistanceSelTrack = seltrackClusterDistanceTest;
		      minEnergyDifferenceSelTrack = selenergyDifference;
		      sigmaPOverPClosestSelTrack=std::sqrt(seltrack->getCovMatrix()[5]) / std::fabs(seltrack->getOmega());
		      PClosestSelTrack=sqrt(seltrackParametersAll.m_momentumAtDca.Get().GetMagnitudeSquared()-seltrackParametersAll.m_momentumAtDca.Get().GetZ()*seltrackParametersAll.m_momentumAtDca.Get().GetZ());
		      PtClosestSelTrack=seltrackParametersAll.m_momentumAtDca.Get().GetMagnitude();
		    }
		  }
		}
	      }
	      if(PClosestSelTrack<0){
		//i.e. no track found
		minDistanceSelTrack=-1;
	      }
	      _phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP->push_back(sigmaPOverPClosestSelTrack);
	      _phClosestSelTrackPReRecoNoTrkSigmaPOverP->push_back(PClosestSelTrack);
	      _phClosestSelTrackPtReRecoNoTrkSigmaPOverP->push_back(PtClosestSelTrack);
	      _phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP->push_back(minDistanceSelTrack);
	    }
	    if(pandorapartReRecoNoTrkSigmaPOverP->getTracks().size()==0){
	      _phT1Chi2_NDOFReRecoNoTrkSigmaPOverP  ->push_back(-1.);
	      _phT1PhiReRecoNoTrkSigmaPOverP        ->push_back(-10.);
	      _phT1NHitsFittedReRecoNoTrkSigmaPOverP->push_back(0);
	      _phT1NHitsVetoReRecoNoTrkSigmaPOverP  ->push_back(0);
	      _phT2Chi2_NDOFReRecoNoTrkSigmaPOverP  ->push_back(-1.);
	      _phT2PhiReRecoNoTrkSigmaPOverP        ->push_back(-10.);
	      _phT2NHitsFittedReRecoNoTrkSigmaPOverP->push_back(0);
	      _phT2NHitsVetoReRecoNoTrkSigmaPOverP  ->push_back(0);
	    }else if(pandorapartReRecoNoTrkSigmaPOverP->getTracks().size()==2){
	      //converted photon: track values play a role now
	      _phT1Chi2_NDOFReRecoNoTrkSigmaPOverP  ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getChi2()/(float)pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getNdf());
	      _phT1PhiReRecoNoTrkSigmaPOverP        ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getPhi());
	      _phT1NHitsFittedReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[0]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[1]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[2]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[3]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[4]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	      _phT1NHitsVetoReRecoNoTrkSigmaPOverP  ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[6]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[7]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[8]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[9]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[10]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[11]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[0]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[1]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[2]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[3]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[4]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[0]->getSubdetectorHitNumbers()[5]);
	      _phT2Chi2_NDOFReRecoNoTrkSigmaPOverP  ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getChi2()/(float)pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getNdf());
	      _phT2PhiReRecoNoTrkSigmaPOverP        ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getPhi());
	      _phT2NHitsFittedReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[0]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[1]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[2]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[3]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[4]
					  +pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	      _phT2NHitsVetoReRecoNoTrkSigmaPOverP  ->push_back(pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[6]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[7]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[8]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[9]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[10]+
					  pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[11]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[0]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[1]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[2]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[3]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[4]
					  -pandorapartReRecoNoTrkSigmaPOverP->getTracks()[1]->getSubdetectorHitNumbers()[5]);
	      if(pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()==0){
		_phClusterEReRecoNoTrkSigmaPOverP->push_back(0.);
		//sudetector 0:ECAL
		_phECALReRecoNoTrkSigmaPOverP->push_back(0.);
		//sudetector 1:HCAL
		_phHCALReRecoNoTrkSigmaPOverP->push_back(0.);
		_phEoverTotReRecoNoTrkSigmaPOverP->push_back(-1.);
		//fraction of the first subcluster -> 0 cluster, obviously 0
		_phf1ClusterEReRecoNoTrkSigmaPOverP->push_back(0.);
		_phf1ECALReRecoNoTrkSigmaPOverP->push_back(0.);
		_phf1HCALReRecoNoTrkSigmaPOverP->push_back(0.);
	      }else if(pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()==2){
		//sudetector 0:ECAL
		_phECALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getSubdetectorEnergies()[0]);
		//sudetector 1:HCAL
		_phHCALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getSubdetectorEnergies()[1]);
		_phEoverTotReRecoNoTrkSigmaPOverP->push_back((pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getSubdetectorEnergies()[0])/(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy()+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getEnergy()));
		//fraction of the first subcluster -> 2 clusters, obviously NOT 100 % 
		_phf1ClusterEReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy()/(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getEnergy()+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getEnergy()));
		_phf1ECALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]/(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[0]+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getSubdetectorEnergies()[0]));
		_phf1HCALReRecoNoTrkSigmaPOverP->push_back(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[1]/(pandorapartReRecoNoTrkSigmaPOverP->getClusters()[0]->getSubdetectorEnergies()[1]+pandorapartReRecoNoTrkSigmaPOverP->getClusters()[1]->getSubdetectorEnergies()[1]));
	      }else if (pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()!=1){
		std::cout<<"converted photon, but unexpected size of clusters (neither 0,1,2) "<<pandorapartReRecoNoTrkSigmaPOverP->getClusters().size()<<std::endl;
	      }
	    }else{
	      std::cout<<"photon unusual track number: neither 0 (unconverted) nor 2 (conversion) "<<pandorapartReRecoNoTrkSigmaPOverP->getTracks().size()<<std::endl;
	    }
	    //fill energy contribution


	  }
	  //for isolation of interesting particles (i.e. photons ==22, electrons | |==11 and muons | |==13
	  if((abs(pandorapartReRecoNoTrkSigmaPOverP->getType())==13 && pandorapartReRecoNoTrkSigmaPOverP->getEnergy()>=_lepEMin) || (abs(pandorapartReRecoNoTrkSigmaPOverP->getType())==11 && pandorapartReRecoNoTrkSigmaPOverP->getEnergy()>=_lepEMin)|| (pandorapartReRecoNoTrkSigmaPOverP->getType()==22 && pandorapartReRecoNoTrkSigmaPOverP->getEnergy()>=_phEMin)){
	    float sum_rereco_DR03_chargedHadrons=0;
	    float sum_rereco_DR03_neutralHadrons=0;
	    float sum_rereco_DR03_photons=0;
	    float sum_rereco_DR03_muons=0;
	    float sum_rereco_DR03_electrons=0;
	    float sum_rereco_DR04_chargedHadrons=0;
	    float sum_rereco_DR04_neutralHadrons=0;
	    float sum_rereco_DR04_photons=0;
	    float sum_rereco_DR04_muons=0;
	    float sum_rereco_DR04_electrons=0;
	    float sum_rereco_cosTheta995_chargedHadrons=0;
	    float sum_rereco_cosTheta995_neutralHadrons=0;
	    float sum_rereco_cosTheta995_photons=0;
	    float sum_rereco_cosTheta995_muons=0;
	    float sum_rereco_cosTheta995_electrons=0;
	    for(int j=0;j<recoPartColReRecoNoTrkSigmaPOverP->getNumberOfElements();j++){
	      if(i!=j){//exclude particle from its own sum
		ReconstructedParticle* pandorapartReRecoNoTrkSigmaPOverP2 = dynamic_cast<ReconstructedParticle*>(recoPartColReRecoNoTrkSigmaPOverP->getElementAt(j));
		TLorentzVector temp2;
		temp2.SetPxPyPzE(pandorapartReRecoNoTrkSigmaPOverP2->getMomentum()[0],pandorapartReRecoNoTrkSigmaPOverP2->getMomentum()[1],pandorapartReRecoNoTrkSigmaPOverP2->getMomentum()[2],pandorapartReRecoNoTrkSigmaPOverP2->getEnergy());
		double deltaR=temp2.DeltaR(temp);
		double cosTheta=(temp.Px()*temp2.Px()+temp.Py()*temp2.Py()+temp.Pz()*temp2.Pz())/(temp.P()*temp2.P());
		if(cosTheta>0.995){
		  if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==211){
		    sum_rereco_cosTheta995_chargedHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		  }else if(pandorapartReRecoNoTrkSigmaPOverP2->getType()==22){
		    sum_rereco_cosTheta995_photons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==13){
		    sum_rereco_cosTheta995_muons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==11){
		    sum_rereco_cosTheta995_electrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())>300){
		    sum_rereco_cosTheta995_neutralHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		  }
		}
		if(deltaR<0.4){
		  if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==211){
		    sum_rereco_DR04_chargedHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_chargedHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    }	      
		  }else if(pandorapartReRecoNoTrkSigmaPOverP2->getType()==22){
		    sum_rereco_DR04_photons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_photons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==13){
		    sum_rereco_DR04_muons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_muons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())==11){
		    sum_rereco_DR04_electrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_electrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    }	      
		  }else if(abs(pandorapartReRecoNoTrkSigmaPOverP2->getType())>300){
		    //so far neutral I found are 310 for K_short, 2112 for neutrons and 3122 for Lambdas --> should include all neutrals
		    //obviously larger than leptons, charged pions and photons
		    sum_rereco_DR04_neutralHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    if(deltaR<0.3){
		      sum_rereco_DR03_neutralHadrons+=pandorapartReRecoNoTrkSigmaPOverP2->getEnergy();
		    }
		  }
		}
	      }
	    }//loop over second recopart collection for isolation calculation
	    if(abs(pandorapartReRecoNoTrkSigmaPOverP->getType())==11){
	      _elCHIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_chargedHadrons);
	      _elCHIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_chargedHadrons);
	      _elNHIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_neutralHadrons);
	      _elNHIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_neutralHadrons);
	      _elPhIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_photons);
	      _elPhIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_photons);
	      _elMuIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_muons);
	      _elMuIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_muons);
	      _elElIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_electrons);
	      _elElIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_electrons);
	      _elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_chargedHadrons);
	      _elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_neutralHadrons);
	      _elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_photons);
	      _elMuIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_muons);
	      _elElIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_electrons);
	    }else if(pandorapartReRecoNoTrkSigmaPOverP->getType()==22){
	      _phCHIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_chargedHadrons);
	      _phCHIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_chargedHadrons);
	      _phNHIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_neutralHadrons);
	      _phNHIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_neutralHadrons);
	      _phPhIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_photons);
	      _phPhIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_photons);
	      _phMuIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_muons);
	      _phMuIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_muons);
	      _phElIsoDR03ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR03_electrons);
	      _phElIsoDR04ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_DR04_electrons);
	      _phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_chargedHadrons);
	      _phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_neutralHadrons);
	      _phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_photons);
	      _phMuIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_muons);
	      _phElIsoCosTheta995ReRecoNoTrkSigmaPOverP->push_back(sum_rereco_cosTheta995_electrons);
	    }
	  }//recopart is either muon, electron or photon
	}
      }
    }

    _tree->Fill();
    //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !
    
    streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
			 << "   in run:  " << evt->getRunNumber() << std::endl ;
    
    
    
    _nEvt ++ ;
}



void MyDemoAnalyzer::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void MyDemoAnalyzer::end(){ 

    //   std::cout << "MyDemoAnalyzer::end()  " << name() 
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;

  if(_outputFile==NULL){
    std::cout<<"stupid file says it is NULL what the hell"<<std::endl;
  }

  _outputFile->Write();
  _outputFile->Close();

}

  // ---- method that builds the tree -------------------------------
void MyDemoAnalyzer::buildTree(){

  _mcpartStatus          = new std::vector<float>();
  _mcpartPx              = new std::vector<float>();
  _mcpartPy              = new std::vector<float>();
  _mcpartPz              = new std::vector<float>();
  _mcpartE               = new std::vector<float>();
  _mcpartPhi             = new std::vector<float>();
  _mcpartTheta           = new std::vector<float>();
  _mcpartMass            = new std::vector<float>();
  _mcpartPDGID           = new std::vector<int>();
  _mcpartParent1PDGID    = new std::vector<int>();
  _mcpartParent2PDGID    = new std::vector<int>();
  _mcpartDaughter1PDGID  = new std::vector<int>();
  _mcpartDaughter2PDGID  = new std::vector<int>();
  _mcpartDaughter3PDGID  = new std::vector<int>();

  //fill only status 1 here
  _elGENPx              = new std::vector<float>();
  _elGENPy              = new std::vector<float>();
  _elGENPz              = new std::vector<float>();
  _elGENE               = new std::vector<float>();
  _elGENPhi             = new std::vector<float>();
  _elGENTheta           = new std::vector<float>();
  _elGENIsoDR03         = new std::vector<float>();
  _elGENIsoDR04         = new std::vector<float>();
  _elGENIsoCosTheta995  = new std::vector<float>();
  _elGENPDGID           = new std::vector<int>();
  _elGENParent1PDGID    = new std::vector<int>();
  _elGENParent2PDGID    = new std::vector<int>();
  _elGENGrandParent1PDGID    = new std::vector<int>();

  _muGENPx              = new std::vector<float>();
  _muGENPy              = new std::vector<float>();
  _muGENPz              = new std::vector<float>();
  _muGENE               = new std::vector<float>();
  _muGENPhi             = new std::vector<float>();
  _muGENTheta           = new std::vector<float>();
  _muGENIsoDR03         = new std::vector<float>();
  _muGENIsoDR04         = new std::vector<float>();
  _muGENIsoCosTheta995  = new std::vector<float>();
  _muGENPDGID           = new std::vector<int>();
  _muGENParent1PDGID    = new std::vector<int>();
  _muGENParent2PDGID    = new std::vector<int>();
  _muGENGrandParent1PDGID    = new std::vector<int>();

  _phGENPx              = new std::vector<float>();
  _phGENPy              = new std::vector<float>();
  _phGENPz              = new std::vector<float>();
  _phGENE               = new std::vector<float>();
  _phGENPhi             = new std::vector<float>();
  _phGENTheta           = new std::vector<float>();
  _phGENIsoDR03         = new std::vector<float>();
  _phGENIsoDR04         = new std::vector<float>();
  _phGENIsoCosTheta995  = new std::vector<float>();
  _phGENParent1PDGID    = new std::vector<int>();
  _phGENParent2PDGID    = new std::vector<int>();
  _phGENGrandParent1PDGID = new std::vector<int>();

  _nuGENPx              = new std::vector<float>();
  _nuGENPy              = new std::vector<float>();
  _nuGENPz              = new std::vector<float>();
  _nuGENE               = new std::vector<float>();
  _nuGENPhi             = new std::vector<float>();
  _nuGENTheta           = new std::vector<float>();
  _nuGENPDGID           = new std::vector<int>();
  _nuGENParent1PDGID    = new std::vector<int>();
  _nuGENParent2PDGID    = new std::vector<int>();
  _nuGENGrandParent1PDGID    = new std::vector<int>();

  _jetPx        = new std::vector<float>(); 
  _jetPy        = new std::vector<float>(); 
  _jetPz        = new std::vector<float>(); 
  _jetE         = new std::vector<float>(); 
  _jetPhi       = new std::vector<float>(); 
  _jetTheta         = new std::vector<float>(); 
  _jetMass      = new std::vector<float>(); 
  _jetNF        = new std::vector<float>();
  _jetLF        = new std::vector<float>();
  _jetKF        = new std::vector<float>();
  _jetCHF        = new std::vector<float>();
  _jetPhF        = new std::vector<float>();
  _jetElF        = new std::vector<float>();
  _jetMuF        = new std::vector<float>();

  _jetNMult        = new std::vector<int>();
  _jetLMult        = new std::vector<int>();
  _jetKMult        = new std::vector<int>();
  _jetCHMult        = new std::vector<int>();
  _jetPhMult        = new std::vector<int>();
  _jetElMult        = new std::vector<int>();
  _jetMuMult        = new std::vector<int>();
  //std::cout<<"vector electron"<<std::endl;
  _elPx        = new std::vector<float>(); 
  _elPy        = new std::vector<float>(); 
  _elPz        = new std::vector<float>(); 
  _elE         = new std::vector<float>(); 
  _elPhi       = new std::vector<float>(); 
  _elTheta         = new std::vector<float>(); 
  _elMass      = new std::vector<float>();
  _elEoverTot  = new std::vector<float>();
  _elClusterE  = new std::vector<float>();
  _elECAL      = new std::vector<float>();
  _elHCAL      = new std::vector<float>();
  _el0D0        = new std::vector<float>();
  _el0Z0        = new std::vector<float>();
  _el0Phi       = new std::vector<float>();
  _el0Chi2_NDOF = new std::vector<float>();
  _el1D0        = new std::vector<float>();
  _el1Z0        = new std::vector<float>();
  _el1Phi       = new std::vector<float>();
  _el1Chi2_NDOF = new std::vector<float>();
  _elCHIsoDR03  = new std::vector<float>();
  _elNHIsoDR03  = new std::vector<float>();
  _elPhIsoDR03  = new std::vector<float>();
  _elCHIsoDR04  = new std::vector<float>();
  _elNHIsoDR04  = new std::vector<float>();
  _elPhIsoDR04  = new std::vector<float>();
  _elCHIsoCosTheta995  = new std::vector<float>();
  _elNHIsoCosTheta995  = new std::vector<float>();
  _elPhIsoCosTheta995  = new std::vector<float>();
  _elMuIsoCosTheta995  = new std::vector<float>();
  _elElIsoCosTheta995  = new std::vector<float>();
  //_elLongShowerProfile      = new std::vector<float>();
  //_elTransShowerProfile     = new std::vector<float>();
  //_elTrackBasedTransProfile = new std::vector<float>();
  _elEnergyEE               = new std::vector<float>();
  _elEnergyEB               = new std::vector<float>();
  _elEnergyEElse            = new std::vector<float>();
  _elEnergyHE               = new std::vector<float>();
  _elEnergyHB               = new std::vector<float>();
  _elEnergyHElse            = new std::vector<float>();
  _elEnergyElseB            = new std::vector<float>();
  _elEnergyElseE            = new std::vector<float>();
  _elEnergyElseElse         = new std::vector<float>();
  _elEnergyEM               = new std::vector<float>();
  _elEnergyHAD              = new std::vector<float>();
  _elEnergyMuon             = new std::vector<float>();
  //_elPeakEnergy             = new std::vector<float>();
  //_elNRadiationLength       = new std::vector<float>();
  _elMaxEInECALLayer        = new std::vector<float>();
  _elMaxEInHCALLayer        = new std::vector<float>();
  _elMinTrackClustDiff      = new std::vector<float>();
  _elSigmaPOverPTrack       = new std::vector<float>();
  _elFirstLayerECAL       = new std::vector<int>();
  _elLastLayerECAL        = new std::vector<int>();
  _elLayerMaxEECAL        = new std::vector<int>();
  _elLayerMaxHitsECAL     = new std::vector<int>();
  _elHitsECAL             = new std::vector<int>();
  _elHitsHCAL             = new std::vector<int>();
  _elMaxHitsPerLayerECAL  = new std::vector<int>();
  _elFirstLayerHCAL       = new std::vector<int>();
  _elLastLayerHCAL        = new std::vector<int>();
  _elMaxHitsPerLayerHCAL  = new std::vector<int>();
  _elNHitsVTX   = new std::vector<int>();
  _elNHitsFTD   = new std::vector<int>();
  _elNHitsSIT   = new std::vector<int>();
  _elNHitsTPC   = new std::vector<int>();
  _elNHitsSET   = new std::vector<int>();
  _elNHitsETD   = new std::vector<int>();
  _elNHitsVTXVeto   = new std::vector<int>();
  _elNHitsFTDVeto   = new std::vector<int>();
  _elNHitsSITVeto   = new std::vector<int>();
  _elNHitsTPCVeto   = new std::vector<int>();
  _elNHitsSETVeto   = new std::vector<int>();
  _elNHitsETDVeto   = new std::vector<int>();
  _el1NHits     = new std::vector<int>();
  _el1NHitsVeto = new std::vector<int>();
  _elID         = new std::vector<int>();
  _elTrackPt    = new std::vector<float>();
  _elTrackP     = new std::vector<float>();
  _elOmega      = new std::vector<float>();
  _elTrackClusterOpeningAngle      = new std::vector<float>(); 

  _muPx        = new std::vector<float>(); 
  _muPy        = new std::vector<float>(); 
  _muPz        = new std::vector<float>(); 
  _muE         = new std::vector<float>(); 
  _muPhi       = new std::vector<float>(); 
  _muTheta         = new std::vector<float>(); 
  _muMass      = new std::vector<float>();
  _muEoverTot  = new std::vector<float>();
  _muClusterE  = new std::vector<float>();
  _muECAL      = new std::vector<float>();
  _muHCAL      = new std::vector<float>();
  _mu0D0        = new std::vector<float>();
  _mu0Z0        = new std::vector<float>();
  _mu0Phi       = new std::vector<float>();
  _mu0Chi2_NDOF = new std::vector<float>();
  _mu1D0        = new std::vector<float>();
  _mu1Z0        = new std::vector<float>();
  _mu1Phi       = new std::vector<float>();
  _mu1Chi2_NDOF = new std::vector<float>();
  _muCHIsoDR03  = new std::vector<float>();
  _muNHIsoDR03  = new std::vector<float>();
  _muPhIsoDR03  = new std::vector<float>();
  _muMuIsoDR03  = new std::vector<float>();
  _muElIsoDR03  = new std::vector<float>();
  _muCHIsoDR04  = new std::vector<float>();
  _muNHIsoDR04  = new std::vector<float>();
  _muPhIsoDR04  = new std::vector<float>();
  _muMuIsoDR04  = new std::vector<float>();
  _muElIsoDR04  = new std::vector<float>();
  _muCHIsoCosTheta995  = new std::vector<float>();
  _muNHIsoCosTheta995  = new std::vector<float>();
  _muPhIsoCosTheta995  = new std::vector<float>();
  _muMuIsoCosTheta995  = new std::vector<float>();
  _muElIsoCosTheta995  = new std::vector<float>();
  _muNHitsVTX   = new std::vector<int>();
  _muNHitsFTD   = new std::vector<int>();
  _muNHitsSIT   = new std::vector<int>();
  _muNHitsTPC   = new std::vector<int>();
  _muNHitsSET   = new std::vector<int>();
  _muNHitsETD   = new std::vector<int>();
  _muNHitsVTXVeto   = new std::vector<int>();
  _muNHitsFTDVeto   = new std::vector<int>();
  _muNHitsSITVeto   = new std::vector<int>();
  _muNHitsTPCVeto   = new std::vector<int>();
  _muNHitsSETVeto   = new std::vector<int>();
  _muNHitsETDVeto   = new std::vector<int>();
  _mu1NHits     = new std::vector<int>();
  _mu1NHitsVeto = new std::vector<int>();
  _muID         = new std::vector<int>();
  _muTrackPt    = new std::vector<float>();
  _muTrackP     = new std::vector<float>();
  _muOmega      = new std::vector<float>();
  _muSigmaPOverPTrack       = new std::vector<float>();

  _phPx        = new std::vector<float>(); 
  _phPy        = new std::vector<float>(); 
  _phPz        = new std::vector<float>(); 
  _phE         = new std::vector<float>(); 
  _phPhi       = new std::vector<float>(); 
  _phTheta         = new std::vector<float>(); 
  _phMass      = new std::vector<float>();
  _phEoverTot  = new std::vector<float>();
  _phClusterE  = new std::vector<float>();
  _phECAL      = new std::vector<float>();
  _phHCAL      = new std::vector<float>();
  _phf1ClusterE  = new std::vector<float>();
  _phf1ECAL      = new std::vector<float>();
  _phf1HCAL      = new std::vector<float>();
  _phCHIsoDR03  = new std::vector<float>();
  _phNHIsoDR03  = new std::vector<float>();
  _phPhIsoDR03  = new std::vector<float>();
  _phMuIsoDR03  = new std::vector<float>();
  _phElIsoDR03  = new std::vector<float>();
  _phCHIsoDR04  = new std::vector<float>();
  _phNHIsoDR04  = new std::vector<float>();
  _phPhIsoDR04  = new std::vector<float>();
  _phMuIsoDR04  = new std::vector<float>();
  _phElIsoDR04  = new std::vector<float>();
  _phCHIsoCosTheta995  = new std::vector<float>();
  _phNHIsoCosTheta995  = new std::vector<float>();
  _phPhIsoCosTheta995  = new std::vector<float>();
  _phMuIsoCosTheta995  = new std::vector<float>();
  _phElIsoCosTheta995  = new std::vector<float>();
  //_phLongShowerProfile      = new std::vector<float>();
  //_phTransShowerProfile     = new std::vector<float>();
  //_phTrackBasedTransProfile = new std::vector<float>();
  _phEnergyEE               = new std::vector<float>();
  _phEnergyEB               = new std::vector<float>();
  _phEnergyEElse            = new std::vector<float>();
  _phEnergyHE               = new std::vector<float>();
  _phEnergyHB               = new std::vector<float>();
  _phEnergyHElse            = new std::vector<float>();
  _phEnergyElseB            = new std::vector<float>();
  _phEnergyElseE            = new std::vector<float>();
  _phEnergyElseElse         = new std::vector<float>();
  _phEnergyEM               = new std::vector<float>();
  _phEnergyHAD              = new std::vector<float>();
  _phEnergyMuon             = new std::vector<float>();
  //_phPeakEnergy             = new std::vector<float>();
  //_phNRadiationLength       = new std::vector<float>();
  _phMaxEInECALLayer        = new std::vector<float>();
  _phMaxEInHCALLayer        = new std::vector<float>();
  _phMinTrackClustDiff      = new std::vector<float>();
  _phSigmaPOverPClosestTrack= new std::vector<float>();
  _phClosestTrackPt         = new std::vector<float>();
  _phClosestTrackP          = new std::vector<float>();
  _phClosestTrackD0         = new std::vector<float>();
  _phClosestTrackZ0         = new std::vector<float>();
  _phClosestTrackPhi        = new std::vector<float>();
  _phClosestTrackChi2_NDOF  = new std::vector<float>();
  _phClosestTrackOmega      = new std::vector<float>();
  _phClosestTrackNHitsVTX   = new std::vector<int>();
  _phClosestTrackNHitsFTD   = new std::vector<int>();
  _phClosestTrackNHitsSIT   = new std::vector<int>();
  _phClosestTrackNHitsTPC   = new std::vector<int>();
  _phClosestTrackNHitsSET   = new std::vector<int>();
  _phClosestTrackNHitsETD   = new std::vector<int>();
  _phClosestTrackNHitsVTXVeto   = new std::vector<int>();
  _phClosestTrackNHitsFTDVeto   = new std::vector<int>();
  _phClosestTrackNHitsSITVeto   = new std::vector<int>();
  _phClosestTrackNHitsTPCVeto   = new std::vector<int>();
  _phClosestTrackNHitsSETVeto   = new std::vector<int>();
  _phClosestTrackNHitsETDVeto   = new std::vector<int>();
  _phMinSelTrackClustDiff      = new std::vector<float>();
  _phSigmaPOverPClosestSelTrack= new std::vector<float>();
  _phClosestSelTrackPt         = new std::vector<float>();
  _phClosestSelTrackP          = new std::vector<float>();
  _phFirstLayerECAL       = new std::vector<int>();
  _phLastLayerECAL        = new std::vector<int>();
  _phLayerMaxEECAL        = new std::vector<int>();
  _phLayerMaxHitsECAL     = new std::vector<int>();
  _phHitsECAL             = new std::vector<int>();
  _phHitsHCAL             = new std::vector<int>();
  _phMaxHitsPerLayerECAL  = new std::vector<int>();
  _phFirstLayerHCAL       = new std::vector<int>();
  _phLastLayerHCAL        = new std::vector<int>();
  _phMaxHitsPerLayerHCAL  = new std::vector<int>();
  _phT1Phi       = new std::vector<float>(); 
  _phT1Chi2_NDOF = new std::vector<float>();
  _phT2Phi       = new std::vector<float>(); 
  _phT2Chi2_NDOF = new std::vector<float>();
  _phT1NHitsFitted = new std::vector<int>();
  _phT1NHitsVeto  = new std::vector<int>();
  _phT2NHitsFitted = new std::vector<int>();
  _phT2NHitsVeto  = new std::vector<int>();
  _phClosestTrackClusterOpeningAngle      = new std::vector<float>(); 
  _phClosestSelTrackClusterOpeningAngle   = new std::vector<float>(); 

  if(_hasrereco_information){
    _jetPxReReco        = new std::vector<float>(); 
    _jetPyReReco        = new std::vector<float>(); 
    _jetPzReReco        = new std::vector<float>(); 
    _jetEReReco         = new std::vector<float>(); 
    _jetPhiReReco       = new std::vector<float>(); 
    _jetThetaReReco         = new std::vector<float>(); 
    _jetMassReReco      = new std::vector<float>(); 
    _jetNFReReco        = new std::vector<float>();
    _jetLFReReco        = new std::vector<float>();
    _jetKFReReco        = new std::vector<float>();
    _jetCHFReReco        = new std::vector<float>();
    _jetPhFReReco        = new std::vector<float>();
    _jetElFReReco        = new std::vector<float>();
    _jetMuFReReco        = new std::vector<float>();
    
    _jetNMultReReco        = new std::vector<int>();
    _jetLMultReReco        = new std::vector<int>();
    _jetKMultReReco        = new std::vector<int>();
    _jetCHMultReReco        = new std::vector<int>();
    _jetPhMultReReco        = new std::vector<int>();
    _jetElMultReReco        = new std::vector<int>();
    _jetMuMultReReco        = new std::vector<int>();
    
    _elPxReReco        = new std::vector<float>(); 
    _elPyReReco        = new std::vector<float>(); 
    _elPzReReco        = new std::vector<float>(); 
    _elEReReco         = new std::vector<float>(); 
    _elPhiReReco       = new std::vector<float>(); 
    _elThetaReReco     = new std::vector<float>(); 
    _elMassReReco      = new std::vector<float>();
    _elEoverTotReReco  = new std::vector<float>();
    _elClusterEReReco  = new std::vector<float>();
    _elECALReReco      = new std::vector<float>();
    _elHCALReReco      = new std::vector<float>();
    _el0D0ReReco        = new std::vector<float>();
    _el0Z0ReReco        = new std::vector<float>();
    _el0PhiReReco       = new std::vector<float>();
    _el0Chi2_NDOFReReco = new std::vector<float>();
    _el1D0ReReco        = new std::vector<float>();
    _el1Z0ReReco        = new std::vector<float>();
    _el1PhiReReco       = new std::vector<float>();
    _el1Chi2_NDOFReReco = new std::vector<float>();
    _elCHIsoDR03ReReco  = new std::vector<float>();
    _elNHIsoDR03ReReco  = new std::vector<float>();
    _elPhIsoDR03ReReco  = new std::vector<float>();
    _elMuIsoDR03ReReco  = new std::vector<float>();
    _elElIsoDR03ReReco  = new std::vector<float>();
    _elCHIsoDR04ReReco  = new std::vector<float>();
    _elNHIsoDR04ReReco  = new std::vector<float>();
    _elPhIsoDR04ReReco  = new std::vector<float>();
    _elMuIsoDR04ReReco  = new std::vector<float>();
    _elElIsoDR04ReReco  = new std::vector<float>();
    _elCHIsoCosTheta995ReReco  = new std::vector<float>();
    _elNHIsoCosTheta995ReReco  = new std::vector<float>();
    _elPhIsoCosTheta995ReReco  = new std::vector<float>();
    _elMuIsoCosTheta995ReReco  = new std::vector<float>();
    _elElIsoCosTheta995ReReco  = new std::vector<float>();
    _elTrackPtReReco    = new std::vector<float>();
    _elTrackPReReco     = new std::vector<float>();
    //_elLongShowerProfileReReco      = new std::vector<float>();
    //_elTransShowerProfileReReco     = new std::vector<float>();
    //_elTrackBasedTransProfileReReco = new std::vector<float>();
    _elEnergyEEReReco               = new std::vector<float>();
    _elEnergyEBReReco               = new std::vector<float>();
    _elEnergyEElseReReco            = new std::vector<float>();
    _elEnergyHEReReco               = new std::vector<float>();
    _elEnergyHBReReco               = new std::vector<float>();
    _elEnergyHElseReReco            = new std::vector<float>();
    _elEnergyElseBReReco            = new std::vector<float>();
    _elEnergyElseEReReco            = new std::vector<float>();
    _elEnergyElseElseReReco         = new std::vector<float>();
    _elEnergyEMReReco               = new std::vector<float>();
    _elEnergyHADReReco              = new std::vector<float>();
    _elEnergyMuonReReco             = new std::vector<float>();
    //_elPeakEnergyReReco             = new std::vector<float>();
    //_elNRadiationLengthReReco       = new std::vector<float>();
    _elMaxEInECALLayerReReco        = new std::vector<float>();
    _elMaxEInHCALLayerReReco        = new std::vector<float>();
    _elMinTrackClustDiffReReco      = new std::vector<float>();
    _elSigmaPOverPTrackReReco       = new std::vector<float>();
    _elFirstLayerECALReReco       = new std::vector<int>();
    _elLastLayerECALReReco        = new std::vector<int>();
    _elLayerMaxHitsECALReReco     = new std::vector<int>();
    _elLayerMaxEECALReReco        = new std::vector<int>();
    _elHitsECALReReco             = new std::vector<int>();
    _elHitsHCALReReco             = new std::vector<int>();
    _elMaxHitsPerLayerECALReReco  = new std::vector<int>();
    _elFirstLayerHCALReReco       = new std::vector<int>();
    _elLastLayerHCALReReco       = new std::vector<int>();
    _elMaxHitsPerLayerHCALReReco  = new std::vector<int>();
    _elNHitsVTXReReco   = new std::vector<int>();
    _elNHitsFTDReReco   = new std::vector<int>();
    _elNHitsSITReReco   = new std::vector<int>();
    _elNHitsTPCReReco   = new std::vector<int>();
    _elNHitsSETReReco   = new std::vector<int>();
    _elNHitsETDReReco   = new std::vector<int>();
    _elNHitsVTXVetoReReco   = new std::vector<int>();
    _elNHitsFTDVetoReReco   = new std::vector<int>();
    _elNHitsSITVetoReReco   = new std::vector<int>();
    _elNHitsTPCVetoReReco   = new std::vector<int>();
    _elNHitsSETVetoReReco   = new std::vector<int>();
    _elNHitsETDVetoReReco   = new std::vector<int>();
    _el1NHitsReReco     = new std::vector<int>();
    _el1NHitsVetoReReco = new std::vector<int>();
    _elIDReReco         = new std::vector<int>();
    _elOmegaReReco      = new std::vector<float>();
    _elTrackClusterOpeningAngleReReco   = new std::vector<float>(); 
    
    _muPxReReco        = new std::vector<float>(); 
    _muPyReReco        = new std::vector<float>(); 
    _muPzReReco        = new std::vector<float>(); 
    _muEReReco         = new std::vector<float>(); 
    _muPhiReReco       = new std::vector<float>(); 
    _muThetaReReco     = new std::vector<float>(); 
    _muMassReReco      = new std::vector<float>();
    _muEoverTotReReco  = new std::vector<float>();
    _muClusterEReReco  = new std::vector<float>();
    _muECALReReco      = new std::vector<float>();
    _muHCALReReco      = new std::vector<float>();
    _mu0D0ReReco        = new std::vector<float>();
    _mu0Z0ReReco        = new std::vector<float>();
    _mu0PhiReReco       = new std::vector<float>();
    _mu0Chi2_NDOFReReco = new std::vector<float>();
    _mu1D0ReReco        = new std::vector<float>();
    _mu1Z0ReReco        = new std::vector<float>();
    _mu1PhiReReco       = new std::vector<float>();
    _mu1Chi2_NDOFReReco = new std::vector<float>();
    _muCHIsoDR03ReReco  = new std::vector<float>();
    _muNHIsoDR03ReReco  = new std::vector<float>();
    _muPhIsoDR03ReReco  = new std::vector<float>();
    _muMuIsoDR03ReReco  = new std::vector<float>();
    _muElIsoDR03ReReco  = new std::vector<float>();
    _muCHIsoDR04ReReco  = new std::vector<float>();
    _muNHIsoDR04ReReco  = new std::vector<float>();
    _muPhIsoDR04ReReco  = new std::vector<float>();
    _muMuIsoDR04ReReco  = new std::vector<float>();
    _muElIsoDR04ReReco  = new std::vector<float>();
    _muCHIsoCosTheta995ReReco  = new std::vector<float>();
    _muNHIsoCosTheta995ReReco  = new std::vector<float>();
    _muPhIsoCosTheta995ReReco  = new std::vector<float>();
    _muMuIsoCosTheta995ReReco  = new std::vector<float>();
    _muElIsoCosTheta995ReReco  = new std::vector<float>();
    _muNHitsVTXReReco   = new std::vector<int>();
    _muNHitsFTDReReco   = new std::vector<int>();
    _muNHitsSITReReco   = new std::vector<int>();
    _muNHitsTPCReReco   = new std::vector<int>();
    _muNHitsSETReReco   = new std::vector<int>();
    _muNHitsETDReReco   = new std::vector<int>();
    _muNHitsVTXVetoReReco   = new std::vector<int>();
    _muNHitsFTDVetoReReco   = new std::vector<int>();
    _muNHitsSITVetoReReco   = new std::vector<int>();
    _muNHitsTPCVetoReReco   = new std::vector<int>();
    _muNHitsSETVetoReReco   = new std::vector<int>();
    _muNHitsETDVetoReReco   = new std::vector<int>();
    _mu1NHitsReReco     = new std::vector<int>();
    _mu1NHitsVetoReReco = new std::vector<int>();
    _muIDReReco         = new std::vector<int>();
    _muSigmaPOverPTrack = new std::vector<float>();
    _muOmegaReReco      = new std::vector<float>();
    _muTrackPtReReco    = new std::vector<float>();
    _muTrackPReReco     = new std::vector<float>();
    _muSigmaPOverPTrackReReco     = new std::vector<float>();
    
    _phPxReReco        = new std::vector<float>(); 
    _phPyReReco        = new std::vector<float>(); 
    _phPzReReco        = new std::vector<float>(); 
    _phEReReco         = new std::vector<float>(); 
    _phPhiReReco       = new std::vector<float>(); 
    _phThetaReReco     = new std::vector<float>(); 
    _phMassReReco      = new std::vector<float>();
    _phEoverTotReReco  = new std::vector<float>();
    _phClusterEReReco  = new std::vector<float>();
    _phECALReReco      = new std::vector<float>();
    _phHCALReReco      = new std::vector<float>();
    _phf1ClusterEReReco  = new std::vector<float>();
    _phf1ECALReReco      = new std::vector<float>();
    _phf1HCALReReco      = new std::vector<float>();
    _phCHIsoDR03ReReco  = new std::vector<float>();
    _phNHIsoDR03ReReco  = new std::vector<float>();
    _phPhIsoDR03ReReco  = new std::vector<float>();
    _phMuIsoDR03ReReco  = new std::vector<float>();
    _phElIsoDR03ReReco  = new std::vector<float>();
    _phCHIsoDR04ReReco  = new std::vector<float>();
    _phNHIsoDR04ReReco  = new std::vector<float>();
    _phPhIsoDR04ReReco  = new std::vector<float>();
    _phMuIsoDR04ReReco  = new std::vector<float>();
    _phElIsoDR04ReReco  = new std::vector<float>();
    _phCHIsoCosTheta995ReReco  = new std::vector<float>();
    _phNHIsoCosTheta995ReReco  = new std::vector<float>();
    _phPhIsoCosTheta995ReReco  = new std::vector<float>();
    _phMuIsoCosTheta995ReReco  = new std::vector<float>();
    _phElIsoCosTheta995ReReco  = new std::vector<float>();
    //_phLongShowerProfileReReco      = new std::vector<float>();
    //_phTransShowerProfileReReco     = new std::vector<float>();
    //_phTrackBasedTransProfileReReco = new std::vector<float>();
    _phEnergyEEReReco               = new std::vector<float>();
    _phEnergyEBReReco               = new std::vector<float>();
    _phEnergyEElseReReco            = new std::vector<float>();
    _phEnergyHEReReco               = new std::vector<float>();
    _phEnergyHBReReco               = new std::vector<float>();
    _phEnergyHElseReReco            = new std::vector<float>();
    _phEnergyElseBReReco            = new std::vector<float>();
    _phEnergyElseEReReco            = new std::vector<float>();
    _phEnergyElseElseReReco         = new std::vector<float>();
    _phEnergyEMReReco               = new std::vector<float>();
    _phEnergyHADReReco              = new std::vector<float>();
    _phEnergyMuonReReco             = new std::vector<float>();
    //_phPeakEnergyReReco             = new std::vector<float>();
    //_phNRadiationLengthReReco       = new std::vector<float>();
    _phMaxEInECALLayerReReco        = new std::vector<float>();
    _phMaxEInHCALLayerReReco        = new std::vector<float>();
    _phMinTrackClustDiffReReco      = new std::vector<float>();
    _phSigmaPOverPClosestTrackReReco= new std::vector<float>();
    _phClosestTrackPtReReco         = new std::vector<float>();
    _phClosestTrackPReReco          = new std::vector<float>();
    _phMinSelTrackClustDiffReReco      = new std::vector<float>();
    _phSigmaPOverPClosestSelTrackReReco= new std::vector<float>();
    _phClosestSelTrackPtReReco         = new std::vector<float>();
    _phClosestSelTrackPReReco          = new std::vector<float>();
    _phClosestTrackD0ReReco         = new std::vector<float>();
    _phClosestTrackZ0ReReco         = new std::vector<float>();
    _phClosestTrackPhiReReco        = new std::vector<float>();
    _phClosestTrackChi2_NDOFReReco  = new std::vector<float>();
    _phClosestTrackOmegaReReco      = new std::vector<float>();
    _phClosestTrackNHitsVTXReReco   = new std::vector<int>();
    _phClosestTrackNHitsFTDReReco   = new std::vector<int>();
    _phClosestTrackNHitsSITReReco   = new std::vector<int>();
    _phClosestTrackNHitsTPCReReco   = new std::vector<int>();
    _phClosestTrackNHitsSETReReco   = new std::vector<int>();
    _phClosestTrackNHitsETDReReco   = new std::vector<int>();
    _phClosestTrackNHitsVTXVetoReReco   = new std::vector<int>();
    _phClosestTrackNHitsFTDVetoReReco   = new std::vector<int>();
    _phClosestTrackNHitsSITVetoReReco   = new std::vector<int>();
    _phClosestTrackNHitsTPCVetoReReco   = new std::vector<int>();
    _phClosestTrackNHitsSETVetoReReco   = new std::vector<int>();
    _phClosestTrackNHitsETDVetoReReco   = new std::vector<int>();
    _phFirstLayerECALReReco       = new std::vector<int>();
    _phLastLayerECALReReco        = new std::vector<int>();
    _phLayerMaxHitsECALReReco     = new std::vector<int>();
    _phLayerMaxEECALReReco       = new std::vector<int>();
    _phHitsECALReReco             = new std::vector<int>();
    _phHitsHCALReReco             = new std::vector<int>();
    _phMaxHitsPerLayerECALReReco  = new std::vector<int>();
    _phFirstLayerHCALReReco       = new std::vector<int>();
    _phLastLayerHCALReReco        = new std::vector<int>();
    _phMaxHitsPerLayerHCALReReco  = new std::vector<int>();
    _phT1PhiReReco       = new std::vector<float>(); 
    _phT1Chi2_NDOFReReco = new std::vector<float>();
    _phT2PhiReReco       = new std::vector<float>(); 
    _phT2Chi2_NDOFReReco = new std::vector<float>();
    _phT1NHitsFittedReReco = new std::vector<int>();
    _phT1NHitsVetoReReco  = new std::vector<int>();
    _phT2NHitsFittedReReco = new std::vector<int>();
    _phT2NHitsVetoReReco  = new std::vector<int>();
    _phClosestTrackClusterOpeningAngleReReco   = new std::vector<float>(); 
    _phClosestSelTrackClusterOpeningAngleReReco   = new std::vector<float>(); 
  }
  if(_hasrerecoNoTrackSigmaPOverP_information){
    _elPxReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _elPyReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _elPzReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _elEReRecoNoTrkSigmaPOverP         = new std::vector<float>(); 
    _elPhiReRecoNoTrkSigmaPOverP       = new std::vector<float>(); 
    _elThetaReRecoNoTrkSigmaPOverP     = new std::vector<float>(); 
    _elMassReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _elEoverTotReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elClusterEReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elECALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _elHCALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _el0D0ReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _el0Z0ReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _el0PhiReRecoNoTrkSigmaPOverP       = new std::vector<float>();
    _el0Chi2_NDOFReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _el1D0ReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _el1Z0ReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _el1PhiReRecoNoTrkSigmaPOverP       = new std::vector<float>();
    _el1Chi2_NDOFReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _elCHIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elNHIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elPhIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elMuIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elElIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elCHIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elNHIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elPhIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elMuIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elElIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elMuIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elElIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _elTrackPtReRecoNoTrkSigmaPOverP    = new std::vector<float>();
    _elTrackPReRecoNoTrkSigmaPOverP     = new std::vector<float>();
    //_elLongShowerProfileReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    //_elTransShowerProfileReRecoNoTrkSigmaPOverP     = new std::vector<float>();
    //_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _elEnergyEEReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _elEnergyEBReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _elEnergyEElseReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _elEnergyHEReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _elEnergyHBReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _elEnergyHElseReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _elEnergyElseBReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _elEnergyElseEReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _elEnergyElseElseReRecoNoTrkSigmaPOverP         = new std::vector<float>();
    _elEnergyEMReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _elEnergyHADReRecoNoTrkSigmaPOverP              = new std::vector<float>();
    _elEnergyMuonReRecoNoTrkSigmaPOverP             = new std::vector<float>();
    //_elPeakEnergyReRecoNoTrkSigmaPOverP             = new std::vector<float>();
    //_elNRadiationLengthReRecoNoTrkSigmaPOverP       = new std::vector<float>();
    _elMaxEInECALLayerReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _elMaxEInHCALLayerReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _elMinTrackClustDiffReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _elSigmaPOverPTrackReRecoNoTrkSigmaPOverP       = new std::vector<float>();
    _elFirstLayerECALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _elLastLayerECALReRecoNoTrkSigmaPOverP        = new std::vector<int>();
    _elLayerMaxHitsECALReRecoNoTrkSigmaPOverP     = new std::vector<int>();
    _elLayerMaxEECALReRecoNoTrkSigmaPOverP        = new std::vector<int>();
    _elHitsECALReRecoNoTrkSigmaPOverP             = new std::vector<int>();
    _elHitsHCALReRecoNoTrkSigmaPOverP             = new std::vector<int>();
    _elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP  = new std::vector<int>();
    _elFirstLayerHCALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _elLastLayerHCALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP  = new std::vector<int>();
    _elNHitsVTXReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsFTDReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsSITReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsTPCReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsSETReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsETDReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsVTXVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsFTDVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsSITVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsTPCVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsSETVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _elNHitsETDVetoReRecoNoTrkSigmaPOverP   = new std::vector<int>();
    _el1NHitsReRecoNoTrkSigmaPOverP     = new std::vector<int>();
    _el1NHitsVetoReRecoNoTrkSigmaPOverP = new std::vector<int>();
    _elIDReRecoNoTrkSigmaPOverP         = new std::vector<int>();
    _elOmegaReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    
    _phPxReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _phPyReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _phPzReRecoNoTrkSigmaPOverP        = new std::vector<float>(); 
    _phEReRecoNoTrkSigmaPOverP         = new std::vector<float>(); 
    _phPhiReRecoNoTrkSigmaPOverP       = new std::vector<float>(); 
    _phThetaReRecoNoTrkSigmaPOverP     = new std::vector<float>(); 
    _phMassReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phEoverTotReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phClusterEReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phECALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phHCALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phf1ClusterEReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phf1ECALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phf1HCALReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phCHIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phNHIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phPhIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phMuIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phElIsoDR03ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phCHIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phNHIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phPhIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phMuIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phElIsoDR04ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phMuIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    _phElIsoCosTheta995ReRecoNoTrkSigmaPOverP  = new std::vector<float>();
    //_phLongShowerProfileReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    //_phTransShowerProfileReRecoNoTrkSigmaPOverP     = new std::vector<float>();
    //_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _phEnergyEEReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _phEnergyEBReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _phEnergyEElseReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _phEnergyHEReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _phEnergyHBReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _phEnergyHElseReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _phEnergyElseBReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _phEnergyElseEReRecoNoTrkSigmaPOverP            = new std::vector<float>();
    _phEnergyElseElseReRecoNoTrkSigmaPOverP         = new std::vector<float>();
    _phEnergyEMReRecoNoTrkSigmaPOverP               = new std::vector<float>();
    _phEnergyHADReRecoNoTrkSigmaPOverP              = new std::vector<float>();
    _phEnergyMuonReRecoNoTrkSigmaPOverP             = new std::vector<float>();
    //_phPeakEnergyReRecoNoTrkSigmaPOverP             = new std::vector<float>();
    //_phNRadiationLengthReRecoNoTrkSigmaPOverP       = new std::vector<float>();
    _phMaxEInECALLayerReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _phMaxEInHCALLayerReRecoNoTrkSigmaPOverP        = new std::vector<float>();
    _phMinTrackClustDiffReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP= new std::vector<float>();
    _phClosestTrackPtReRecoNoTrkSigmaPOverP         = new std::vector<float>();
    _phClosestTrackPReRecoNoTrkSigmaPOverP          = new std::vector<float>();
    _phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP      = new std::vector<float>();
    _phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP= new std::vector<float>();
    _phClosestSelTrackPtReRecoNoTrkSigmaPOverP         = new std::vector<float>();
    _phClosestSelTrackPReRecoNoTrkSigmaPOverP          = new std::vector<float>();
    _phFirstLayerECALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _phLastLayerECALReRecoNoTrkSigmaPOverP        = new std::vector<int>();
    _phLayerMaxHitsECALReRecoNoTrkSigmaPOverP     = new std::vector<int>();
    _phLayerMaxEECALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _phHitsECALReRecoNoTrkSigmaPOverP             = new std::vector<int>();
    _phHitsHCALReRecoNoTrkSigmaPOverP             = new std::vector<int>();
    _phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP  = new std::vector<int>();
    _phFirstLayerHCALReRecoNoTrkSigmaPOverP       = new std::vector<int>();
    _phLastLayerHCALReRecoNoTrkSigmaPOverP        = new std::vector<int>();
    _phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP  = new std::vector<int>();
    _phT1PhiReRecoNoTrkSigmaPOverP       = new std::vector<float>(); 
    _phT1Chi2_NDOFReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _phT2PhiReRecoNoTrkSigmaPOverP       = new std::vector<float>(); 
    _phT2Chi2_NDOFReRecoNoTrkSigmaPOverP = new std::vector<float>();
    _phT1NHitsFittedReRecoNoTrkSigmaPOverP = new std::vector<int>();
    _phT1NHitsVetoReRecoNoTrkSigmaPOverP  = new std::vector<int>();
    _phT2NHitsFittedReRecoNoTrkSigmaPOverP = new std::vector<int>();
    _phT2NHitsVetoReRecoNoTrkSigmaPOverP  = new std::vector<int>();
  }

  _tree->Branch("runNum",  &_numRun,"runNum/I");
  _tree->Branch("eventNum",&_numEvt,"eventNum/I");
  _tree->Branch("weight",  &_weight,"weight/D");
  _tree->Branch("MEx",     &_MEx,   "MEx/F");
  _tree->Branch("MEy",     &_MEy,   "MEy/F");
  _tree->Branch("MEz",     &_MEz,   "MEz/F");
  _tree->Branch("MET",     &_MET,   "MET/F");//actually MPt
  _tree->Branch("SumE",    &_SumE,  "SumE/F");
  _tree->Branch("SumPt",   &_SumPt, "SumPt/F");
  _tree->Branch("SumChargedHadronE",   &_SumChargedHadronE, "SumChargedHadronE/F");
  _tree->Branch("SumChargedHadronPt",  &_SumChargedHadronPt,"SumChargedHadronPt/F");
  _tree->Branch("SumNeutralHadronE",   &_SumNeutralHadronE, "SumNeutralHadronE/F");
  _tree->Branch("SumNeutralHadronPt",  &_SumNeutralHadronPt,"SumNeutralHadronPt/F");
  _tree->Branch("SumPhotonE",          &_SumPhotonE,   "SumPhotonE/F");
  _tree->Branch("SumPhotonPt",         &_SumPhotonPt,  "SumPhotonPt/F");
  _tree->Branch("SumMuonE",            &_SumMuonE,     "SumMuonE/F");
  _tree->Branch("SumMuonPt",           &_SumMuonPt,    "SumMuonPt/F");
  _tree->Branch("SumElectronE",        &_SumElectronE, "SumElectronE/F");
  _tree->Branch("SumElectronPt",       &_SumElectronPt,"SumElectronPt/F");
  _tree->Branch("MExGEN",  &_MExGEN,   "MExGEN/F");
  _tree->Branch("MEyGEN",  &_MEyGEN,   "MEyGEN/F");
  _tree->Branch("MEzGEN",  &_MEzGEN,   "MEzGEN/F");
  _tree->Branch("METGEN",  &_METGEN,   "METGEN/F");//actually MPt
  _tree->Branch("MEGEN",   &_MEGEN, "MEGEN/F");//actually MPt
  _tree->Branch("SumEGEN",    &_SumEGEN,  "SumEGEN/F");
  _tree->Branch("SumPtGEN",   &_SumPtGEN, "SumPtGEN/F");
  _tree->Branch("SumChargedHadronEGEN",   &_SumChargedHadronEGEN, "SumChargedHadronEGEN/F");
  _tree->Branch("SumChargedHadronPtGEN",  &_SumChargedHadronPtGEN,"SumChargedHadronPtGEN/F");
  _tree->Branch("SumNeutralHadronEGEN",   &_SumNeutralHadronEGEN, "SumNeutralHadronEGEN/F");
  _tree->Branch("SumNeutralHadronPtGEN",  &_SumNeutralHadronPtGEN,"SumNeutralHadronPtGEN/F");
  _tree->Branch("SumPhotonEGEN",          &_SumPhotonEGEN,   "SumPhotonEGEN/F");
  _tree->Branch("SumPhotonPtGEN",         &_SumPhotonPtGEN,  "SumPhotonPtGEN/F");
  _tree->Branch("SumMuonEGEN",            &_SumMuonEGEN,     "SumMuonEGEN/F");
  _tree->Branch("SumMuonPtGEN",           &_SumMuonPtGEN,    "SumMuonPtGEN/F");
  _tree->Branch("SumElectronEGEN",        &_SumElectronEGEN, "SumElectronEGEN/F");
  _tree->Branch("SumElectronPtGEN",       &_SumElectronPtGEN,"SumElectronPtGEN/F");

  //after all isolation against WHAT?
  _tree->Branch("mcpartStatus",          "vector<float>",       &_mcpartStatus);
  _tree->Branch("mcpartPx",              "vector<float>",       &_mcpartPx);
  _tree->Branch("mcpartPy",              "vector<float>",       &_mcpartPy);
  _tree->Branch("mcpartPz",              "vector<float>",       &_mcpartPz);
  _tree->Branch("mcpartE",               "vector<float>",       &_mcpartE);
  _tree->Branch("mcpartPhi",             "vector<float>",       &_mcpartPhi);
  _tree->Branch("mcpartTheta",           "vector<float>",       &_mcpartTheta);
  _tree->Branch("mcpartMass",            "vector<float>",       &_mcpartMass);
  _tree->Branch("mcpartPDGID",           "vector<int>",         &_mcpartPDGID);
  _tree->Branch("mcpartParent1PDGID",    "vector<int>",         &_mcpartParent1PDGID);
  _tree->Branch("mcpartParent2PDGID",    "vector<int>",         &_mcpartParent2PDGID);
  _tree->Branch("mcpartDaughter1PDGID",  "vector<int>",         &_mcpartDaughter1PDGID);
  _tree->Branch("mcpartDaughter2PDGID",  "vector<int>",         &_mcpartDaughter2PDGID);
  _tree->Branch("mcpartDaughter3PDGID",  "vector<int>",         &_mcpartDaughter3PDGID);

  //fill only status 1 here
  _tree->Branch("elGENPx",              "vector<float>",        &_elGENPx);
  _tree->Branch("elGENPy",              "vector<float>",        &_elGENPy);
  _tree->Branch("elGENPz",              "vector<float>",        &_elGENPz);
  _tree->Branch("elGENE",               "vector<float>",        &_elGENE);
  _tree->Branch("elGENPhi",             "vector<float>",        &_elGENPhi);
  _tree->Branch("elGENTheta",           "vector<float>",        &_elGENTheta);
  _tree->Branch("elGENIsoDR03",         "vector<float>",        &_elGENIsoDR03);
  _tree->Branch("elGENIsoDR04",         "vector<float>",        &_elGENIsoDR04);
  _tree->Branch("elGENIsoCosTheta995",  "vector<float>",        &_elGENIsoCosTheta995);
  _tree->Branch("elGENPDGID",           "vector<int>",          &_elGENPDGID);
  _tree->Branch("elGENParent1PDGID",    "vector<int>",          &_elGENParent1PDGID);
  _tree->Branch("elGENParent2PDGID",    "vector<int>",          &_elGENParent2PDGID);
  _tree->Branch("elGENGrandParent1PDGID","vector<int>",         &_elGENGrandParent1PDGID);

  _tree->Branch("muGENPx",              "vector<float>",        &_muGENPx);
  _tree->Branch("muGENPy",              "vector<float>",        &_muGENPz);
  _tree->Branch("muGENPz",              "vector<float>",        &_muGENE);
  _tree->Branch("muGENE",               "vector<float>",        &_muGENPhi);
  _tree->Branch("muGENPhi",             "vector<float>",        &_muGENPhi);
  _tree->Branch("muGENTheta",           "vector<float>",        &_muGENTheta);
  _tree->Branch("muGENIsoDR03",         "vector<float>",        &_muGENIsoDR03);
  _tree->Branch("muGENIsoDR04",         "vector<float>",        &_muGENIsoDR04);
  _tree->Branch("muGENIsoCosTheta995",  "vector<float>",        &_muGENIsoCosTheta995);
  _tree->Branch("muGENPDGID",           "vector<int>",          &_muGENPDGID);
  _tree->Branch("muGENParent1PDGID",    "vector<int>",          &_muGENParent1PDGID);
  _tree->Branch("muGENParent2PDGID",    "vector<int>",          &_muGENParent2PDGID);
  _tree->Branch("muGENGrandParent1PDGID","vector<int>",         &_muGENGrandParent1PDGID);

  _tree->Branch("phGENPx",              "vector<float>",        &_phGENPx);
  _tree->Branch("phGENPy",              "vector<float>",        &_phGENPy);
  _tree->Branch("phGENPz",              "vector<float>",        &_phGENPz);
  _tree->Branch("phGENE",               "vector<float>",        &_phGENE);
  _tree->Branch("phGENPhi",             "vector<float>",        &_phGENPhi);
  _tree->Branch("phGENTheta",           "vector<float>",        &_phGENTheta);
  _tree->Branch("phGENIsoDR03",         "vector<float>",        &_phGENIsoDR03);
  _tree->Branch("phGENIsoDR04",         "vector<float>",        &_phGENIsoDR04);
  _tree->Branch("phGENIsoCosTheta995",  "vector<float>",        &_phGENIsoCosTheta995);
  _tree->Branch("phGENParent1PDGID",    "vector<int>",          &_phGENParent1PDGID);
  _tree->Branch("phGENParent2PDGID",    "vector<int>",          &_phGENParent2PDGID);
  _tree->Branch("phGENGrandParent1PDGID","vector<int>",         &_phGENGrandParent1PDGID);

  _tree->Branch("nuGENPx",              "vector<float>",        &_nuGENPx);
  _tree->Branch("nuGENPy",              "vector<float>",        &_nuGENPy);
  _tree->Branch("nuGENPz",              "vector<float>",        &_nuGENPz);
  _tree->Branch("nuGENE",               "vector<float>",        &_nuGENE);
  _tree->Branch("nuGENPhi",             "vector<float>",        &_nuGENPhi);
  _tree->Branch("nuGENTheta",           "vector<float>",        &_nuGENTheta);
  _tree->Branch("nuGENPDGID",           "vector<int>",          &_nuGENPDGID);
  _tree->Branch("nuGENParent1PDGID",    "vector<int>",          &_nuGENParent1PDGID);
  _tree->Branch("nuGENParent2PDGID",    "vector<int>",          &_nuGENParent2PDGID);
  _tree->Branch("nuGENGrandParent1PDGID","vector<int>",         &_nuGENGrandParent1PDGID);

  _tree->Branch("jetPx",   "vector<float>",   &_jetPx);
  _tree->Branch("jetPy",   "vector<float>",   &_jetPy);
  _tree->Branch("jetPz",   "vector<float>",   &_jetPz);
  _tree->Branch("jetE",    "vector<float>",   &_jetE);
  _tree->Branch("jetPhi",  "vector<float>",   &_jetPhi);
  _tree->Branch("jetTheta", "vector<float>",  &_jetTheta);
  _tree->Branch("jetMass", "vector<float>",   &_jetMass);
  _tree->Branch("jetNF",   "vector<float>",   &_jetNF);     
  _tree->Branch("jetLF",   "vector<float>",   &_jetLF);   
  _tree->Branch("jetKF",   "vector<float>",   &_jetKF);   
  _tree->Branch("jetCHF",  "vector<float>",   &_jetCHF);  
  _tree->Branch("jetPhF",  "vector<float>",   &_jetPhF);  
  _tree->Branch("jetMuF",  "vector<float>",   &_jetElF);  
  _tree->Branch("jetMuF",  "vector<float>",   &_jetMuF);
  _tree->Branch("jetNMult","vector<int>",     &_jetNMult);     
  _tree->Branch("jetLMult","vector<int>",     &_jetLMult);   
  _tree->Branch("jetKMult","vector<int>",     &_jetKMult);   
  _tree->Branch("jetCHMult","vector<int>",    &_jetCHMult);  
  _tree->Branch("jetPhMult","vector<int>",    &_jetPhMult);  
  _tree->Branch("jetElMult","vector<int>",    &_jetElMult);  
  _tree->Branch("jetMuMult","vector<int>",    &_jetMuMult);
   
  _tree->Branch("elPx",  "vector<float>",   &_elPx); 
  _tree->Branch("elPy",  "vector<float>",   &_elPy); 
  _tree->Branch("elPz",  "vector<float>",   &_elPz); 
  _tree->Branch("elE",   "vector<float>",   &_elE); 
  _tree->Branch("elPhi", "vector<float>",   &_elPhi); 
  _tree->Branch("elTheta","vector<float>",  &_elTheta); 
  _tree->Branch("elMass","vector<float>",   &_elMass);
  _tree->Branch("elClusterE",  "vector<float>",   &_elClusterE);
  _tree->Branch("elEoverTot",  "vector<float>",   &_elEoverTot);
  _tree->Branch("elECAL",      "vector<float>",   &_elECAL);
  _tree->Branch("elHCAL",      "vector<float>",   &_elHCAL);
  _tree->Branch("el0D0",   "vector<float>",  &_el0D0); 
  _tree->Branch("el0Z0",   "vector<float>",  &_el0Z0); 
  _tree->Branch("el0Phi",  "vector<float>",  &_el0Phi); 
  _tree->Branch("el0Chi2_NDOF","vector<float>",&_el0Chi2_NDOF);
  _tree->Branch("el1D0",   "vector<float>",  &_el1D0); 
  _tree->Branch("el1Z0",   "vector<float>",  &_el1Z0); 
  _tree->Branch("el1Phi",  "vector<float>",  &_el1Phi); 
  _tree->Branch("el1Chi2_NDOF","vector<float>",&_el1Chi2_NDOF);
  _tree->Branch("elCHIsoDR03",  "vector<float>",   &_elCHIsoDR03);
  _tree->Branch("elNHIsoDR03",  "vector<float>",   &_elNHIsoDR03);
  _tree->Branch("elPhIsoDR03",  "vector<float>",   &_elPhIsoDR03);
  _tree->Branch("elMuIsoDR03",  "vector<float>",   &_elMuIsoDR03);
  _tree->Branch("elElIsoDR03",  "vector<float>",   &_elElIsoDR03);
  _tree->Branch("elCHIsoDR04",  "vector<float>",   &_elCHIsoDR04);
  _tree->Branch("elNHIsoDR04",  "vector<float>",   &_elNHIsoDR04);
  _tree->Branch("elPhIsoDR04",  "vector<float>",   &_elPhIsoDR04);
  _tree->Branch("elMuIsoDR04",  "vector<float>",   &_elMuIsoDR04);
  _tree->Branch("elElIsoDR04",  "vector<float>",   &_elElIsoDR04);
  _tree->Branch("elNHitsVTX",   "vector<int>",     &_elNHitsVTX);
  _tree->Branch("elNHitsFTD",   "vector<int>",     &_elNHitsFTD);
  _tree->Branch("elNHitsSIT",   "vector<int>",     &_elNHitsSIT);
  _tree->Branch("elNHitsTPC",   "vector<int>",     &_elNHitsTPC);
  _tree->Branch("elNHitsSET",   "vector<int>",     &_elNHitsSET);
  _tree->Branch("elNHitsETD",   "vector<int>",     &_elNHitsETD);
  _tree->Branch("elNHitsVTXVeto",   "vector<int>",     &_elNHitsVTXVeto);
  _tree->Branch("elNHitsFTDVeto",   "vector<int>",     &_elNHitsFTDVeto);
  _tree->Branch("elNHitsSITVeto",   "vector<int>",     &_elNHitsSITVeto);
  _tree->Branch("elNHitsTPCVeto",   "vector<int>",     &_elNHitsTPCVeto);
  _tree->Branch("elNHitsSETVeto",   "vector<int>",     &_elNHitsSETVeto);
  _tree->Branch("elNHitsETDVeto",   "vector<int>",     &_elNHitsETDVeto);
  _tree->Branch("el1NHit",       "vector<int>",     &_el1NHits);
  _tree->Branch("el1NHitsVeto",  "vector<int>",     &_el1NHitsVeto);
  _tree->Branch("elCHIsoCosTheta995",  "vector<float>",   &_elCHIsoCosTheta995);
  _tree->Branch("elNHIsoCosTheta995",  "vector<float>",   &_elNHIsoCosTheta995);
  _tree->Branch("elPhIsoCosTheta995",  "vector<float>",   &_elPhIsoCosTheta995);
  _tree->Branch("elMuIsoCosTheta995",  "vector<float>",   &_elMuIsoCosTheta995);
  _tree->Branch("elElIsoCosTheta995",  "vector<float>",   &_elElIsoCosTheta995);
  //_tree->Branch("elLongShowerProfile",      "vector<float>",    &_elLongShowerProfile);
  //_tree->Branch("elTransShowerProfile",     "vector<float>",    &_elTransShowerProfile);
  //_tree->Branch("elTrackBasedTransProfile", "vector<float>",    &_elTrackBasedTransProfile);
  _tree->Branch("elTrackPt",                "vector<float>",    &_elTrackPt);
  _tree->Branch("elTrackP",                 "vector<float>",    &_elTrackP);
  _tree->Branch("elEnergyEE",               "vector<float>",    &_elEnergyEE);
  _tree->Branch("elEnergyEB",               "vector<float>",    &_elEnergyEB);
  _tree->Branch("elEnergyEElse",            "vector<float>",    &_elEnergyEElse);
  _tree->Branch("elEnergyHE",               "vector<float>",    &_elEnergyHE);
  _tree->Branch("elEnergyHB",               "vector<float>",    &_elEnergyHB);
  _tree->Branch("elEnergyHElse",            "vector<float>",    &_elEnergyHElse);
  _tree->Branch("elEnergyElseB",            "vector<float>",    &_elEnergyElseB);
  _tree->Branch("elEnergyElseE",            "vector<float>",    &_elEnergyElseE);
  _tree->Branch("elEnergyElseElse",         "vector<float>",    &_elEnergyElseElse);
  _tree->Branch("elEnergyEM",               "vector<float>",    &_elEnergyEM);
  _tree->Branch("elEnergyHAD",              "vector<float>",    &_elEnergyHAD);
  _tree->Branch("elEnergyMuon",             "vector<float>",    &_elEnergyMuon);
  //_tree->Branch("elPeakEnergy",             "vector<float>",    &_elPeakEnergy);
  //_tree->Branch("elNRadiationLength",       "vector<float>",    &_elNRadiationLength);
  _tree->Branch("elMaxEInECALLayer",        "vector<float>",    &_elMaxEInECALLayer);
  _tree->Branch("elMaxEInHCALLayer",        "vector<float>",    &_elMaxEInHCALLayer);
  _tree->Branch("elMinTrackClustDiff",      "vector<float>",    &_elMinTrackClustDiff);
  _tree->Branch("elSigmaPOverPTrack",       "vector<float>",    &_elSigmaPOverPTrack);
  _tree->Branch("elFirstLayerECAL",       "vector<int>",    &_elFirstLayerECAL);
  _tree->Branch("elLastLayerECAL",        "vector<int>",    &_elLastLayerECAL);
  _tree->Branch("elLayerMaxEECAL",        "vector<int>",    &_elLayerMaxEECAL);
  _tree->Branch("elLayerMaxHitsECAL",     "vector<int>",    &_elLayerMaxHitsECAL);
  _tree->Branch("elHitsECAL",             "vector<int>",    &_elHitsECAL);
  _tree->Branch("elHitsHCAL",             "vector<int>",    &_elHitsHCAL);
  _tree->Branch("elMaxHitsPerLayerECAL",  "vector<int>",    &_elMaxHitsPerLayerECAL);
  _tree->Branch("elFirstLayerHCAL",       "vector<int>",    &_elFirstLayerHCAL);
  _tree->Branch("elLastLayerHCAL",        "vector<int>",    &_elFirstLayerHCAL);
  _tree->Branch("elMaxHitsPerLayerHCAL",  "vector<int>",    &_elMaxHitsPerLayerHCAL);
  _tree->Branch("elID",  "vector<int>",   &_elID);
  _tree->Branch("elOmega","vector<float>",&_elOmega);
  _tree->Branch("elTrackClusterOpeningAngle" ,"vector<float>",  &_elTrackClusterOpeningAngle);

  _tree->Branch("muPx",  "vector<float>",   &_muPx); 
  _tree->Branch("muPy",  "vector<float>",   &_muPy); 
  _tree->Branch("muPz",  "vector<float>",   &_muPz); 
  _tree->Branch("muE",  "vector<float>",    &_muE); 
  _tree->Branch("muPhi",  "vector<float>",  &_muPhi); 
  _tree->Branch("muTheta","vector<float>",  &_muTheta); 
  _tree->Branch("muMass",  "vector<float>", &_muMass);
  _tree->Branch("muEoverTot",  "vector<float>",   &_muEoverTot);
  _tree->Branch("muClusterE",  "vector<float>",   &_muClusterE);
  _tree->Branch("muECAL",      "vector<float>",   &_muECAL);
  _tree->Branch("muHCAL",      "vector<float>",   &_muHCAL);
  _tree->Branch("mu0D0",   "vector<float>",  &_mu0D0); 
  _tree->Branch("mu0Z0",   "vector<float>",  &_mu0Z0); 
  _tree->Branch("mu0Phi",  "vector<float>",  &_mu0Phi); 
  _tree->Branch("mu0Chi2_NDOF","vector<float>",&_mu0Chi2_NDOF);
  _tree->Branch("muD0",   "vector<float>",  &_mu1D0); 
  _tree->Branch("mu1Z0",   "vector<float>",  &_mu1Z0); 
  _tree->Branch("mu1Phi",  "vector<float>",  &_mu1Phi); 
  _tree->Branch("mu1Chi2_NDOF","vector<float>",&_mu1Chi2_NDOF);
  _tree->Branch("muCHIsoDR03",  "vector<float>",   &_muCHIsoDR03);
  _tree->Branch("muNHIsoDR03",  "vector<float>",   &_muNHIsoDR03);
  _tree->Branch("muPhIsoDR03",  "vector<float>",   &_muPhIsoDR03);
  _tree->Branch("muMuIsoDR03",  "vector<float>",   &_muMuIsoDR03);
  _tree->Branch("muElIsoDR03",  "vector<float>",   &_muElIsoDR03);
  _tree->Branch("muCHIsoDR04",  "vector<float>",   &_muCHIsoDR04);
  _tree->Branch("muNHIsoDR04",  "vector<float>",   &_muNHIsoDR04);
  _tree->Branch("muPhIsoDR04",  "vector<float>",   &_muPhIsoDR04);
  _tree->Branch("muMuIsoDR04",  "vector<float>",   &_muMuIsoDR04);
  _tree->Branch("muElIsoDR04",  "vector<float>",   &_muElIsoDR04);
  _tree->Branch("muTrackPt",    "vector<float>",   &_muTrackPt);
  _tree->Branch("muTrackP",     "vector<float>",   &_muTrackP);
  _tree->Branch("muSigmaPOverPTrack",       "vector<float>",    &_muSigmaPOverPTrack);
  _tree->Branch("muNHitsVTX",   "vector<int>",     &_muNHitsVTX);
  _tree->Branch("muNHitsFTD",   "vector<int>",     &_muNHitsFTD);
  _tree->Branch("muNHitsSIT",   "vector<int>",     &_muNHitsSIT);
  _tree->Branch("muNHitsTPC",   "vector<int>",     &_muNHitsTPC);
  _tree->Branch("muNHitsSET",   "vector<int>",     &_muNHitsSET);
  _tree->Branch("muNHitsETD",   "vector<int>",     &_muNHitsETD);
  _tree->Branch("muNHitsVTXVeto",   "vector<int>",     &_muNHitsVTXVeto);
  _tree->Branch("muNHitsFTDVeto",   "vector<int>",     &_muNHitsFTDVeto);
  _tree->Branch("muNHitsSITVeto",   "vector<int>",     &_muNHitsSITVeto);
  _tree->Branch("muNHitsTPCVeto",   "vector<int>",     &_muNHitsTPCVeto);
  _tree->Branch("muNHitsSETVeto",   "vector<int>",     &_muNHitsSETVeto);
  _tree->Branch("muNHitsETDVeto",   "vector<int>",     &_muNHitsETDVeto);
  _tree->Branch("mu1NHit",       "vector<int>",     &_mu1NHits);
  _tree->Branch("mu1NHitsVeto",  "vector<int>",     &_mu1NHitsVeto);
  _tree->Branch("muCHIsoCosTheta995",  "vector<float>",   &_muCHIsoCosTheta995);
  _tree->Branch("muNHIsoCosTheta995",  "vector<float>",   &_muNHIsoCosTheta995);
  _tree->Branch("muPhIsoCosTheta995",  "vector<float>",   &_muPhIsoCosTheta995);
  _tree->Branch("muMuIsoCosTheta995",  "vector<float>",   &_muMuIsoCosTheta995);
  _tree->Branch("muElIsoCosTheta995",  "vector<float>",   &_muElIsoCosTheta995);
  _tree->Branch("muID",  "vector<int>",   &_muID);
  _tree->Branch("muOmega","vector<float>",&_muOmega);
  
  _tree->Branch("phPx",  "vector<float>",   &_phPx); 
  _tree->Branch("phPy",  "vector<float>",   &_phPy); 
  _tree->Branch("phPz",  "vector<float>",   &_phPz); 
  _tree->Branch("phE",   "vector<float>",   &_phE); 
  _tree->Branch("phPhi", "vector<float>",   &_phPhi); 
  _tree->Branch("phTheta","vector<float>",  &_phTheta); 
  _tree->Branch("phMass","vector<float>",   &_phMass);
  _tree->Branch("phEoverTot",  "vector<float>",   &_phEoverTot);
  _tree->Branch("phClusterE",  "vector<float>",   &_phClusterE);
  _tree->Branch("phECAL",      "vector<float>",   &_phECAL);
  _tree->Branch("phHCAL",      "vector<float>",   &_phHCAL);
  _tree->Branch("phf1ClusterE",  "vector<float>",   &_phf1ClusterE);
  _tree->Branch("phf1ECAL",      "vector<float>",   &_phf1ECAL);
  _tree->Branch("phf1HCAL",      "vector<float>",   &_phf1HCAL);
  _tree->Branch("phT1Phi",      "vector<float>",  &_phT1Phi);
  _tree->Branch("phT1Chi2_NDOF","vector<float>",  &_phT1Chi2_NDOF);
  _tree->Branch("phT2Phi",      "vector<float>",  &_phT2Phi);
  _tree->Branch("phT2Chi2_NDOF","vector<float>",  &_phT2Chi2_NDOF);
  _tree->Branch("phCHIsoDR03",  "vector<float>",   &_phCHIsoDR03);
  _tree->Branch("phNHIsoDR03",  "vector<float>",   &_phNHIsoDR03);
  _tree->Branch("phPhIsoDR03",  "vector<float>",   &_phPhIsoDR03);
  _tree->Branch("phMuIsoDR03",  "vector<float>",   &_phMuIsoDR03);
  _tree->Branch("phElIsoDR03",  "vector<float>",   &_phElIsoDR03);
  _tree->Branch("phCHIsoDR04",  "vector<float>",   &_phCHIsoDR04);
  _tree->Branch("phNHIsoDR04",  "vector<float>",   &_phNHIsoDR04);
  _tree->Branch("phPhIsoDR04",  "vector<float>",   &_phPhIsoDR04);
  _tree->Branch("phMuIsoDR04",  "vector<float>",   &_phMuIsoDR04);
  _tree->Branch("phElIsoDR04",  "vector<float>",   &_phElIsoDR04);
  _tree->Branch("phT1NHitsFitted",  "vector<int>",     &_phT1NHitsFitted);
  _tree->Branch("phT1NHitsVeto",    "vector<int>",     &_phT1NHitsVeto);
  _tree->Branch("phT2NHitsFitted",  "vector<int>",     &_phT2NHitsFitted);
  _tree->Branch("ph21NHitsVeto",    "vector<int>",     &_phT2NHitsVeto);
  _tree->Branch("phCHIsoCosTheta995",  "vector<float>",   &_phCHIsoCosTheta995);
  _tree->Branch("phNHIsoCosTheta995",  "vector<float>",   &_phNHIsoCosTheta995);
  _tree->Branch("phPhIsoCosTheta995",  "vector<float>",   &_phPhIsoCosTheta995);
  _tree->Branch("phElIsoCosTheta995",  "vector<float>",   &_phElIsoCosTheta995);
  _tree->Branch("phMuIsoCosTheta995",  "vector<float>",   &_phMuIsoCosTheta995);
  //_tree->Branch("phLongShowerProfile",      "vector<float>",    &_phLongShowerProfile);
  //_tree->Branch("phTransShowerProfile",     "vector<float>",    &_phTransShowerProfile);
  //_tree->Branch("phTrackBasedTransProfile", "vector<float>",    &_phTrackBasedTransProfile);
  _tree->Branch("phEnergyEE",               "vector<float>",    &_phEnergyEE);
  _tree->Branch("phEnergyEB",               "vector<float>",    &_phEnergyEB);
  _tree->Branch("phEnergyEElse",            "vector<float>",    &_phEnergyEElse);
  _tree->Branch("phEnergyHE",               "vector<float>",    &_phEnergyHE);
  _tree->Branch("phEnergyHB",               "vector<float>",    &_phEnergyHB);
  _tree->Branch("phEnergyHElse",            "vector<float>",    &_phEnergyHElse);
  _tree->Branch("phEnergyElseB",            "vector<float>",    &_phEnergyElseB);
  _tree->Branch("phEnergyElseE",            "vector<float>",    &_phEnergyElseE);
  _tree->Branch("phEnergyElseElse",         "vector<float>",    &_phEnergyElseElse);
  _tree->Branch("phEnergyEM",               "vector<float>",    &_phEnergyEM);
  _tree->Branch("phEnergyHAD",              "vector<float>",    &_phEnergyHAD);
  _tree->Branch("phEnergyMuon",             "vector<float>",    &_phEnergyMuon);
  //_tree->Branch("phPeakEnergy",             "vector<float>",    &_phPeakEnergy);
  //_tree->Branch("phNRadiationLength",       "vector<float>",    &_phNRadiationLength);
  _tree->Branch("phMaxEInECALLayer",        "vector<float>",    &_phMaxEInECALLayer);
  _tree->Branch("phMaxEInHCALLayer",        "vector<float>",    &_phMaxEInHCALLayer);
  _tree->Branch("phMinTrackClustDiff",      "vector<float>",    &_phMinTrackClustDiff);
  _tree->Branch("phSigmaPOverPClosestTrack","vector<float>",    &_phSigmaPOverPClosestTrack);
  _tree->Branch("phClosestTrackPt",         "vector<float>",    &_phClosestTrackPt);
  _tree->Branch("phClosestTrackP",          "vector<float>",    &_phClosestTrackP);
  _tree->Branch("phMinSelTrackClustDiff",      "vector<float>",    &_phMinSelTrackClustDiff);
  _tree->Branch("phSigmaPOverPClosestSelTrack","vector<float>",    &_phSigmaPOverPClosestSelTrack);
  _tree->Branch("phClosestSelTrackPt",         "vector<float>",    &_phClosestSelTrackPt);
  _tree->Branch("phClosestSelTrackP",          "vector<float>",    &_phClosestSelTrackP);
  _tree->Branch("phFirstLayerECAL",       "vector<int>",    &_phFirstLayerECAL);
  _tree->Branch("phLastLayerECAL",        "vector<int>",    &_phLastLayerECAL);
  _tree->Branch("phLayerMaxHitsECAL",     "vector<int>",    &_phLayerMaxHitsECAL);
  _tree->Branch("phLayerMaxEECAL",        "vector<int>",    &_phLayerMaxEECAL);
  _tree->Branch("phHitsECAL",             "vector<int>",    &_phHitsECAL);
  _tree->Branch("phHitsHCAL",             "vector<int>",    &_phHitsHCAL);
  _tree->Branch("phMaxHitsPerLayerECAL",  "vector<int>",    &_phMaxHitsPerLayerECAL);
  _tree->Branch("phFirstLayerHCAL",       "vector<int>",    &_phFirstLayerHCAL);
  _tree->Branch("phLastLayerHCAL",         "vector<int>",   &_phLastLayerHCAL);
  _tree->Branch("phMaxHitsPerLayerHCAL",  "vector<int>",    &_phMaxHitsPerLayerHCAL);
  _tree->Branch("phClosestTrackD0",         "vector<float>",      &_phClosestTrackD0);
  _tree->Branch("phClosestTrackZ0",         "vector<float>",      &_phClosestTrackZ0);
  _tree->Branch("phClosestTrackPhi",        "vector<float>",      &_phClosestTrackPhi);
  _tree->Branch("phClosestTrackChi2_NDOF",  "vector<float>",      &_phClosestTrackChi2_NDOF);
  _tree->Branch("phClosestTrackOmega",      "vector<float>",      &_phClosestTrackOmega);
  _tree->Branch("phClosestTrackNHitsVTX",   "vector<int>",      &_phClosestTrackNHitsVTX);
  _tree->Branch("phClosestTrackNHitsFTD",   "vector<int>",      &_phClosestTrackNHitsFTD);
  _tree->Branch("phClosestTrackNHitsSIT",   "vector<int>",      &_phClosestTrackNHitsSIT);
  _tree->Branch("phClosestTrackNHitsTPC",   "vector<int>",      &_phClosestTrackNHitsTPC);
  _tree->Branch("phClosestTrackNHitsSET",   "vector<int>",      &_phClosestTrackNHitsSET);
  _tree->Branch("phClosestTrackNHitsETD",   "vector<int>",      &_phClosestTrackNHitsETD);
  _tree->Branch("phClosestTrackNHitsVTXVeto",   "vector<int>",      &_phClosestTrackNHitsVTXVeto);
  _tree->Branch("phClosestTrackNHitsFTDVeto",   "vector<int>",      &_phClosestTrackNHitsFTDVeto);
  _tree->Branch("phClosestTrackNHitsSITVeto",   "vector<int>",      &_phClosestTrackNHitsSITVeto);
  _tree->Branch("phClosestTrackNHitsTPCVeto",   "vector<int>",      &_phClosestTrackNHitsTPCVeto);
  _tree->Branch("phClosestTrackNHitsSETVeto",   "vector<int>",      &_phClosestTrackNHitsSETVeto);
  _tree->Branch("phClosestTrackNHitsETDVeto",   "vector<int>",      &_phClosestTrackNHitsETDVeto); 
  _tree->Branch("phClosestTrackClusterOpeningAngle",    "vector<float>",  &_phClosestTrackClusterOpeningAngle);
  _tree->Branch("phClosestSelTrackClusterOpeningAngle", "vector<float>",  &_phClosestSelTrackClusterOpeningAngle);
  
  if(_hasrereco_information){
    _tree->Branch("MExReReco",     &_MExReReco,   "MExReReco/F");
    _tree->Branch("MEyReReco",     &_MEyReReco,   "MEyReReco/F");
    _tree->Branch("MEzReReco",     &_MEzReReco,   "MEzReReco/F");
    _tree->Branch("METReReco",     &_METReReco,   "METReReco/F");//actually MPt
    _tree->Branch("SumEReReco",    &_SumEReReco,  "SumEReReco/F");
    _tree->Branch("SumPtReReco",   &_SumPtReReco, "SumPtReReco/F");
    _tree->Branch("SumChargedHadronEReReco",   &_SumChargedHadronEReReco, "SumChargedHadronEReReco/F");
    _tree->Branch("SumChargedHadronPtReReco",  &_SumChargedHadronPtReReco,"SumChargedHadronPtReReco/F");
    _tree->Branch("SumNeutralHadronEReReco",   &_SumNeutralHadronEReReco, "SumNeutralHadronEReReco/F");
    _tree->Branch("SumNeutralHadronPtReReco",  &_SumNeutralHadronPtReReco,"SumNeutralHadronPtReReco/F");
    _tree->Branch("SumPhotonEReReco",          &_SumPhotonEReReco,   "SumPhotonEReReco/F");
    _tree->Branch("SumPhotonPtReReco",         &_SumPhotonPtReReco,  "SumPhotonPtReReco/F");
    _tree->Branch("SumMuonEReReco",            &_SumMuonEReReco,     "SumMuonEReReco/F");
    _tree->Branch("SumMuonPtReReco",           &_SumMuonPtReReco,    "SumMuonPtReReco/F");
    _tree->Branch("SumElectronEReReco",        &_SumElectronEReReco, "SumElectronEReReco/F");
    _tree->Branch("SumElectronPtReReco",       &_SumElectronPtReReco,"SumElectronPtReReco/F");

    _tree->Branch("jetPxReReco",   "vector<float>",   &_jetPxReReco);
    _tree->Branch("jetPyReReco",   "vector<float>",   &_jetPyReReco);
    _tree->Branch("jetPzReReco",   "vector<float>",   &_jetPzReReco);
    _tree->Branch("jetEReReco",    "vector<float>",   &_jetEReReco);
    _tree->Branch("jetPhiReReco",  "vector<float>",   &_jetPhiReReco);
    _tree->Branch("jetThetaReReco", "vector<float>",  &_jetThetaReReco);
    _tree->Branch("jetMassReReco", "vector<float>",   &_jetMassReReco);
    _tree->Branch("jetNFReReco",   "vector<float>",   &_jetNFReReco);     
    _tree->Branch("jetLFReReco",   "vector<float>",   &_jetLFReReco);   
    _tree->Branch("jetKFReReco",   "vector<float>",   &_jetKFReReco);   
    _tree->Branch("jetCHFReReco",  "vector<float>",   &_jetCHFReReco);  
    _tree->Branch("jetPhFReReco",  "vector<float>",   &_jetPhFReReco);  
    _tree->Branch("jetMuFReReco",  "vector<float>",   &_jetElFReReco);  
    _tree->Branch("jetMuFReReco",  "vector<float>",   &_jetMuFReReco);
    _tree->Branch("jetNMultReReco","vector<int>",     &_jetNMultReReco);     
    _tree->Branch("jetLMultReReco","vector<int>",     &_jetLMultReReco);   
    _tree->Branch("jetKMultReReco","vector<int>",     &_jetKMultReReco);   
    _tree->Branch("jetCHMultReReco","vector<int>",    &_jetCHMultReReco);  
    _tree->Branch("jetPhMultReReco","vector<int>",    &_jetPhMultReReco);  
    _tree->Branch("jetElMultReReco","vector<int>",    &_jetElMultReReco);  
    _tree->Branch("jetMuMultReReco","vector<int>",    &_jetMuMultReReco);
   
    _tree->Branch("elPxReReco",  "vector<float>",   &_elPxReReco); 
    _tree->Branch("elPyReReco",  "vector<float>",   &_elPyReReco); 
    _tree->Branch("elPzReReco",  "vector<float>",   &_elPzReReco); 
    _tree->Branch("elEReReco",   "vector<float>",   &_elEReReco); 
    _tree->Branch("elPhiReReco", "vector<float>",   &_elPhiReReco); 
    _tree->Branch("elThetaReReco","vector<float>",  &_elThetaReReco); 
    _tree->Branch("elMassReReco","vector<float>",   &_elMassReReco);
    _tree->Branch("elClusterEReReco",  "vector<float>",   &_elClusterEReReco);
    _tree->Branch("elEoverTotReReco",  "vector<float>",   &_elEoverTotReReco);
    _tree->Branch("elECALReReco",      "vector<float>",   &_elECALReReco);
    _tree->Branch("elHCALReReco",      "vector<float>",   &_elHCALReReco);
    _tree->Branch("el0D0ReReco",   "vector<float>",  &_el0D0ReReco); 
    _tree->Branch("el0Z0ReReco",   "vector<float>",  &_el0Z0ReReco); 
    _tree->Branch("el0PhiReReco",  "vector<float>",  &_el0PhiReReco); 
    _tree->Branch("el0Chi2_NDOFReReco","vector<float>",&_el0Chi2_NDOFReReco);
    _tree->Branch("el1D0ReReco",   "vector<float>",  &_el1D0ReReco); 
    _tree->Branch("el1Z0ReReco",   "vector<float>",  &_el1Z0ReReco); 
    _tree->Branch("el1PhiReReco",  "vector<float>",  &_el1PhiReReco); 
    _tree->Branch("el1Chi2_NDOFReReco","vector<float>",&_el1Chi2_NDOFReReco);
    _tree->Branch("elCHIsoDR03ReReco",  "vector<float>",   &_elCHIsoDR03ReReco);
    _tree->Branch("elNHIsoDR03ReReco",  "vector<float>",   &_elNHIsoDR03ReReco);
    _tree->Branch("elPhIsoDR03ReReco",  "vector<float>",   &_elPhIsoDR03ReReco);
    _tree->Branch("elMuIsoDR03ReReco",  "vector<float>",   &_elMuIsoDR03ReReco);
    _tree->Branch("elElIsoDR03ReReco",  "vector<float>",   &_elElIsoDR03ReReco);
    _tree->Branch("elCHIsoDR04ReReco",  "vector<float>",   &_elCHIsoDR04ReReco);
    _tree->Branch("elNHIsoDR04ReReco",  "vector<float>",   &_elNHIsoDR04ReReco);
    _tree->Branch("elPhIsoDR04ReReco",  "vector<float>",   &_elPhIsoDR04ReReco);
    _tree->Branch("elMuIsoDR04ReReco",  "vector<float>",   &_elMuIsoDR04ReReco);
    _tree->Branch("elElIsoDR04ReReco",  "vector<float>",   &_elElIsoDR04ReReco);
    _tree->Branch("elNHitsVTXReReco",   "vector<int>",     &_elNHitsVTXReReco);
    _tree->Branch("elNHitsFTDReReco",   "vector<int>",     &_elNHitsFTDReReco);
    _tree->Branch("elNHitsSITReReco",   "vector<int>",     &_elNHitsSITReReco);
    _tree->Branch("elNHitsTPCReReco",   "vector<int>",     &_elNHitsTPCReReco);
    _tree->Branch("elNHitsSETReReco",   "vector<int>",     &_elNHitsSETReReco);
    _tree->Branch("elNHitsETDReReco",   "vector<int>",     &_elNHitsETDReReco);
    _tree->Branch("elNHitsVTXVetoReReco",   "vector<int>",     &_elNHitsVTXVetoReReco);
    _tree->Branch("elNHitsFTDVetoReReco",   "vector<int>",     &_elNHitsFTDVetoReReco);
    _tree->Branch("elNHitsSITVetoReReco",   "vector<int>",     &_elNHitsSITVetoReReco);
    _tree->Branch("elNHitsTPCVetoReReco",   "vector<int>",     &_elNHitsTPCVetoReReco);
    _tree->Branch("elNHitsSETVetoReReco",   "vector<int>",     &_elNHitsSETVetoReReco);
    _tree->Branch("elNHitsETDVetoReReco",   "vector<int>",     &_elNHitsETDVetoReReco);
    _tree->Branch("el1NHitReReco",       "vector<int>",     &_el1NHitsReReco);
    _tree->Branch("el1NHitsVetoReReco",  "vector<int>",     &_el1NHitsVetoReReco);
    _tree->Branch("elCHIsoCosTheta995ReReco",  "vector<float>",   &_elCHIsoCosTheta995ReReco);
    _tree->Branch("elNHIsoCosTheta995ReReco",  "vector<float>",   &_elNHIsoCosTheta995ReReco);
    _tree->Branch("elPhIsoCosTheta995ReReco",  "vector<float>",   &_elPhIsoCosTheta995ReReco);
    _tree->Branch("elMuIsoCosTheta995ReReco",  "vector<float>",   &_elMuIsoCosTheta995ReReco);
    _tree->Branch("elElIsoCosTheta995ReReco",  "vector<float>",   &_elElIsoCosTheta995ReReco);
    //_tree->Branch("elLongShowerProfileReReco",      "vector<float>",    &_elLongShowerProfileReReco);
    //_tree->Branch("elTransShowerProfileReReco",     "vector<float>",    &_elTransShowerProfileReReco);
    //_tree->Branch("elTrackBasedTransProfileReReco", "vector<float>",    &_elTrackBasedTransProfileReReco);
    _tree->Branch("elEnergyEEReReco",               "vector<float>",    &_elEnergyEEReReco);
    _tree->Branch("elEnergyEBReReco",               "vector<float>",    &_elEnergyEBReReco);
    _tree->Branch("elEnergyEElseReReco",            "vector<float>",    &_elEnergyEElseReReco);
    _tree->Branch("elEnergyHEReReco",               "vector<float>",    &_elEnergyHEReReco);
    _tree->Branch("elEnergyHBReReco",               "vector<float>",    &_elEnergyHBReReco);
    _tree->Branch("elEnergyHElseReReco",            "vector<float>",    &_elEnergyHElseReReco);
    _tree->Branch("elEnergyElseBReReco",            "vector<float>",    &_elEnergyElseBReReco);
    _tree->Branch("elEnergyElseEReReco",            "vector<float>",    &_elEnergyElseEReReco);
    _tree->Branch("elEnergyElseElseReReco",         "vector<float>",    &_elEnergyElseElseReReco);
    _tree->Branch("elEnergyEMReReco",               "vector<float>",    &_elEnergyEMReReco);
    _tree->Branch("elEnergyHADReReco",              "vector<float>",    &_elEnergyHADReReco);
    _tree->Branch("elEnergyMuonReReco",             "vector<float>",    &_elEnergyMuonReReco);
    //_tree->Branch("elPeakEnergyReReco",             "vector<float>",    &_elPeakEnergyReReco);
    //_tree->Branch("elNRadiationLengthReReco",       "vector<float>",    &_elNRadiationLengthReReco);
    _tree->Branch("elMaxEInECALLayerReReco",        "vector<float>",    &_elMaxEInECALLayerReReco);
    _tree->Branch("elMaxEInHCALLayerReReco",        "vector<float>",    &_elMaxEInHCALLayerReReco);
    _tree->Branch("elMinTrackClustDiffReReco",      "vector<float>",    &_elMinTrackClustDiffReReco);
    _tree->Branch("elSigmaPOverPTrackReReco",       "vector<float>",    &_elSigmaPOverPTrackReReco);
    _tree->Branch("elTrackPtReReco",                "vector<float>",    &_elTrackPtReReco);
    _tree->Branch("elTrackPReReco",                 "vector<float>",    &_elTrackPReReco);
    _tree->Branch("elFirstLayerECALReReco",       "vector<int>",    &_elFirstLayerECALReReco);
    _tree->Branch("elLastLayerECALReReco",        "vector<int>",    &_elLastLayerECALReReco);
    _tree->Branch("elLayerMaxHitsECALReReco",     "vector<int>",    &_elLayerMaxHitsECALReReco);
    _tree->Branch("elLayerMaxEECALReReco",        "vector<int>",    &_elLayerMaxEECALReReco);
    _tree->Branch("elHitsECALReReco",             "vector<int>",    &_elHitsECALReReco);
    _tree->Branch("elHitsHCALReReco",             "vector<int>",    &_elHitsHCALReReco);
    _tree->Branch("elMaxHitsPerLayerECALReReco",  "vector<int>",    &_elMaxHitsPerLayerECALReReco);
    _tree->Branch("elFirstLayerHCALReReco",       "vector<int>",    &_elFirstLayerHCALReReco);
    _tree->Branch("elLastLayerHCALReReco",        "vector<int>",    &_elLastLayerHCALReReco);
    _tree->Branch("elMaxHitsPerLayerHCALReReco",  "vector<int>",    &_elMaxHitsPerLayerHCALReReco);
    _tree->Branch("elIDReReco",  "vector<int>",   &_elIDReReco);
    _tree->Branch("elOmegaReReco","vector<float>",&_elOmegaReReco);
    _tree->Branch("elTrackClusterOpeningAngleReReco" ,"vector<float>",  &_elTrackClusterOpeningAngleReReco);
    
    _tree->Branch("muPxReReco",  "vector<float>",   &_muPxReReco); 
    _tree->Branch("muPyReReco",  "vector<float>",   &_muPyReReco); 
    _tree->Branch("muPzReReco",  "vector<float>",   &_muPzReReco); 
    _tree->Branch("muEReReco",  "vector<float>",    &_muEReReco); 
    _tree->Branch("muPhiReReco",  "vector<float>",  &_muPhiReReco); 
    _tree->Branch("muThetaReReco","vector<float>",  &_muThetaReReco); 
    _tree->Branch("muMassReReco",  "vector<float>", &_muMassReReco);
    _tree->Branch("muEoverTotReReco",  "vector<float>",   &_muEoverTotReReco);
    _tree->Branch("muClusterEReReco",  "vector<float>",   &_muClusterEReReco);
    _tree->Branch("muECALReReco",      "vector<float>",   &_muECALReReco);
    _tree->Branch("muHCALReReco",      "vector<float>",   &_muHCALReReco);
    _tree->Branch("mu0D0ReReco",   "vector<float>",  &_mu0D0ReReco); 
    _tree->Branch("mu0Z0ReReco",   "vector<float>",  &_mu0Z0ReReco); 
    _tree->Branch("mu0PhiReReco",  "vector<float>",  &_mu0PhiReReco); 
    _tree->Branch("mu0Chi2_NDOFReReco","vector<float>",&_mu0Chi2_NDOFReReco);
    _tree->Branch("muD0ReReco",   "vector<float>",  &_mu1D0ReReco); 
    _tree->Branch("mu1Z0ReReco",   "vector<float>",  &_mu1Z0ReReco); 
    _tree->Branch("mu1PhiReReco",  "vector<float>",  &_mu1PhiReReco); 
    _tree->Branch("mu1Chi2_NDOFReReco","vector<float>",&_mu1Chi2_NDOFReReco);
    _tree->Branch("muCHIsoDR03ReReco",  "vector<float>",   &_muCHIsoDR03ReReco);
    _tree->Branch("muNHIsoDR03ReReco",  "vector<float>",   &_muNHIsoDR03ReReco);
    _tree->Branch("muPhIsoDR03ReReco",  "vector<float>",   &_muPhIsoDR03ReReco);
    _tree->Branch("muMuIsoDR03ReReco",  "vector<float>",   &_muMuIsoDR03ReReco);
    _tree->Branch("muElIsoDR03ReReco",  "vector<float>",   &_muElIsoDR03ReReco);
    _tree->Branch("muCHIsoDR04ReReco",  "vector<float>",   &_muCHIsoDR04ReReco);
    _tree->Branch("muNHIsoDR04ReReco",  "vector<float>",   &_muNHIsoDR04ReReco);
    _tree->Branch("muPhIsoDR04ReReco",  "vector<float>",   &_muPhIsoDR04ReReco);
    _tree->Branch("muMuIsoDR04ReReco",  "vector<float>",   &_muMuIsoDR04ReReco);
    _tree->Branch("muElIsoDR04ReReco",  "vector<float>",   &_muElIsoDR04ReReco);
    _tree->Branch("muNHitsVTXReReco",   "vector<int>",     &_muNHitsVTXReReco);
    _tree->Branch("muNHitsFTDReReco",   "vector<int>",     &_muNHitsFTDReReco);
    _tree->Branch("muNHitsSITReReco",   "vector<int>",     &_muNHitsSITReReco);
    _tree->Branch("muNHitsTPCReReco",   "vector<int>",     &_muNHitsTPCReReco);
    _tree->Branch("muNHitsSETReReco",   "vector<int>",     &_muNHitsSETReReco);
    _tree->Branch("muNHitsETDReReco",   "vector<int>",     &_muNHitsETDReReco);
    _tree->Branch("muNHitsVTXVetoReReco",   "vector<int>",     &_muNHitsVTXVetoReReco);
    _tree->Branch("muNHitsFTDVetoReReco",   "vector<int>",     &_muNHitsFTDVetoReReco);
    _tree->Branch("muNHitsSITVetoReReco",   "vector<int>",     &_muNHitsSITVetoReReco);
    _tree->Branch("muNHitsTPCVetoReReco",   "vector<int>",     &_muNHitsTPCVetoReReco);
    _tree->Branch("muNHitsSETVetoReReco",   "vector<int>",     &_muNHitsSETVetoReReco);
    _tree->Branch("muNHitsETDVetoReReco",   "vector<int>",     &_muNHitsETDVetoReReco);
    _tree->Branch("mu1NHitReReco",       "vector<int>",     &_mu1NHitsReReco);
    _tree->Branch("mu1NHitsVetoReReco",  "vector<int>",     &_mu1NHitsVetoReReco);
    _tree->Branch("muCHIsoCosTheta995ReReco",  "vector<float>",   &_muCHIsoCosTheta995ReReco);
    _tree->Branch("muNHIsoCosTheta995ReReco",  "vector<float>",   &_muNHIsoCosTheta995ReReco);
    _tree->Branch("muPhIsoCosTheta995ReReco",  "vector<float>",   &_muPhIsoCosTheta995ReReco);
    _tree->Branch("muMuIsoCosTheta995ReReco",  "vector<float>",   &_muMuIsoCosTheta995ReReco);
    _tree->Branch("muElIsoCosTheta995ReReco",  "vector<float>",   &_muElIsoCosTheta995ReReco);
    _tree->Branch("muTrackPtReReco",          "vector<float>",    &_muTrackPtReReco);
    _tree->Branch("muTrackPReReco",           "vector<float>",    &_muTrackPReReco);
    _tree->Branch("muSigmaPOverPTrackReReco",       "vector<float>",    &_muSigmaPOverPTrackReReco);
    _tree->Branch("muIDReReco",  "vector<int>",   &_muIDReReco);
    _tree->Branch("muOmegaReReco","vector<float>",&_muOmegaReReco);
    
    _tree->Branch("phPxReReco",  "vector<float>",   &_phPxReReco); 
    _tree->Branch("phPyReReco",  "vector<float>",   &_phPyReReco); 
    _tree->Branch("phPzReReco",  "vector<float>",   &_phPzReReco); 
    _tree->Branch("phEReReco",   "vector<float>",   &_phEReReco); 
    _tree->Branch("phPhiReReco", "vector<float>",   &_phPhiReReco); 
    _tree->Branch("phThetaReReco","vector<float>",  &_phThetaReReco); 
    _tree->Branch("phMassReReco","vector<float>",   &_phMassReReco);
    _tree->Branch("phEoverTotReReco",  "vector<float>",   &_phEoverTotReReco);
    _tree->Branch("phClusterEReReco",  "vector<float>",   &_phClusterEReReco);
    _tree->Branch("phECALReReco",      "vector<float>",   &_phECALReReco);
    _tree->Branch("phHCALReReco",      "vector<float>",   &_phHCALReReco);
    _tree->Branch("phf1ClusterEReReco",  "vector<float>",   &_phf1ClusterEReReco);
    _tree->Branch("phf1ECALReReco",      "vector<float>",   &_phf1ECALReReco);
    _tree->Branch("phf1HCALReReco",      "vector<float>",   &_phf1HCALReReco);
    _tree->Branch("phT1PhiReReco",      "vector<float>",  &_phT1PhiReReco);
    _tree->Branch("phT1Chi2_NDOFReReco","vector<float>",  &_phT1Chi2_NDOFReReco);
    _tree->Branch("phT2PhiReReco",      "vector<float>",  &_phT2PhiReReco);
    _tree->Branch("phT2Chi2_NDOFReReco","vector<float>",  &_phT2Chi2_NDOFReReco);
    _tree->Branch("phCHIsoDR03ReReco",  "vector<float>",   &_phCHIsoDR03ReReco);
    _tree->Branch("phNHIsoDR03ReReco",  "vector<float>",   &_phNHIsoDR03ReReco);
    _tree->Branch("phPhIsoDR03ReReco",  "vector<float>",   &_phPhIsoDR03ReReco);
    _tree->Branch("phMuIsoDR03ReReco",  "vector<float>",   &_phMuIsoDR03ReReco);
    _tree->Branch("phElIsoDR03ReReco",  "vector<float>",   &_phElIsoDR03ReReco);
    _tree->Branch("phCHIsoDR04ReReco",  "vector<float>",   &_phCHIsoDR04ReReco);
    _tree->Branch("phNHIsoDR04ReReco",  "vector<float>",   &_phNHIsoDR04ReReco);
    _tree->Branch("phPhIsoDR04ReReco",  "vector<float>",   &_phPhIsoDR04ReReco);
    _tree->Branch("phMuIsoDR04ReReco",  "vector<float>",   &_phMuIsoDR04ReReco);
    _tree->Branch("phElIsoDR04ReReco",  "vector<float>",   &_phElIsoDR04ReReco);
    _tree->Branch("phT1NHitsFittedReReco",  "vector<int>",     &_phT1NHitsFittedReReco);
    _tree->Branch("phT1NHitsVetoReReco",    "vector<int>",     &_phT1NHitsVetoReReco);
    _tree->Branch("phT2NHitsFittedReReco",  "vector<int>",     &_phT2NHitsFittedReReco);
    _tree->Branch("ph21NHitsVetoReReco",    "vector<int>",     &_phT2NHitsVetoReReco);
    _tree->Branch("phCHIsoCosTheta995ReReco",  "vector<float>",   &_phCHIsoCosTheta995ReReco);
    _tree->Branch("phNHIsoCosTheta995ReReco",  "vector<float>",   &_phNHIsoCosTheta995ReReco);
    _tree->Branch("phPhIsoCosTheta995ReReco",  "vector<float>",   &_phPhIsoCosTheta995ReReco);   
    _tree->Branch("phMuIsoCosTheta995ReReco",  "vector<float>",   &_phMuIsoCosTheta995ReReco);
    _tree->Branch("phElIsoCosTheta995ReReco",  "vector<float>",   &_phElIsoCosTheta995ReReco);  
    //_tree->Branch("phLongShowerProfileReReco",      "vector<float>",    &_phLongShowerProfileReReco);
    //_tree->Branch("phTransShowerProfileReReco",     "vector<float>",    &_phTransShowerProfileReReco);
    //_tree->Branch("phTrackBasedTransProfileReReco", "vector<float>",    &_phTrackBasedTransProfileReReco);
    _tree->Branch("phEnergyEEReReco",               "vector<float>",    &_phEnergyEEReReco);
    _tree->Branch("phEnergyEBReReco",               "vector<float>",    &_phEnergyEBReReco);
    _tree->Branch("phEnergyEElseReReco",            "vector<float>",    &_phEnergyEElseReReco);
    _tree->Branch("phEnergyHEReReco",               "vector<float>",    &_phEnergyHEReReco);
    _tree->Branch("phEnergyHBReReco",               "vector<float>",    &_phEnergyHBReReco);
    _tree->Branch("phEnergyHElseReReco",            "vector<float>",    &_phEnergyHElseReReco);
    _tree->Branch("phEnergyElseBReReco",            "vector<float>",    &_phEnergyElseBReReco);
    _tree->Branch("phEnergyElseEReReco",            "vector<float>",    &_phEnergyElseEReReco);
    _tree->Branch("phEnergyElseElseReReco",         "vector<float>",    &_phEnergyElseElseReReco);
    _tree->Branch("phEnergyEMReReco",               "vector<float>",    &_phEnergyEMReReco);
    _tree->Branch("phEnergyHADReReco",              "vector<float>",    &_phEnergyHADReReco);
    _tree->Branch("phEnergyMuonReReco",             "vector<float>",    &_phEnergyMuonReReco);
    //_tree->Branch("phPeakEnergyReReco",             "vector<float>",    &_phPeakEnergyReReco);
    //_tree->Branch("phNRadiationLengthReReco",       "vector<float>",    &_phNRadiationLengthReReco);
    _tree->Branch("phMaxEInECALLayerReReco",        "vector<float>",    &_phMaxEInECALLayerReReco);
    _tree->Branch("phMaxEInHCALLayerReReco",        "vector<float>",    &_phMaxEInHCALLayerReReco);
    _tree->Branch("phMinTrackClustDiffReReco",      "vector<float>",    &_phMinTrackClustDiffReReco);
    _tree->Branch("phSigmaPOverPClosestTrackReReco","vector<float>",    &_phSigmaPOverPClosestTrackReReco);
    _tree->Branch("phClosestTrackPtReReco",          "vector<float>",    &_phClosestTrackPtReReco);
    _tree->Branch("phClosestTrackPReReco",           "vector<float>",    &_phClosestTrackPReReco);
    _tree->Branch("phMinSelTrackClustDiffReReco",      "vector<float>",    &_phMinSelTrackClustDiffReReco);
    _tree->Branch("phSigmaPOverPClosestSelTrackReReco","vector<float>",    &_phSigmaPOverPClosestSelTrackReReco);
    _tree->Branch("phClosestSelTrackPtReReco",         "vector<float>",    &_phClosestSelTrackPtReReco);
    _tree->Branch("phClosestSelTrackPReReco",          "vector<float>",    &_phClosestSelTrackPReReco);
    _tree->Branch("phFirstLayerECALReReco",       "vector<int>",    &_phFirstLayerECALReReco);
    _tree->Branch("phLastLayerECALReReco",        "vector<int>",    &_phLastLayerECALReReco);
    _tree->Branch("phLayerMaxEECALReReco",        "vector<int>",    &_phLayerMaxEECALReReco);
    _tree->Branch("phLayerMaxHitsECALReReco",     "vector<int>",    &_phLayerMaxHitsECALReReco);
    _tree->Branch("phHitsECALReReco",             "vector<int>",    &_phHitsECALReReco);
    _tree->Branch("phHitsHCALReReco",             "vector<int>",    &_phHitsHCALReReco);
    _tree->Branch("phMaxHitsPerLayerECALReReco",  "vector<int>",    &_phMaxHitsPerLayerECALReReco);
    _tree->Branch("phFirstLayerHCALReReco",       "vector<int>",    &_phFirstLayerHCALReReco);
    _tree->Branch("phLastLayerHCALReReco",        "vector<int>",    &_phLastLayerHCALReReco);
    _tree->Branch("phMaxHitsPerLayerHCALReReco",  "vector<int>",    &_phMaxHitsPerLayerHCALReReco);
    _tree->Branch("phClosestTrackD0ReReco",         "vector<float>",      &_phClosestTrackD0ReReco);
    _tree->Branch("phClosestTrackZ0ReReco",         "vector<float>",      &_phClosestTrackZ0ReReco);
    _tree->Branch("phClosestTrackPhiReReco",        "vector<float>",      &_phClosestTrackPhiReReco);
    _tree->Branch("phClosestTrackChi2_NDOFReReco",  "vector<float>",      &_phClosestTrackChi2_NDOFReReco);
    _tree->Branch("phClosestTrackOmegaReReco",      "vector<float>",      &_phClosestTrackOmegaReReco);
    _tree->Branch("phClosestTrackNHitsVTXReReco",   "vector<int>",      &_phClosestTrackNHitsVTXReReco);
    _tree->Branch("phClosestTrackNHitsFTDReReco",   "vector<int>",      &_phClosestTrackNHitsFTDReReco);
    _tree->Branch("phClosestTrackNHitsSITReReco",   "vector<int>",      &_phClosestTrackNHitsSITReReco);
    _tree->Branch("phClosestTrackNHitsTPCReReco",   "vector<int>",      &_phClosestTrackNHitsTPCReReco);
    _tree->Branch("phClosestTrackNHitsSETReReco",   "vector<int>",      &_phClosestTrackNHitsSETReReco);
    _tree->Branch("phClosestTrackNHitsETDReReco",   "vector<int>",      &_phClosestTrackNHitsETDReReco);
    _tree->Branch("phClosestTrackNHitsVTXVetoReReco",   "vector<int>",      &_phClosestTrackNHitsVTXVetoReReco);
    _tree->Branch("phClosestTrackNHitsFTDVetoReReco",   "vector<int>",      &_phClosestTrackNHitsFTDVetoReReco);
    _tree->Branch("phClosestTrackNHitsSITVetoReReco",   "vector<int>",      &_phClosestTrackNHitsSITVetoReReco);
    _tree->Branch("phClosestTrackNHitsTPCVetoReReco",   "vector<int>",      &_phClosestTrackNHitsTPCVetoReReco);
    _tree->Branch("phClosestTrackNHitsSETVetoReReco",   "vector<int>",      &_phClosestTrackNHitsSETVetoReReco);
    _tree->Branch("phClosestTrackNHitsETDVetoReReco",   "vector<int>",      &_phClosestTrackNHitsETDVetoReReco); 
    _tree->Branch("phClosestTrackClusterOpeningAngleReReco",    "vector<float>",  &_phClosestTrackClusterOpeningAngleReReco);
    _tree->Branch("phClosestSelTrackClusterOpeningAngleReReco", "vector<float>",  &_phClosestSelTrackClusterOpeningAngleReReco);
  }
  if(_hasrerecoNoTrackSigmaPOverP_information){
    _tree->Branch("elPxReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPxReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elPyReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPyReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elPzReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPzReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elEReRecoNoTrkSigmaPOverP",   "vector<float>",   &_elEReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elPhiReRecoNoTrkSigmaPOverP", "vector<float>",   &_elPhiReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elThetaReRecoNoTrkSigmaPOverP","vector<float>",  &_elThetaReRecoNoTrkSigmaPOverP); 
    _tree->Branch("elMassReRecoNoTrkSigmaPOverP","vector<float>",   &_elMassReRecoNoTrkSigmaPOverP);
    _tree->Branch("elClusterEReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elClusterEReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEoverTotReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elEoverTotReRecoNoTrkSigmaPOverP);
    _tree->Branch("elECALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_elECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elHCALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_elHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("el0D0ReRecoNoTrkSigmaPOverP",   "vector<float>",  &_el0D0ReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el0Z0ReRecoNoTrkSigmaPOverP",   "vector<float>",  &_el0Z0ReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el0PhiReRecoNoTrkSigmaPOverP",  "vector<float>",  &_el0PhiReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el0Chi2_NDOFReRecoNoTrkSigmaPOverP","vector<float>",&_el0Chi2_NDOFReRecoNoTrkSigmaPOverP);
    _tree->Branch("el1D0ReRecoNoTrkSigmaPOverP",   "vector<float>",  &_el1D0ReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el1Z0ReRecoNoTrkSigmaPOverP",   "vector<float>",  &_el1Z0ReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el1PhiReRecoNoTrkSigmaPOverP",  "vector<float>",  &_el1PhiReRecoNoTrkSigmaPOverP); 
    _tree->Branch("el1Chi2_NDOFReRecoNoTrkSigmaPOverP","vector<float>",&_el1Chi2_NDOFReRecoNoTrkSigmaPOverP);
    _tree->Branch("elCHIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elCHIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elNHIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elPhIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPhIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMuIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elMuIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elElIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elElIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elCHIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elCHIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elNHIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elPhIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPhIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMuIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elMuIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elElIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elElIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsVTXReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsVTXReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsFTDReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsFTDReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsSITReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsSITReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsTPCReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsTPCReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsSETReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsSETReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsETDReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsETDReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsVTXVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsVTXVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsFTDVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsFTDVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsSITVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsSITVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsTPCVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsTPCVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsSETVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsSETVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHitsETDVetoReRecoNoTrkSigmaPOverP",   "vector<int>",     &_elNHitsETDVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("el1NHitReRecoNoTrkSigmaPOverP",       "vector<int>",     &_el1NHitsReRecoNoTrkSigmaPOverP);
    _tree->Branch("el1NHitsVetoReRecoNoTrkSigmaPOverP",  "vector<int>",     &_el1NHitsVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMuIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elMuIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("elElIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_elElIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    //_tree->Branch("elLongShowerProfileReRecoNoTrkSigmaPOverP",      "vector<float>",    &_elLongShowerProfileReRecoNoTrkSigmaPOverP);
    //_tree->Branch("elTransShowerProfileReRecoNoTrkSigmaPOverP",     "vector<float>",    &_elTransShowerProfileReRecoNoTrkSigmaPOverP);
    //_tree->Branch("elTrackBasedTransProfileReRecoNoTrkSigmaPOverP", "vector<float>",    &_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyEEReRecoNoTrkSigmaPOverP",               "vector<float>",    &_elEnergyEEReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyEBReRecoNoTrkSigmaPOverP",               "vector<float>",    &_elEnergyEBReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyEElseReRecoNoTrkSigmaPOverP",            "vector<float>",    &_elEnergyEElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyHEReRecoNoTrkSigmaPOverP",               "vector<float>",    &_elEnergyHEReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyHBReRecoNoTrkSigmaPOverP",               "vector<float>",    &_elEnergyHBReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyHElseReRecoNoTrkSigmaPOverP",            "vector<float>",    &_elEnergyHElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyElseBReRecoNoTrkSigmaPOverP",            "vector<float>",    &_elEnergyElseBReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyElseEReRecoNoTrkSigmaPOverP",            "vector<float>",    &_elEnergyElseEReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyElseElseReRecoNoTrkSigmaPOverP",         "vector<float>",    &_elEnergyElseElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyEMReRecoNoTrkSigmaPOverP",               "vector<float>",    &_elEnergyEMReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyHADReRecoNoTrkSigmaPOverP",              "vector<float>",    &_elEnergyHADReRecoNoTrkSigmaPOverP);
    _tree->Branch("elEnergyMuonReRecoNoTrkSigmaPOverP",             "vector<float>",    &_elEnergyMuonReRecoNoTrkSigmaPOverP);
    //_tree->Branch("elPeakEnergyReRecoNoTrkSigmaPOverP",             "vector<float>",    &_elPeakEnergyReRecoNoTrkSigmaPOverP);
    //_tree->Branch("elNRadiationLengthReRecoNoTrkSigmaPOverP",       "vector<float>",    &_elNRadiationLengthReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMaxEInECALLayerReRecoNoTrkSigmaPOverP",        "vector<float>",    &_elMaxEInECALLayerReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMaxEInHCALLayerReRecoNoTrkSigmaPOverP",        "vector<float>",    &_elMaxEInHCALLayerReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMinTrackClustDiffReRecoNoTrkSigmaPOverP",      "vector<float>",    &_elMinTrackClustDiffReRecoNoTrkSigmaPOverP);
    _tree->Branch("elSigmaPOverPTrackReRecoNoTrkSigmaPOverP",       "vector<float>",    &_elSigmaPOverPTrackReRecoNoTrkSigmaPOverP);
    _tree->Branch("elTrackPtReRecoNoTrkSigmaPOverP",                "vector<float>",    &_elTrackPtReRecoNoTrkSigmaPOverP);
    _tree->Branch("elTrackPReRecoNoTrkSigmaPOverP",                 "vector<float>",    &_elTrackPReRecoNoTrkSigmaPOverP);
    _tree->Branch("elFirstLayerECALReRecoNoTrkSigmaPOverP",       "vector<int>",    &_elFirstLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elLastLayerECALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_elLastLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elLayerMaxHitsECALReRecoNoTrkSigmaPOverP",     "vector<int>",    &_elLayerMaxHitsECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elLayerMaxEECALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_elLayerMaxEECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elHitsECALReRecoNoTrkSigmaPOverP",             "vector<int>",    &_elHitsECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elHitsHCALReRecoNoTrkSigmaPOverP",             "vector<int>",    &_elHitsHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP",  "vector<int>",    &_elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elFirstLayerHCALReRecoNoTrkSigmaPOverP",       "vector<int>",    &_elFirstLayerHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elLastLayerHCALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_elLastLayerHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP",  "vector<int>",    &_elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("elIDReRecoNoTrkSigmaPOverP",  "vector<int>",   &_elIDReRecoNoTrkSigmaPOverP);
    _tree->Branch("elOmegaReRecoNoTrkSigmaPOverP","vector<float>",&_elOmegaReRecoNoTrkSigmaPOverP);
    
    _tree->Branch("phPxReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPxReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phPyReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPyReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phPzReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPzReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phEReRecoNoTrkSigmaPOverP",   "vector<float>",   &_phEReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phPhiReRecoNoTrkSigmaPOverP", "vector<float>",   &_phPhiReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phThetaReRecoNoTrkSigmaPOverP","vector<float>",  &_phThetaReRecoNoTrkSigmaPOverP); 
    _tree->Branch("phMassReRecoNoTrkSigmaPOverP","vector<float>",   &_phMassReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEoverTotReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phEoverTotReRecoNoTrkSigmaPOverP);
    _tree->Branch("phClusterEReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phClusterEReRecoNoTrkSigmaPOverP);
    _tree->Branch("phECALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_phECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phHCALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_phHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phf1ClusterEReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phf1ClusterEReRecoNoTrkSigmaPOverP);
    _tree->Branch("phf1ECALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_phf1ECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phf1HCALReRecoNoTrkSigmaPOverP",      "vector<float>",   &_phf1HCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT1PhiReRecoNoTrkSigmaPOverP",      "vector<float>",  &_phT1PhiReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT1Chi2_NDOFReRecoNoTrkSigmaPOverP","vector<float>",  &_phT1Chi2_NDOFReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT2PhiReRecoNoTrkSigmaPOverP",      "vector<float>",  &_phT2PhiReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT2Chi2_NDOFReRecoNoTrkSigmaPOverP","vector<float>",  &_phT2Chi2_NDOFReRecoNoTrkSigmaPOverP);
    _tree->Branch("phCHIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phCHIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phNHIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phNHIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phPhIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPhIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMuIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phMuIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phElIsoDR03ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phElIsoDR03ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phCHIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phCHIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phNHIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phNHIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phPhIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPhIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMuIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phMuIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phElIsoDR04ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phElIsoDR04ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT1NHitsFittedReRecoNoTrkSigmaPOverP",  "vector<int>",     &_phT1NHitsFittedReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT1NHitsVetoReRecoNoTrkSigmaPOverP",    "vector<int>",     &_phT1NHitsVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("phT2NHitsFittedReRecoNoTrkSigmaPOverP",  "vector<int>",     &_phT2NHitsFittedReRecoNoTrkSigmaPOverP);
    _tree->Branch("ph21NHitsVetoReRecoNoTrkSigmaPOverP",    "vector<int>",     &_phT2NHitsVetoReRecoNoTrkSigmaPOverP);
    _tree->Branch("phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP);   
    _tree->Branch("phMuIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phMuIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    _tree->Branch("phElIsoCosTheta995ReRecoNoTrkSigmaPOverP",  "vector<float>",   &_phElIsoCosTheta995ReRecoNoTrkSigmaPOverP);
    //_tree->Branch("phLongShowerProfileReRecoNoTrkSigmaPOverP",      "vector<float>",    &_phLongShowerProfileReRecoNoTrkSigmaPOverP);
    //_tree->Branch("phTransShowerProfileReRecoNoTrkSigmaPOverP",     "vector<float>",    &_phTransShowerProfileReRecoNoTrkSigmaPOverP);
    //_tree->Branch("phTrackBasedTransProfileReRecoNoTrkSigmaPOverP", "vector<float>",    &_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyEEReRecoNoTrkSigmaPOverP",               "vector<float>",    &_phEnergyEEReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyEBReRecoNoTrkSigmaPOverP",               "vector<float>",    &_phEnergyEBReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyEElseReRecoNoTrkSigmaPOverP",            "vector<float>",    &_phEnergyEElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyHEReRecoNoTrkSigmaPOverP",               "vector<float>",    &_phEnergyHEReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyHBReRecoNoTrkSigmaPOverP",               "vector<float>",    &_phEnergyHBReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyHElseReRecoNoTrkSigmaPOverP",            "vector<float>",    &_phEnergyHElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyElseBReRecoNoTrkSigmaPOverP",            "vector<float>",    &_phEnergyElseBReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyElseEReRecoNoTrkSigmaPOverP",            "vector<float>",    &_phEnergyElseEReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyElseElseReRecoNoTrkSigmaPOverP",         "vector<float>",    &_phEnergyElseElseReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyEMReRecoNoTrkSigmaPOverP",               "vector<float>",    &_phEnergyEMReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyHADReRecoNoTrkSigmaPOverP",              "vector<float>",    &_phEnergyHADReRecoNoTrkSigmaPOverP);
    _tree->Branch("phEnergyMuonReRecoNoTrkSigmaPOverP",             "vector<float>",    &_phEnergyMuonReRecoNoTrkSigmaPOverP);
    //_tree->Branch("phPeakEnergyReRecoNoTrkSigmaPOverP",             "vector<float>",    &_phPeakEnergyReRecoNoTrkSigmaPOverP);
    //_tree->Branch("phNRadiationLengthReRecoNoTrkSigmaPOverP",       "vector<float>",    &_phNRadiationLengthReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMaxEInECALLayerReRecoNoTrkSigmaPOverP",        "vector<float>",    &_phMaxEInECALLayerReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMaxEInHCALLayerReRecoNoTrkSigmaPOverP",        "vector<float>",    &_phMaxEInHCALLayerReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMinTrackClustDiffReRecoNoTrkSigmaPOverP",      "vector<float>",    &_phMinTrackClustDiffReRecoNoTrkSigmaPOverP);
    _tree->Branch("phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP","vector<float>",    &_phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP);
    _tree->Branch("phClosestTrackPtReRecoNoTrkSigmaPOverP",          "vector<float>",    &_phClosestTrackPtReRecoNoTrkSigmaPOverP);
    _tree->Branch("phClosestTrackPReRecoNoTrkSigmaPOverP",           "vector<float>",    &_phClosestTrackPReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP",      "vector<float>",    &_phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP);
    _tree->Branch("phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP","vector<float>",    &_phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP);
    _tree->Branch("phClosestSelTrackPtReRecoNoTrkSigmaPOverP",         "vector<float>",    &_phClosestSelTrackPtReRecoNoTrkSigmaPOverP);
    _tree->Branch("phClosestSelTrackPReRecoNoTrkSigmaPOverP",          "vector<float>",    &_phClosestSelTrackPReRecoNoTrkSigmaPOverP);
    _tree->Branch("phFirstLayerECALReRecoNoTrkSigmaPOverP",       "vector<int>",    &_phFirstLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phLastLayerECALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_phLastLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phLayerMaxEECALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_phLayerMaxEECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phLayerMaxHitsECALReRecoNoTrkSigmaPOverP",     "vector<int>",    &_phLayerMaxHitsECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phHitsECALReRecoNoTrkSigmaPOverP",             "vector<int>",    &_phHitsECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phHitsHCALReRecoNoTrkSigmaPOverP",             "vector<int>",    &_phHitsHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP",  "vector<int>",    &_phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phFirstLayerHCALReRecoNoTrkSigmaPOverP",       "vector<int>",    &_phFirstLayerHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phLastLayerHCALReRecoNoTrkSigmaPOverP",        "vector<int>",    &_phLastLayerHCALReRecoNoTrkSigmaPOverP);
    _tree->Branch("phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP",  "vector<int>",    &_phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP);
  }

}


// ---- method that re-initializes the tree branches --------------
void MyDemoAnalyzer::clearTree(){

  _numRun=-1;
  _numEvt=-1;
  _weight=-1;

  _mcpartStatus          ->clear();
  _mcpartPx              ->clear();
  _mcpartPy              ->clear();
  _mcpartPz              ->clear();
  _mcpartE               ->clear();
  _mcpartPhi             ->clear();
  _mcpartTheta           ->clear();
  _mcpartMass            ->clear();
  _mcpartPDGID           ->clear();
  //check for status to find out what is going on
  _mcpartParent1PDGID    ->clear();
  _mcpartParent2PDGID    ->clear();
  _mcpartDaughter1PDGID  ->clear();
  _mcpartDaughter2PDGID  ->clear();
  _mcpartDaughter3PDGID  ->clear();

  //fill only status 1 here
  _elGENPx              ->clear();
  _elGENPy              ->clear();
  _elGENPz              ->clear();
  _elGENE               ->clear();
  _elGENPhi             ->clear();
  _elGENTheta           ->clear();
  _elGENIsoDR03         ->clear();
  _elGENIsoDR04         ->clear();
  _elGENIsoCosTheta995  ->clear();
  _elGENPDGID           ->clear();
  _elGENParent1PDGID    ->clear();
  _elGENParent2PDGID    ->clear();
  _elGENGrandParent1PDGID->clear();

  _muGENPx              ->clear();
  _muGENPy              ->clear();
  _muGENPz              ->clear();
  _muGENE               ->clear();
  _muGENPhi             ->clear();
  _muGENTheta           ->clear();
  _muGENIsoDR03         ->clear();
  _muGENIsoDR04         ->clear();
  _muGENIsoCosTheta995  ->clear();
  _muGENPDGID           ->clear();
  _muGENParent1PDGID    ->clear();
  _muGENParent2PDGID    ->clear();
  _muGENGrandParent1PDGID->clear();

  _phGENPx              ->clear();
  _phGENPy              ->clear();
  _phGENPz              ->clear();
  _phGENE               ->clear();
  _phGENPhi             ->clear();
  _phGENTheta           ->clear();
  _phGENIsoDR03         ->clear();
  _phGENIsoDR04         ->clear();
  _phGENIsoCosTheta995  ->clear();
  _phGENParent1PDGID    ->clear();
  _phGENParent2PDGID    ->clear();
  _phGENGrandParent1PDGID  ->clear();     

  _nuGENPx              ->clear();
  _nuGENPy              ->clear();
  _nuGENPz              ->clear();
  _nuGENE               ->clear();
  _nuGENPhi             ->clear();
  _nuGENTheta           ->clear();
  _nuGENPDGID           ->clear();
  _nuGENParent1PDGID    ->clear();
  _nuGENParent2PDGID    ->clear();
  _nuGENGrandParent1PDGID  ->clear();     

  _MExGEN=0;
  _MEyGEN=0;
  _MEzGEN=0;
  _METGEN=0;//actually missing pt due to neutrinos/LSPs
  _MEGEN=0;
  _SumEGEN=0;
  _SumPtGEN=0;
  _SumChargedHadronEGEN=0;
  _SumPhotonEGEN=0;
  _SumNeutralHadronEGEN=0;
  _SumChargedHadronPtGEN=0;
  _SumPhotonPtGEN=0;
  _SumNeutralHadronPtGEN=0;
  _SumMuonEGEN=0;
  _SumMuonPtGEN=0;
  _SumElectronEGEN=0;
  _SumElectronPtGEN=0;

  _MEx=0;
  _MEy=0;
  _MEz=0;
  _MET=0;
  _SumE=0;
  _SumPt=0;
  _SumChargedHadronE=0;
  _SumPhotonE=0;
  _SumNeutralHadronE=0;
  _SumChargedHadronPt=0;
  _SumPhotonPt=0;
  _SumNeutralHadronPt=0;
  _SumMuonE=0;
  _SumMuonPt=0;
  _SumElectronE=0;
  _SumElectronPt=0;
  
  _jetPy        ->clear(); 
  _jetPz        ->clear(); 
  _jetE         ->clear(); 
  _jetPhi       ->clear(); 
  _jetTheta     ->clear(); 
  _jetMass      ->clear(); 
  _jetNF        ->clear();
  _jetLF        ->clear();
  _jetKF        ->clear();
  _jetCHF       ->clear();
  _jetPhF       ->clear();
  _jetElF       ->clear();
  _jetMuF       ->clear();
  _jetNMult     ->clear();
  _jetLMult     ->clear();
  _jetKMult     ->clear();
  _jetCHMult    ->clear();
  _jetPhMult    ->clear();
  _jetElMult    ->clear();
  _jetMuMult    ->clear();
 
  _elPx         ->clear(); 
  _elPy         ->clear(); 
  _elPz         ->clear(); 
  _elE          ->clear(); 
  _elPhi        ->clear(); 
  _elTheta      ->clear(); 
  _elMass       ->clear();
  _elEoverTot   ->clear();
  _elECAL       ->clear();
  _elHCAL       ->clear();
  _elClusterE   ->clear();
  _el0D0        ->clear();
  _el0Z0        ->clear();
  _el0Phi       ->clear();
  _el0Chi2_NDOF ->clear();
  _el1D0        ->clear();
  _el1Z0        ->clear();
  _el1Phi       ->clear();
  _el1Chi2_NDOF ->clear();
  _elCHIsoDR03  ->clear();
  _elNHIsoDR03  ->clear();
  _elPhIsoDR03  ->clear();
  _elCHIsoDR04  ->clear();
  _elNHIsoDR04  ->clear();
  _elPhIsoDR04  ->clear();
  _elCHIsoCosTheta995  ->clear();
  _elNHIsoCosTheta995  ->clear();
  _elPhIsoCosTheta995  ->clear();
  //_elLongShowerProfile       ->clear();
  //_elTransShowerProfile      ->clear();
  //_elTrackBasedTransProfile  ->clear();
  _elEnergyEE                ->clear();
  _elEnergyEB                ->clear();
  _elEnergyEElse             ->clear();
  _elEnergyHE                ->clear();
  _elEnergyHB                ->clear();
  _elEnergyHElse             ->clear();
  _elEnergyElseB             ->clear();
  _elEnergyElseE             ->clear();
  _elEnergyElseElse          ->clear();
  _elEnergyEM                ->clear();
  _elEnergyHAD               ->clear();
  _elEnergyMuon              ->clear();
  //_elPeakEnergy              ->clear();
  //_elNRadiationLength        ->clear();
  _elFirstLayerECAL          ->clear();
  _elLastLayerECAL           ->clear();
  _elLayerMaxHitsECAL        ->clear();
  _elHitsECAL                ->clear();
  _elHitsHCAL                ->clear();
  _elMaxHitsPerLayerECAL     ->clear();
  _elFirstLayerHCAL          ->clear();
  _elLastLayerHCAL           ->clear();
  _elLayerMaxEECAL           ->clear();
  _elMaxHitsPerLayerHCAL     ->clear();
  _elMaxEInECALLayer         ->clear();
  _elMaxEInHCALLayer         ->clear();
  _elMinTrackClustDiff       ->clear();
  _elSigmaPOverPTrack        ->clear();
  _elTrackP                  ->clear(); 
  _elTrackPt                 ->clear();
  _elNHitsVTX   ->clear();
  _elNHitsFTD   ->clear();
  _elNHitsSIT   ->clear();
  _elNHitsTPC   ->clear();
  _elNHitsSET   ->clear();
  _elNHitsETD   ->clear();
  _elNHitsVTXVeto   ->clear();
  _elNHitsFTDVeto   ->clear();
  _elNHitsSITVeto   ->clear();
  _elNHitsTPCVeto   ->clear();
  _elNHitsSETVeto   ->clear();
  _elNHitsETDVeto   ->clear();
  _el1NHitsVeto ->clear();
  _el1NHits     ->clear();
  _elID         ->clear();
  _elOmega      ->clear();
  _elTrackClusterOpeningAngle  ->clear();

  _muPx         ->clear(); 
  _muPy         ->clear(); 
  _muPz         ->clear(); 
  _muE          ->clear(); 
  _muPhi        ->clear(); 
  _muTheta      ->clear(); 
  _muMass       ->clear();
  _muEoverTot   ->clear();
  _muECAL       ->clear();
  _muHCAL       ->clear();
  _muClusterE   ->clear();
  _mu0D0         ->clear();
  _mu0Z0         ->clear();
  _mu0Phi        ->clear();
  _mu0Chi2_NDOF  ->clear();
  _mu1D0         ->clear();
  _mu1Z0         ->clear();
  _mu1Phi        ->clear();
  _mu1Chi2_NDOF  ->clear();
  _muCHIsoDR03  ->clear();
  _muNHIsoDR03  ->clear();
  _muPhIsoDR03  ->clear();
  _muCHIsoDR04  ->clear();
  _muNHIsoDR04  ->clear();
  _muPhIsoDR04  ->clear();
  _muCHIsoCosTheta995  ->clear();
  _muNHIsoCosTheta995  ->clear();
  _muPhIsoCosTheta995  ->clear();
  _muNHitsVTX   ->clear();
  _muNHitsFTD   ->clear();
  _muNHitsSIT   ->clear();
  _muNHitsTPC   ->clear();
  _muNHitsSET   ->clear();
  _muNHitsETD   ->clear();
  _muNHitsVTXVeto   ->clear();
  _muNHitsFTDVeto   ->clear();
  _muNHitsSITVeto   ->clear();
  _muNHitsTPCVeto   ->clear();
  _muNHitsSETVeto   ->clear();
  _muNHitsETDVeto   ->clear();
  _muTrackP         ->clear(); 
  _muTrackPt        ->clear();
  _muSigmaPOverPTrack     ->clear();
  _mu1NHits     ->clear();
  _mu1NHitsVeto ->clear();
  _muID         ->clear();
  _muOmega      ->clear();
  
  _phPx         ->clear(); 
  _phPy         ->clear(); 
  _phPz         ->clear(); 
  _phE          ->clear(); 
  _phPhi        ->clear(); 
  _phTheta       ->clear(); 
  _phMass       ->clear();
  _phEoverTot   ->clear();
  _phECAL       ->clear();
  _phHCAL       ->clear();
  _phClusterE   ->clear();
  _phf1ECAL       ->clear();
  _phf1HCAL       ->clear();
  _phf1ClusterE   ->clear();
  _phCHIsoDR03  ->clear();
  _phNHIsoDR03  ->clear();
  _phPhIsoDR03  ->clear();
  _phCHIsoDR04  ->clear();
  _phNHIsoDR04  ->clear();
  _phPhIsoDR04  ->clear();
  _phCHIsoCosTheta995  ->clear();
  _phNHIsoCosTheta995  ->clear();
  _phPhIsoCosTheta995  ->clear();
  //_phLongShowerProfile       ->clear();
  //_phTransShowerProfile      ->clear();
  //_phTrackBasedTransProfile  ->clear();
  _phEnergyEE                ->clear();
  _phEnergyEB                ->clear();
  _phEnergyEElse             ->clear();
  _phEnergyHE                ->clear();
  _phEnergyHB                ->clear();
  _phEnergyHElse             ->clear();
  _phEnergyElseB             ->clear();
  _phEnergyElseE             ->clear();
  _phEnergyElseElse          ->clear();
  _phEnergyEM                ->clear();
  _phEnergyHAD               ->clear();
  _phEnergyMuon              ->clear();
  //_phPeakEnergy              ->clear();
  //_phNRadiationLength        ->clear();
  _phFirstLayerECAL          ->clear();
  _phLastLayerECAL           ->clear();
  _phLayerMaxHitsECAL        ->clear();
  _phHitsECAL                ->clear();
  _phHitsHCAL                ->clear();
  _phMaxHitsPerLayerECAL     ->clear();
  _phFirstLayerHCAL          ->clear();
  _phLastLayerHCAL           ->clear();
  _phLayerMaxEECAL           ->clear();
  _phMaxHitsPerLayerHCAL     ->clear();
  _phMaxEInECALLayer         ->clear();
  _phMaxEInHCALLayer         ->clear();
  _phMinTrackClustDiff       ->clear();
  _phSigmaPOverPClosestTrack ->clear();
  _phClosestTrackP              ->clear(); 
  _phClosestTrackPt             ->clear();
  _phClosestTrackD0        ->clear();
  _phClosestTrackZ0        ->clear();
  _phClosestTrackPhi       ->clear();
  _phClosestTrackChi2_NDOF ->clear();    
  _phClosestTrackOmega ->clear();  
  _phClosestTrackNHitsVTX   ->clear();
  _phClosestTrackNHitsFTD   ->clear();
  _phClosestTrackNHitsSIT   ->clear();
  _phClosestTrackNHitsTPC   ->clear();
  _phClosestTrackNHitsSET   ->clear();
  _phClosestTrackNHitsETD   ->clear();
  _phClosestTrackNHitsVTXVeto   ->clear();
  _phClosestTrackNHitsFTDVeto   ->clear();
  _phClosestTrackNHitsSITVeto   ->clear();
  _phClosestTrackNHitsTPCVeto   ->clear();
  _phClosestTrackNHitsSETVeto   ->clear();
  _phClosestTrackNHitsETDVeto   ->clear();
  _phMinSelTrackClustDiff       ->clear();
  _phSigmaPOverPClosestSelTrack ->clear();
  _phClosestSelTrackP           ->clear(); 
  _phClosestSelTrackPt          ->clear();
  _phT1Chi2_NDOF  ->clear();
  _phT1Phi        ->clear();
  _phT1NHitsFitted->clear();
  _phT1NHitsVeto  ->clear();
  _phT2Chi2_NDOF  ->clear();
  _phT2Phi        ->clear();
  _phT2NHitsFitted->clear();
  _phT2NHitsVeto  ->clear();
  _phClosestTrackClusterOpeningAngle  ->clear();
  _phClosestSelTrackClusterOpeningAngle  ->clear();

  if(_hasrereco_information){
    _MExReReco=0;
    _MEyReReco=0;
    _MEzReReco=0;
    _METReReco=0;
    _SumEReReco=0;
    _SumPtReReco=0;
    _SumChargedHadronEReReco=0;
    _SumPhotonEReReco=0;
    _SumNeutralHadronEReReco=0;
    _SumChargedHadronPtReReco=0;
    _SumPhotonPtReReco=0;
    _SumNeutralHadronPtReReco=0;
    _SumMuonEReReco=0;
    _SumMuonPtReReco=0;
    _SumElectronEReReco=0;
    _SumElectronPtReReco=0;
    
    _jetPyReReco        ->clear(); 
    _jetPzReReco        ->clear(); 
    _jetEReReco         ->clear(); 
    _jetPhiReReco       ->clear(); 
    _jetThetaReReco     ->clear(); 
    _jetMassReReco      ->clear(); 
    _jetNFReReco        ->clear();
    _jetLFReReco        ->clear();
    _jetKFReReco        ->clear();
    _jetCHFReReco       ->clear();
    _jetPhFReReco       ->clear();
    _jetElFReReco       ->clear();
    _jetMuFReReco       ->clear();
    _jetNMultReReco     ->clear();
    _jetLMultReReco     ->clear();
    _jetKMultReReco     ->clear();
    _jetCHMultReReco    ->clear();
    _jetPhMultReReco    ->clear();
    _jetElMultReReco    ->clear();
    _jetMuMultReReco    ->clear();
  
    _elPxReReco         ->clear(); 
    _elPyReReco         ->clear(); 
    _elPzReReco         ->clear(); 
    _elEReReco          ->clear(); 
    _elPhiReReco        ->clear(); 
    _elThetaReReco      ->clear(); 
    _elMassReReco       ->clear();
    _elEoverTotReReco   ->clear();
    _elECALReReco       ->clear();
    _elHCALReReco       ->clear();
    _elClusterEReReco   ->clear();
    _el0D0ReReco        ->clear();
    _el0Z0ReReco        ->clear();
    _el0PhiReReco       ->clear();
    _el0Chi2_NDOFReReco ->clear();
    _el1D0ReReco        ->clear();
    _el1Z0ReReco        ->clear();
    _el1PhiReReco       ->clear();
    _el1Chi2_NDOFReReco ->clear();
    _elCHIsoDR03ReReco  ->clear();
    _elNHIsoDR03ReReco  ->clear();
    _elPhIsoDR03ReReco  ->clear();
    _elCHIsoDR04ReReco  ->clear();
    _elNHIsoDR04ReReco  ->clear();
    _elPhIsoDR04ReReco  ->clear();
    _elCHIsoCosTheta995ReReco  ->clear();
    _elNHIsoCosTheta995ReReco  ->clear();
    _elPhIsoCosTheta995ReReco  ->clear();
    //_elLongShowerProfileReReco       ->clear();
    //_elTransShowerProfileReReco      ->clear();
    //_elTrackBasedTransProfileReReco  ->clear();
    _elEnergyEEReReco                ->clear();
    _elEnergyEBReReco                ->clear();
    _elEnergyEElseReReco             ->clear();
    _elEnergyHEReReco                ->clear();
    _elEnergyHBReReco                ->clear();
    _elEnergyHElseReReco             ->clear();
    _elEnergyElseBReReco             ->clear();
    _elEnergyElseEReReco             ->clear();
    _elEnergyElseElseReReco          ->clear();
    _elEnergyEMReReco                ->clear();
    _elEnergyHADReReco               ->clear();
    _elEnergyMuonReReco              ->clear();
    //_elPeakEnergyReReco              ->clear();
    //_elNRadiationLengthReReco        ->clear();
    _elFirstLayerECALReReco          ->clear();
    _elLastLayerECALReReco           ->clear();
    _elLayerMaxEECALReReco           ->clear();
    _elLayerMaxHitsECALReReco        ->clear();
    _elHitsECALReReco                ->clear();
    _elHitsHCALReReco                ->clear();
    _elMaxHitsPerLayerECALReReco     ->clear();
    _elFirstLayerHCALReReco          ->clear();
    _elLastLayerHCALReReco           ->clear();
    _elMaxHitsPerLayerHCALReReco     ->clear();
    _elMaxEInECALLayerReReco         ->clear();
    _elMaxEInHCALLayerReReco         ->clear();
    _elMinTrackClustDiffReReco       ->clear();
    _elSigmaPOverPTrackReReco        ->clear();
    _elTrackPReReco                  ->clear(); 
    _elTrackPtReReco                 ->clear() ;
    _elNHitsVTXReReco   ->clear();
    _elNHitsFTDReReco   ->clear();
    _elNHitsSITReReco   ->clear();
    _elNHitsTPCReReco   ->clear();
    _elNHitsSETReReco   ->clear();
    _elNHitsETDReReco   ->clear();
    _elNHitsVTXVetoReReco   ->clear();
    _elNHitsFTDVetoReReco   ->clear();
    _elNHitsSITVetoReReco   ->clear();
    _elNHitsTPCVetoReReco   ->clear();
    _elNHitsSETVetoReReco   ->clear();
    _elNHitsETDVetoReReco   ->clear();
    _el1NHitsReReco     ->clear();
    _el1NHitsVetoReReco ->clear();
    _elIDReReco         ->clear();
    _elOmegaReReco      ->clear();
    _elTrackClusterOpeningAngleReReco  ->clear();
    
    _muPxReReco         ->clear(); 
    _muPyReReco         ->clear(); 
    _muPzReReco         ->clear(); 
    _muEReReco          ->clear(); 
    _muPhiReReco        ->clear(); 
    _muThetaReReco      ->clear(); 
    _muMassReReco       ->clear();
    _muEoverTotReReco   ->clear();
    _muECALReReco       ->clear();
    _muHCALReReco       ->clear();
    _muClusterEReReco   ->clear();
    _mu0D0ReReco         ->clear();
    _mu0Z0ReReco         ->clear();
    _mu0PhiReReco        ->clear();
    _mu0Chi2_NDOFReReco  ->clear();
    _mu1D0ReReco         ->clear();
    _mu1Z0ReReco         ->clear();
    _mu1PhiReReco        ->clear();
    _mu1Chi2_NDOFReReco  ->clear();
    _muCHIsoDR03ReReco  ->clear();
    _muNHIsoDR03ReReco  ->clear();
    _muPhIsoDR03ReReco  ->clear();
    _muCHIsoDR04ReReco  ->clear();
    _muNHIsoDR04ReReco  ->clear();
    _muPhIsoDR04ReReco  ->clear();
    _muCHIsoCosTheta995ReReco  ->clear();
    _muNHIsoCosTheta995ReReco  ->clear();
    _muPhIsoCosTheta995ReReco  ->clear();
    _muNHitsVTXReReco   ->clear();
    _muNHitsFTDReReco   ->clear();
    _muNHitsSITReReco   ->clear();
    _muNHitsTPCReReco   ->clear();
    _muNHitsSETReReco   ->clear();
    _muNHitsETDReReco   ->clear();
    _muNHitsVTXVetoReReco   ->clear();
    _muNHitsFTDVetoReReco   ->clear();
    _muNHitsSITVetoReReco   ->clear();
    _muNHitsTPCVetoReReco   ->clear();
    _muNHitsSETVetoReReco   ->clear();
    _muNHitsETDVetoReReco   ->clear();
    _mu1NHitsReReco     ->clear();
    _mu1NHitsVetoReReco ->clear();
    _muIDReReco         ->clear();
    _muOmegaReReco      ->clear();
    _muTrackPReReco     ->clear(); 
    _muTrackPtReReco    ->clear();
    _muSigmaPOverPTrackReReco     ->clear();

    _phPxReReco         ->clear(); 
    _phPyReReco         ->clear(); 
    _phPzReReco         ->clear(); 
    _phEReReco          ->clear(); 
    _phPhiReReco        ->clear(); 
    _phThetaReReco       ->clear(); 
    _phMassReReco       ->clear();
    _phEoverTotReReco   ->clear();
    _phECALReReco       ->clear();
    _phHCALReReco       ->clear();
    _phClusterEReReco   ->clear();
    _phf1ECALReReco       ->clear();
    _phf1HCALReReco       ->clear();
    _phf1ClusterEReReco   ->clear();
    _phCHIsoDR03ReReco  ->clear();
    _phNHIsoDR03ReReco  ->clear();
    _phPhIsoDR03ReReco  ->clear();
    _phCHIsoDR04ReReco  ->clear();
    _phNHIsoDR04ReReco  ->clear();
    _phPhIsoDR04ReReco  ->clear();
    _phCHIsoCosTheta995ReReco  ->clear();
    _phNHIsoCosTheta995ReReco  ->clear();
    _phPhIsoCosTheta995ReReco  ->clear();
    //_phLongShowerProfileReReco       ->clear();
    //_phTransShowerProfileReReco      ->clear();
    //_phTrackBasedTransProfileReReco  ->clear();
    _phEnergyEEReReco                ->clear();
    _phEnergyEBReReco                ->clear();
    _phEnergyEElseReReco             ->clear();
    _phEnergyHEReReco                ->clear();
    _phEnergyHBReReco                ->clear();
    _phEnergyHElseReReco             ->clear();
    _phEnergyElseBReReco             ->clear();
    _phEnergyElseEReReco             ->clear();
    _phEnergyElseElseReReco          ->clear();
    _phEnergyEMReReco                ->clear();
    _phEnergyHADReReco               ->clear();
    _phEnergyMuonReReco              ->clear();
    //_phPeakEnergyReReco              ->clear();
    //_phNRadiationLengthReReco        ->clear();
    _phFirstLayerECALReReco          ->clear();
    _phLastLayerECALReReco           ->clear();
    _phLayerMaxEECALReReco           ->clear();
    _phLayerMaxHitsECALReReco        ->clear();
    _phHitsECALReReco                ->clear();
    _phHitsHCALReReco                ->clear();
    _phMaxHitsPerLayerECALReReco     ->clear();
    _phFirstLayerHCALReReco          ->clear();
    _phLastLayerHCALReReco           ->clear();
    _phMaxHitsPerLayerHCALReReco     ->clear();
    _phMaxEInECALLayerReReco         ->clear();
    _phMaxEInHCALLayerReReco         ->clear();
    _phMinTrackClustDiffReReco       ->clear();
    _phSigmaPOverPClosestTrackReReco ->clear();
    _phClosestTrackPReReco           ->clear(); 
    _phClosestTrackPtReReco          ->clear();
    _phClosestTrackD0ReReco        ->clear();
    _phClosestTrackZ0ReReco        ->clear();
    _phClosestTrackPhiReReco       ->clear();
    _phClosestTrackChi2_NDOFReReco ->clear();
    _phClosestTrackOmegaReReco ->clear();  
    _phClosestTrackNHitsVTXReReco   ->clear();
    _phClosestTrackNHitsFTDReReco   ->clear();
    _phClosestTrackNHitsSITReReco   ->clear();
    _phClosestTrackNHitsTPCReReco   ->clear();
    _phClosestTrackNHitsSETReReco   ->clear();
    _phClosestTrackNHitsETDReReco   ->clear();
    _phClosestTrackNHitsVTXVetoReReco   ->clear();
    _phClosestTrackNHitsFTDVetoReReco   ->clear();
    _phClosestTrackNHitsSITVetoReReco   ->clear();
    _phClosestTrackNHitsTPCVetoReReco   ->clear();
    _phClosestTrackNHitsSETVetoReReco   ->clear();
    _phClosestTrackNHitsETDVetoReReco   ->clear();
    _phMinSelTrackClustDiffReReco        ->clear();
    _phSigmaPOverPClosestSelTrackReReco  ->clear();
    _phClosestSelTrackPReReco            ->clear(); 
    _phClosestSelTrackPtReReco           ->clear();
    _phT1Chi2_NDOFReReco  ->clear();
    _phT1PhiReReco        ->clear();
    _phT1NHitsFittedReReco->clear();
    _phT1NHitsVetoReReco  ->clear();
    _phT2Chi2_NDOFReReco  ->clear();
    _phT2PhiReReco        ->clear();
    _phT2NHitsFittedReReco->clear();
    _phT2NHitsVetoReReco  ->clear();
    _phClosestSelTrackClusterOpeningAngleReReco  ->clear();
    _phClosestTrackClusterOpeningAngleReReco  ->clear();
  }

  //now for no track sigma P over P stuff


  if(_hasrerecoNoTrackSigmaPOverP_information){
    _elPxReRecoNoTrkSigmaPOverP         ->clear(); 
    _elPyReRecoNoTrkSigmaPOverP         ->clear(); 
    _elPzReRecoNoTrkSigmaPOverP         ->clear(); 
    _elEReRecoNoTrkSigmaPOverP          ->clear(); 
    _elPhiReRecoNoTrkSigmaPOverP        ->clear(); 
    _elThetaReRecoNoTrkSigmaPOverP      ->clear(); 
    _elMassReRecoNoTrkSigmaPOverP       ->clear();
    _elEoverTotReRecoNoTrkSigmaPOverP   ->clear();
    _elECALReRecoNoTrkSigmaPOverP       ->clear();
    _elHCALReRecoNoTrkSigmaPOverP       ->clear();
    _elClusterEReRecoNoTrkSigmaPOverP   ->clear();
    _el0D0ReRecoNoTrkSigmaPOverP        ->clear();
    _el0Z0ReRecoNoTrkSigmaPOverP        ->clear();
    _el0PhiReRecoNoTrkSigmaPOverP       ->clear();
    _el0Chi2_NDOFReRecoNoTrkSigmaPOverP ->clear();
    _el1D0ReRecoNoTrkSigmaPOverP        ->clear();
    _el1Z0ReRecoNoTrkSigmaPOverP        ->clear();
    _el1PhiReRecoNoTrkSigmaPOverP       ->clear();
    _el1Chi2_NDOFReRecoNoTrkSigmaPOverP ->clear();
    _elCHIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _elNHIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _elPhIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _elCHIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _elNHIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _elPhIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _elCHIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    _elNHIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    _elPhIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    //_elLongShowerProfileReRecoNoTrkSigmaPOverP       ->clear();
    //_elTransShowerProfileReRecoNoTrkSigmaPOverP      ->clear();
    //_elTrackBasedTransProfileReRecoNoTrkSigmaPOverP  ->clear();
    _elEnergyEEReRecoNoTrkSigmaPOverP                ->clear();
    _elEnergyEBReRecoNoTrkSigmaPOverP                ->clear();
    _elEnergyEElseReRecoNoTrkSigmaPOverP             ->clear();
    _elEnergyHEReRecoNoTrkSigmaPOverP                ->clear();
    _elEnergyHBReRecoNoTrkSigmaPOverP                ->clear();
    _elEnergyHElseReRecoNoTrkSigmaPOverP             ->clear();
    _elEnergyElseBReRecoNoTrkSigmaPOverP             ->clear();
    _elEnergyElseEReRecoNoTrkSigmaPOverP             ->clear();
    _elEnergyElseElseReRecoNoTrkSigmaPOverP          ->clear();
    _elEnergyEMReRecoNoTrkSigmaPOverP                ->clear();
    _elEnergyHADReRecoNoTrkSigmaPOverP               ->clear();
    _elEnergyMuonReRecoNoTrkSigmaPOverP              ->clear();
    //_elPeakEnergyReRecoNoTrkSigmaPOverP              ->clear();
    //_elNRadiationLengthReRecoNoTrkSigmaPOverP        ->clear();
    _elFirstLayerECALReRecoNoTrkSigmaPOverP          ->clear();
    _elLastLayerECALReRecoNoTrkSigmaPOverP           ->clear();
    _elLayerMaxEECALReRecoNoTrkSigmaPOverP           ->clear();
    _elLayerMaxHitsECALReRecoNoTrkSigmaPOverP        ->clear();
    _elHitsECALReRecoNoTrkSigmaPOverP                ->clear();
    _elHitsHCALReRecoNoTrkSigmaPOverP                ->clear();
    _elMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP     ->clear();
    _elFirstLayerHCALReRecoNoTrkSigmaPOverP          ->clear();
    _elLastLayerHCALReRecoNoTrkSigmaPOverP           ->clear();
    _elMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP     ->clear();
    _elMaxEInECALLayerReRecoNoTrkSigmaPOverP         ->clear();
    _elMaxEInHCALLayerReRecoNoTrkSigmaPOverP         ->clear();
    _elMinTrackClustDiffReRecoNoTrkSigmaPOverP       ->clear();
    _elSigmaPOverPTrackReRecoNoTrkSigmaPOverP        ->clear();
    _elTrackPReRecoNoTrkSigmaPOverP                  ->clear(); 
    _elTrackPtReRecoNoTrkSigmaPOverP                 ->clear() ;
    _elNHitsVTXReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsFTDReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsSITReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsTPCReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsSETReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsETDReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsVTXVetoReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsFTDVetoReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsSITVetoReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsTPCVetoReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsSETVetoReRecoNoTrkSigmaPOverP   ->clear();
    _elNHitsETDVetoReRecoNoTrkSigmaPOverP   ->clear();
    _el1NHitsReRecoNoTrkSigmaPOverP     ->clear();
    _el1NHitsVetoReRecoNoTrkSigmaPOverP ->clear();
    _elIDReRecoNoTrkSigmaPOverP         ->clear();
    _elOmegaReRecoNoTrkSigmaPOverP      ->clear();

    _phPxReRecoNoTrkSigmaPOverP         ->clear(); 
    _phPyReRecoNoTrkSigmaPOverP         ->clear(); 
    _phPzReRecoNoTrkSigmaPOverP         ->clear(); 
    _phEReRecoNoTrkSigmaPOverP          ->clear(); 
    _phPhiReRecoNoTrkSigmaPOverP        ->clear(); 
    _phThetaReRecoNoTrkSigmaPOverP       ->clear(); 
    _phMassReRecoNoTrkSigmaPOverP       ->clear();
    _phEoverTotReRecoNoTrkSigmaPOverP   ->clear();
    _phECALReRecoNoTrkSigmaPOverP       ->clear();
    _phHCALReRecoNoTrkSigmaPOverP       ->clear();
    _phClusterEReRecoNoTrkSigmaPOverP   ->clear();
    _phf1ECALReRecoNoTrkSigmaPOverP       ->clear();
    _phf1HCALReRecoNoTrkSigmaPOverP       ->clear();
    _phf1ClusterEReRecoNoTrkSigmaPOverP   ->clear();
    _phCHIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _phNHIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _phPhIsoDR03ReRecoNoTrkSigmaPOverP  ->clear();
    _phCHIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _phNHIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _phPhIsoDR04ReRecoNoTrkSigmaPOverP  ->clear();
    _phCHIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    _phNHIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    _phPhIsoCosTheta995ReRecoNoTrkSigmaPOverP  ->clear();
    //_phLongShowerProfileReRecoNoTrkSigmaPOverP       ->clear();
    //_phTransShowerProfileReRecoNoTrkSigmaPOverP      ->clear();
    //_phTrackBasedTransProfileReRecoNoTrkSigmaPOverP  ->clear();
    _phEnergyEEReRecoNoTrkSigmaPOverP                ->clear();
    _phEnergyEBReRecoNoTrkSigmaPOverP                ->clear();
    _phEnergyEElseReRecoNoTrkSigmaPOverP             ->clear();
    _phEnergyHEReRecoNoTrkSigmaPOverP                ->clear();
    _phEnergyHBReRecoNoTrkSigmaPOverP                ->clear();
    _phEnergyHElseReRecoNoTrkSigmaPOverP             ->clear();
    _phEnergyElseBReRecoNoTrkSigmaPOverP             ->clear();
    _phEnergyElseEReRecoNoTrkSigmaPOverP             ->clear();
    _phEnergyElseElseReRecoNoTrkSigmaPOverP          ->clear();
    _phEnergyEMReRecoNoTrkSigmaPOverP                ->clear();
    _phEnergyHADReRecoNoTrkSigmaPOverP               ->clear();
    _phEnergyMuonReRecoNoTrkSigmaPOverP              ->clear();
    //_phPeakEnergyReRecoNoTrkSigmaPOverP              ->clear();
    //_phNRadiationLengthReRecoNoTrkSigmaPOverP        ->clear();
    _phFirstLayerECALReRecoNoTrkSigmaPOverP          ->clear();
    _phLastLayerECALReRecoNoTrkSigmaPOverP           ->clear();
    _phLayerMaxEECALReRecoNoTrkSigmaPOverP           ->clear();
    _phLayerMaxHitsECALReRecoNoTrkSigmaPOverP        ->clear();
    _phHitsECALReRecoNoTrkSigmaPOverP                ->clear();
    _phHitsHCALReRecoNoTrkSigmaPOverP                ->clear();
    _phMaxHitsPerLayerECALReRecoNoTrkSigmaPOverP     ->clear();
    _phFirstLayerHCALReRecoNoTrkSigmaPOverP          ->clear();
    _phLastLayerHCALReRecoNoTrkSigmaPOverP           ->clear();
    _phMaxHitsPerLayerHCALReRecoNoTrkSigmaPOverP     ->clear();
    _phMaxEInECALLayerReRecoNoTrkSigmaPOverP         ->clear();
    _phMaxEInHCALLayerReRecoNoTrkSigmaPOverP         ->clear();
    _phMinTrackClustDiffReRecoNoTrkSigmaPOverP       ->clear();
    _phSigmaPOverPClosestTrackReRecoNoTrkSigmaPOverP ->clear();
    _phClosestTrackPReRecoNoTrkSigmaPOverP           ->clear(); 
    _phClosestTrackPtReRecoNoTrkSigmaPOverP          ->clear();
    _phMinSelTrackClustDiffReRecoNoTrkSigmaPOverP        ->clear();
    _phSigmaPOverPClosestSelTrackReRecoNoTrkSigmaPOverP  ->clear();
    _phClosestSelTrackPReRecoNoTrkSigmaPOverP            ->clear(); 
    _phClosestSelTrackPtReRecoNoTrkSigmaPOverP           ->clear();
    _phT1Chi2_NDOFReRecoNoTrkSigmaPOverP  ->clear();
    _phT1PhiReRecoNoTrkSigmaPOverP        ->clear();
    _phT1NHitsFittedReRecoNoTrkSigmaPOverP->clear();
    _phT1NHitsVetoReRecoNoTrkSigmaPOverP  ->clear();
    _phT2Chi2_NDOFReRecoNoTrkSigmaPOverP  ->clear();
    _phT2PhiReRecoNoTrkSigmaPOverP        ->clear();
    _phT2NHitsFittedReRecoNoTrkSigmaPOverP->clear();
    _phT2NHitsVetoReRecoNoTrkSigmaPOverP  ->clear();
  }

}

//------------------------------------------------------------------------------------------------------------------------------------------



//do the 0.100 GeV cut on pt
void MyDemoAnalyzer::GetTrackStatesOld(const EVENT::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    pandora::Helix *pHelixFit = new pandora::Helix(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_bField);
    trackParameters.m_momentumAtDca = pHelixFit->GetMomentum();    
    const EVENT::TrackerHitVec &trackerHitvec(pTrack->getTrackerHits());
    float zMin(std::numeric_limits<float>::max()), zMax(-std::numeric_limits<float>::max());

    for (int iz = 0, nTrackHits = trackerHitvec.size(); iz < nTrackHits - 1; ++iz)
    {
        const float hitZ(trackerHitvec[iz]->getPosition()[2]);

        if (hitZ > zMax)
            zMax = hitZ;

        if (hitZ < zMin)
            zMin = hitZ;
    }

    const int signPz((pHelixFit->GetMomentum().GetZ() > 0.f) ? 1 : -1);
    const float zStart((signPz > 0) ? zMin : zMax);
    const float zEnd((signPz > 0) ? zMax : zMin);

    pandora::CartesianVector startPosition(0.f, 0.f, 0.f), startMomentum(0.f, 0.f, 0.f);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(zStart, pHelixFit->GetReferencePoint(), startPosition));
    startMomentum = pHelixFit->GetExtrapolatedMomentum(startPosition);
    trackParameters.m_trackStateAtStart = pandora::TrackState(startPosition, startMomentum);
    
    pandora::CartesianVector endPosition(0.f, 0.f, 0.f), endMomentum(0.f, 0.f, 0.f);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pHelixFit->GetPointInZ(zEnd, pHelixFit->GetReferencePoint(), endPosition));
    endMomentum = pHelixFit->GetExtrapolatedMomentum(endPosition);
    trackParameters.m_trackStateAtEnd = pandora::TrackState(endPosition, endMomentum);

    this->GetECalProjectionOld(pHelixFit, signPz, trackParameters);
      
    delete pHelixFit;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void MyDemoAnalyzer::GetECalProjectionOld(const pandora::Helix *const pHelix, const int signPz, PandoraApi::Track::Parameters &trackParameters) const
{
  
    const pandora::CartesianVector &referencePoint(pHelix->GetReferencePoint());

    // First project to endcap
    float minGenericTime(std::numeric_limits<float>::max());
    bool isProjectedToEndCap(true);


    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    (void) pHelix->GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);


    // Then project to barrel surface(s)
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);

    if (m_eCalBarrelInnerSymmetry > 0)
    {
        // Polygon
        float twopi_n = 2. * M_PI / (static_cast<float>(m_eCalBarrelInnerSymmetry));


        for (int i = 0; i < m_eCalBarrelInnerSymmetry; ++i)
        {
            float genericTime(std::numeric_limits<float>::max());
            const float phi(twopi_n * static_cast<float>(i) + m_eCalBarrelInnerPhi0);

	    const pandora::StatusCode statusCode(pHelix->GetPointInXY(m_eCalBarrelInnerR * std::cos(phi), m_eCalBarrelInnerR * std::sin(phi),
                std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));



            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
	      minGenericTime = genericTime;
                isProjectedToEndCap = false;
                bestECalProjection = barrelProjection;
            }
        }
    }
    else
      {	  
        // Cylinder
        float genericTime(std::numeric_limits<float>::max());
        const pandora::StatusCode statusCode(pHelix->GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));


        if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
        {
            minGenericTime = genericTime;
            isProjectedToEndCap = false;
            bestECalProjection = barrelProjection;
        }
    }

    if (pandora::CartesianVector(0.f, 0.f, 0.f) == bestECalProjection)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(bestECalProjection, pHelix->GetExtrapolatedMomentum(bestECalProjection));
    trackParameters.m_isProjectedToEndCap = isProjectedToEndCap;



    // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
    //set it per default to Pions
    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);
    const float particleMass(trackParameters.m_mass.Get());
    const float particleEnergy(std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
    trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 299.792f;

}

//modify function a bit to get 
void MyDemoAnalyzer::GetTrackClusterDistance(const PandoraApi::Track::Parameters trackParameters, const EVENT::Cluster *const pCluster, const unsigned int maxSearchLayer, const float parallelDistanceCut, const float minTrackClusterCosAngle, float &trackClusterDistance) const
{
  const pandora::CartesianVector &trackPosition=trackParameters.m_trackStateAtCalorimeter.Get().GetPosition();
  const pandora::CartesianVector &trackDirection=trackParameters.m_trackStateAtCalorimeter.Get().GetMomentum().GetUnitVector(); 


  //in pandora cluster initialdirection is the unit-vector of the sum of associated calohit position unit vectors
  //define CartesianVector of the initial direction as unit vector of the  cluster position
  const float *pClusterPosition(pCluster->getPosition());
  double clusterposmag=sqrt(pClusterPosition[0]*pClusterPosition[0]+pClusterPosition[1]*pClusterPosition[1]+pClusterPosition[2]*pClusterPosition[2]);
  //position is NOT normalized but doesn't matter for the opening angle definition
  const pandora::CartesianVector clusterInitialDirection(pClusterPosition[0]/clusterposmag, pClusterPosition[1]/clusterposmag, pClusterPosition[2]/clusterposmag);
  if (trackDirection.GetCosOpeningAngle(clusterInitialDirection) >= minTrackClusterCosAngle){   
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());
    for(unsigned int i=0;i<pCluster->getCalorimeterHits().size();i++){
      static const int fCaloID   =    10 ;
      static const int fLayout   =  1000 ;
      static const int fLayer    = 10000 ;
      unsigned int type=pCluster->getCalorimeterHits()[i]->getType();
      unsigned int caloLayer=type/fLayer;
      int caloLayout=(type-caloLayer*fLayer)/fLayout;
      int caloID=(type-caloLayer*fLayer-caloLayout*fLayout)/fCaloID;
      if(caloID==1 && (caloLayer<=maxSearchLayer)){//caloID 1 is ecal
	if(caloLayer<=maxSearchLayer){
	  const pandora::CartesianVector calohitposition(pCluster->getCalorimeterHits()[i]->getPosition()[0],pCluster->getCalorimeterHits()[i]->getPosition()[1],pCluster->getCalorimeterHits()[i]->getPosition()[2]);
	  const pandora::CartesianVector positionDifference(calohitposition - trackPosition);
	  
	  if (std::fabs(trackDirection.GetDotProduct(positionDifference)) > parallelDistanceCut){
	    continue;
	  }
	  const float perpendicularDistanceSquared((trackDirection.GetCrossProduct(positionDifference)).GetMagnitudeSquared());
	  if (perpendicularDistanceSquared < minDistanceSquared)
	    {
	      minDistanceSquared = perpendicularDistanceSquared;
	      distanceFound = true;
	    }
	}
      }
    }
    if(distanceFound){
      trackClusterDistance = std::sqrt(minDistanceSquared);
    }
  }
}

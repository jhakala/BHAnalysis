// -*- C++ -*-
//
// Package:    BHAnalyzerTLBSM
// Class:      BHAnalyzerTLBSM
// 
/**\class BHAnalyzerTLBSM BHAnalyzerTLBSM.cc AnalysisCodeTLBSM/BHAnalyzerTLBSM/src/BHAnalyzerTLBSM.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Alexey Ferapontov,8 R-021,+41227676332,
//         Created:  Mon Apr 2 11:25:01 CEST 2010
// $Id: BHAnalyzerTLBSM.cc,v 1.1 2012/04/03 03:08:26 aferapon Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"

// ECAL spike cleaning - Swiss cross
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//

using namespace edm;
using namespace std;

class BHAnalyzerTLBSM : public edm::EDAnalyzer {
   public:
      explicit BHAnalyzerTLBSM(const edm::ParameterSet&);
      ~BHAnalyzerTLBSM();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;  
      
      void createHistogram(const std::string& folderName);
      bool check(std::string process, std::string pCheck);
      void init(const edm::TriggerResults &, const edm::TriggerNames & HLTNames);

      std::map<std::string, unsigned int> prescales;
      std::map<std::string, unsigned int> prescale_counter; 
      std::map<std::string, unsigned int> trigger_indices;

      edm::InputTag eleLabel_;
      edm::InputTag muoLabel_;
      edm::InputTag jetLabel_;
      edm::InputTag tauLabel_;
      edm::InputTag metLabel_;
      edm::InputTag phoLabel_;
      edm::InputTag rechitBLabel_;
      edm::InputTag rechitELabel_;
      edm::InputTag pvSrc_;
      edm::InputTag triggerLabel_;

      // tools for clusters
      // std::auto_ptr<EcalClusterLazyTools> lazyTools_;
              
      bool isMCBH;

      edm::InputTag ebRecHitsLabel_;
      edm::InputTag eeRecHitsLabel_;

      edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsToken_;
      edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsToken_;

      reco::TrackBase::TrackQuality _trackQuality;
 
      std::vector< std::map<std::string,TH1*> > histos_; 
      std::vector<std::string> cutNames_;
      //TH1F* h_norm;

      //TTree
      TTree* tree;
      float JetE[25];
      float JetPx[25];
      float JetPy[25];
      float JetPz[25];
      float JetPt[25];
      float JetEta[25];
      float JetPhi[25];      
      
      // JetEMF was made from emEnergyFraction()--this does not work with our PAT jets
      // float JetEMF[25];
      
      float EleE[25];
      float ElePx[25];
      float ElePy[25];
      float ElePz[25];
      float ElePt[25];
      float EleEta[25];
      float ElePhi[25];      

      float PhE[25];
      float PhPx[25];
      float PhPy[25];
      float PhPz[25];
      float PhPt[25];
      float PhEta[25];
      float PhPhi[25];
      float PhSwissCross[25];      

      float MuE[25];
      float MuPx[25];
      float MuPy[25];
      float MuPz[25];
      float MuPt[25];
      float MuEta[25];
      float MuPhi[25];    
      float MuDxy[25];  
      
      float ST;
      float mBH;
      float Met;
      float MetE;
      float MetPx;
      float MetPy;
      float MetPz;
      float MetPt;      
      float MetPhi;      
      float Sphericity;
      float JetArr[4];      
      float EleArr[4];
      float MuArr[4];      
      float PhArr[4];
      int NPV;
      int NTracks;
      int NJets;
      int NElectrons;
      int NPhotons;      
      int NMuons; 
      bool NoScrap;     
      
      int Multiplicity; 
      bool isLeptonPhoton;
      bool isEleChannel;
      bool isMuChannel;
      bool isPhChannel;

      //
      float ResJetEta;
      float ResJetPhi;      
      float ResJetM;
      float ResJetPt;
      float ResEleEta;
      float ResElePhi;      
      float ResEleM;
      float ResElePt;      
      float ResPhEta;
      float ResPhPhi;      
      float ResPhM;
      float ResPhPt;
      float ResMuEta;
      float ResMuPhi;      
      float ResMuM;
      float ResMuPt;
      float ResLepEta;
      float ResLepPhi;      
      float ResLepM;
      float ResLepPt;
      float ResObjEta;
      float ResObjPhi;      
      float ResObjM;
      float ResObjPt;
      
      int runno;
      int evtno;
      int lumiblock;
      int isRealData;
      float muon_d0;
                                                                                      
      float LeadingArr[4];    

      //HLT info (Jets and MET)
      bool firedHLT_L1Jet6U;
      bool firedHLT_L1Jet10U;
      bool firedHLT_Jet15U;
      bool firedHLT_Jet30U;
      bool firedHLT_Jet50U;
      bool firedHLT_DiJetAve15U_8E29;
      bool firedHLT_DiJetAve30U_8E29;
      bool firedHLT_MET45;
      bool firedHLT_MET100;
      bool firedHLT_HT100U;
      bool firedHLT_HT140U;
      bool firedHLT_HT140U_Eta3_v1;
      bool firedHLT_HT160U_v1;
      bool firedHLT_HT200U; 
      bool firedHLT_HT200U_v1;     
      bool firedHLT_FwdJet20U;
      bool firedHLT_QuadJet15U;
      bool firedHLT_L1MET20;
      
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BHAnalyzerTLBSM::BHAnalyzerTLBSM(const edm::ParameterSet& iConfig):
  eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  phoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag")),
  pvSrc_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertex")),
  triggerLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag")),  
  isMCBH(iConfig.getUntrackedParameter<bool>("MCLabel",false)),
  ebRecHitsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("ebRecHitTag")),
  eeRecHitsLabel_(iConfig.getUntrackedParameter<edm::InputTag>("eeRecHitTag")),
  ebRecHitsToken_(consumes<EcalRecHitCollection>(ebRecHitsLabel_)),
  eeRecHitsToken_(consumes<EcalRecHitCollection>(eeRecHitsLabel_))
{
   //now do what ever initialization is needed

}


BHAnalyzerTLBSM::~BHAnalyzerTLBSM()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BHAnalyzerTLBSM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   runno = iEvent.id().run();
   evtno  = iEvent.id().event();
   lumiblock = iEvent.luminosityBlock();
   isRealData = iEvent.isRealData();
   
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // first: get all objects from the event.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   
   edm::Handle<edm::View<pat::Muon> > muonHandle;
   iEvent.getByLabel(muoLabel_,muonHandle);
   const edm::View<pat::Muon> & muons = *muonHandle;   // const ... &, we don't make a copy of it!
   
   edm::Handle<edm::View<pat::Jet> > jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);
   const edm::View<pat::Jet> & jets = *jetHandle;
   
   edm::Handle<edm::View<pat::Electron> > electronHandle;
   iEvent.getByLabel(eleLabel_,electronHandle);
   const edm::View<pat::Electron> & electrons = *electronHandle;
   
   edm::Handle<edm::View<pat::MET> > metHandle;
   iEvent.getByLabel(metLabel_,metHandle);
   const edm::View<pat::MET> & mets = *metHandle;
   
   edm::Handle<edm::View<pat::Photon> > phoHandle;
   iEvent.getByLabel(phoLabel_,phoHandle);
   const edm::View<pat::Photon> & photons = *phoHandle;

   edm::Handle<edm::View<pat::Tau> > tauHandle;
   iEvent.getByLabel(tauLabel_,tauHandle);
   //const edm::View<pat::Tau> & taus = *tauHandle;

   // Swiss cross - ECAL cleaning
   edm::View<pat::Photon>::const_iterator photon;
   Handle<EcalRecHitCollection> Brechit;//barrel
   Handle<EcalRecHitCollection> Erechit;//endcap
   iEvent.getByLabel(rechitBLabel_,Brechit);
   iEvent.getByLabel(rechitELabel_,Erechit);
   //const EcalRecHitCollection* barrelRecHits= Brechit.product();
   //const EcalRecHitCollection* endcapRecHits= Erechit.product();

   edm::ESHandle<CaloTopology> pTopology;
   edm::ESHandle<CaloTopology> theCaloTopo_;
   iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
   //const CaloTopology *topology = theCaloTopo_.product();
   
   // Loop over all objects and evaluate BH properties: ST, mBH, multiplicity ...
   ST=0.;
   
   //Lorentz vectors for jets, egamma, muons, leptons, and all objects
   math::XYZTLorentzVectorF pJet;
   math::XYZTLorentzVectorF pEle;      
   math::XYZTLorentzVectorF pMu;
   math::XYZTLorentzVectorF pPh;
   math::XYZTLorentzVectorF pLep;
   math::XYZTLorentzVectorF pObj;           

   //Beam spot
   math::XYZPoint bs;
       
   // MET
   math::XYZTLorentzVectorF pBH;
   pBH += mets[0].p4();
   Met = mets[0].pt();
   MetPhi = mets[0].phi();
   ST += Met;
   
   MetE = mets[0].energy();
   MetPx = mets[0].px();
   MetPy = mets[0].py();
   MetPz = mets[0].pz();
   MetPt = mets[0].pt();
      
   //Sphericity
   float sumPx2=mets[0].px()*mets[0].px();
   float sumPy2=mets[0].py()*mets[0].py();
   float sumPxPy=mets[0].px()*mets[0].py();   
   
   //Leading objects
   std::vector<double> leadingJets (4,-1);
   std::vector<double> leadingElectrons (4,-1);
   std::vector<double> leadingPhotons (4,-1);
   std::vector<double> leadingMuons (4,-1);
   std::vector<double> leadingObjects (4,-1);
   std::vector<double> leadingLeptons (4,-1);
      
   float leadingElePt = 0;
   float leadingMuPt  = 0;
   float leadingPhPt  = 0;

   int jetcnt         = -1;
   int elecnt         = -1;
   //int mucnt        = -1;
   //int phcnt        = -1;
   int ngoodmuons     = 0;
   int ngoodphotons   = 0;
   int ngoodelectrons = 0;
   int ngoodjets      = 0;   
   int nPVcount       = 0;
   
   //------ Primary Vertices     
   edm::Handle< reco::VertexCollection > PVCollection; 
   if (iEvent.getByLabel(pvSrc_, PVCollection )) {
     for (reco::VertexCollection::const_iterator pv = PVCollection->begin(); pv != PVCollection->end(); ++pv) {
       //--- vertex selection
       //std::cout<<" PV "<<pv->x()<<" "<<pv->y()<<" "pv->z()<<std::endl;
       //std::cout<<"PV parameters "<<pv->isFake()<<" "<<fabs(pv->z())<<" "<<pv->position().Rho()<<std::endl;
       if (!pv->isFake() && pv->ndof() > 4 && fabs(pv->z()) <= 15. && pv->position().Rho() <= 2.) ++nPVcount;       
     }
   }
   
   // No scraping
   bool noscrap = true;
   
   // HLT results   
   firedHLT_L1Jet6U = false;
   firedHLT_L1Jet10U = false;
   firedHLT_Jet15U = false;
   firedHLT_Jet30U = false;
   firedHLT_Jet50U = false;
   firedHLT_DiJetAve15U_8E29 = false;
   firedHLT_DiJetAve30U_8E29 = false;
   firedHLT_MET45 = false;
   firedHLT_MET100 = false;
   firedHLT_HT100U = false;
   firedHLT_HT140U = false;
   firedHLT_HT140U_Eta3_v1 = false;
   firedHLT_HT160U_v1 = false;
   firedHLT_HT200U = false; 
   firedHLT_FwdJet20U = false;
   firedHLT_QuadJet15U = false;
   firedHLT_L1MET20 = false;
            
   TriggerResults tr;
   Handle<TriggerResults> h_trigRes;
   iEvent.getByLabel(triggerLabel_, h_trigRes);

   tr = *h_trigRes;
  
   std::vector<string> triggerList;
   Service<service::TriggerNamesService> tns;
   bool foundNames = tns->getTrigPaths(tr,triggerList);
   if (!foundNames) std::cout << "Could not get trigger names!\n";
   if (tr.size()!=triggerList.size()) std::cout << "ERROR: length of names and paths not the same: " << triggerList.size() << "," << tr.size() << endl;
   // dump trigger list at first event
   for (unsigned int i=0; i< tr.size(); i++) {
     if ( !tr[i].accept() == 1 ) continue;
     if( triggerList[i] == "HLT_L1Jet6U") { firedHLT_L1Jet6U = true; }
     if( triggerList[i] == "HLT_L1Jet10U") { firedHLT_L1Jet10U = true; }    
     if( triggerList[i] == "HLT_Jet15U") { firedHLT_Jet15U = true; }
     if( triggerList[i] == "HLT_Jet30U") { firedHLT_Jet30U = true; }
     if( triggerList[i] == "HLT_Jet50U") { firedHLT_Jet50U = true; }  
     if( triggerList[i] == "HLT_DiJetAve15U_8E29") { firedHLT_DiJetAve15U_8E29 = true; }   
     if( triggerList[i] == "HLT_DiJetAve30U_8E29") { firedHLT_DiJetAve30U_8E29 = true; }
     if( triggerList[i] == "HLT_MET45") { firedHLT_MET45 = true; }  
     if( triggerList[i] == "HLT_MET100") { firedHLT_MET100 = true; }  
     if( triggerList[i] == "HLT_HT100U") { firedHLT_HT100U = true; }
     if( triggerList[i] == "HLT_HT140U") { firedHLT_HT140U = true; }
     if( triggerList[i] == "HLT_HT140U_Eta3_v1") { firedHLT_HT140U_Eta3_v1 = true; }
     if( triggerList[i] == "HLT_HT160U_v1") { firedHLT_HT160U_v1 = true; }
     if( triggerList[i] == "HLT_HT200U") { firedHLT_HT200U = true; }
     if( triggerList[i] == "HLT_HT200U_v1") { firedHLT_HT200U_v1 = true; }
     if( triggerList[i] == "HLT_FwdJet20U") { firedHLT_FwdJet20U = true; }
     if( triggerList[i] == "HLT_QuadJet15U") { firedHLT_QuadJet15U = true; }
     if( triggerList[i] == "HLT_L1MET20") { firedHLT_L1MET20 = true; }

   }

   for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){     
     
     ++jetcnt; 
     pBH += jet->p4();
     ST += jet->pt();
     sumPx2 += jet->px()*jet->px();
     sumPy2 += jet->py()*jet->py();
     sumPxPy += jet->px()*jet->py();      
      
     leadingJets.push_back(jet->pt());
     leadingObjects.push_back(jet->pt());

     pJet += jet->p4();
     pObj += jet->p4(); 
           
     JetE[jetcnt] = jet->energy();
     JetPx[jetcnt] = jet->px();
     JetPy[jetcnt] = jet->py();
     JetPz[jetcnt] = jet->pz();     
     JetPt[jetcnt] = jet->pt();
     JetEta[jetcnt] = jet->eta();      
     JetPhi[jetcnt] = jet->phi();
    // Calling emEnergyFraction() throws the error "This PAT jet was not made from a CaloJet." 
    // JetEMF[jetcnt] = jet->emEnergyFraction();
     
     ++ngoodjets;
     
   }   

   
   for(edm::View<pat::Electron>::const_iterator e = electrons.begin(); e!=electrons.end(); ++e){

     ++elecnt;
          
     pBH += e->p4();
     ST += e->pt();
     sumPx2 += e->px()*e->px();
     sumPy2 += e->py()*e->py();
     sumPxPy += e->px()*e->py();
     if (e->pt() > leadingElePt) leadingElePt = e->pt();
     
     leadingElectrons.push_back(e->pt());
     leadingLeptons.push_back(e->pt());
     leadingObjects.push_back(e->pt());                       
     
     pEle += e->p4();
     pLep += e->p4();      
     pObj += e->p4(); 

     EleE[elecnt] = e->energy();
     ElePx[elecnt] = e->px();
     ElePy[elecnt] = e->py();
     ElePz[elecnt] = e->pz();     
     ElePt[elecnt] = e->pt();
     EleEta[elecnt] = e->eta();      
     ElePhi[elecnt] = e->phi(); 
     ++ngoodelectrons;
                
   }  
   // ecal information
   
   // lazyTools_ = std::auto_ptr<EcalClusterLazyTools>( new EcalClusterLazyTools(iEvent,iSetup,ebRecHitsToken_,eeRecHitsToken_) );

   // get ecal barrel recHits for spike rejection
   edm::Handle<EcalRecHitCollection> recHitsEB_h;
   iEvent.getByLabel(ebRecHitsLabel_, recHitsEB_h );
   // const EcalRecHitCollection * recHitsEB = 0;
   // if ( ! recHitsEB_h.isValid() ) {
   //   LogError("ExoDiPhotonAnalyzer") << " ECAL Barrel RecHit Collection not available !"; return;
   // } else {
   //   recHitsEB = recHitsEB_h.product();
   // }
         
   //PAT photons     
   for(edm::View<pat::Photon>::const_iterator ph = photons.begin(); ph!=photons.end(); ++ph){

     /*
     double e4 = 	lazyTools_->eTop(*((*ph).superCluster()->seed())) 	+ 
     			lazyTools_->eBottom(*((*ph).superCluster()->seed())) 	+
     			lazyTools_->eLeft(*((*ph).superCluster()->seed()))	+
     			lazyTools_->eRight(*((*ph).superCluster()->seed()))	;
     			
     double eMax = ph->maxEnergyXtal();  
     double spikeSelector = 1-e4/eMax;
      
     if(spikeSelector > 0.95 && fabs(ph->eta())<1.560) {
       std::cout <<"Swiss cross variable "<<spikeSelector<<" run "<<runno<<" event "<<evtno<<std::endl;
       std::cout <<"This photon candidate is an ECAL spike identified by Swiss Cross algorithm. Removing it from photons"<<std::endl;
       continue;
     }
     */

     if (abs(ph->eta()) > 2.4) continue;
          	
     //If not a spike, increment # photons
     ++ngoodphotons;
     
     pBH += ph->p4();
     ST += ph->pt();
     sumPx2 += ph->px()*ph->px();
     sumPy2 += ph->py()*ph->py();
     sumPxPy += ph->px()*ph->py();
     if (ph->pt() > leadingPhPt) leadingPhPt = ph->pt();
     
     leadingPhotons.push_back(ph->pt());
     leadingObjects.push_back(ph->pt());      
     
     pPh += ph->p4();    
     pObj += ph->p4();      
     
     PhE[ngoodphotons-1] = ph->energy();
     PhPx[ngoodphotons-1] = ph->px();
     PhPy[ngoodphotons-1] = ph->py();
     PhPz[ngoodphotons-1] = ph->pz();
     PhPt[ngoodphotons-1] = ph->pt();
     PhEta[ngoodphotons-1] = ph->eta();      
     PhPhi[ngoodphotons-1] = ph->phi();
     //PhSwissCross[ngoodphotons-1] = spikeSelector;
     
   }
  
   for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
     
     // Muons don't have an innerTrack() method in miniAOD?
     // muon_d0 = fabs(mu->innerTrack()->dxy(bs));
     // if (fabs(mu->innerTrack()->dxy(bs))>0.2) continue;

     // If a good muon (config muon ID + dxy), increment # muons
     ++ngoodmuons;
     
     pBH += mu->p4();
     ST += mu->pt();
     sumPx2 += mu->px()*mu->px();
     sumPy2 += mu->py()*mu->py();
     sumPxPy += mu->px()*mu->py();
     if (mu->pt() > leadingMuPt) leadingMuPt = mu->pt();
     
     leadingMuons.push_back(mu->pt());
     leadingLeptons.push_back(mu->pt());
     leadingObjects.push_back(mu->pt());
     
     pMu += mu->p4();
     pLep += mu->p4();      
     pObj += mu->p4(); 
     
     MuE[ngoodmuons-1] = mu->energy();
     MuPx[ngoodmuons-1] = mu->px();
     MuPy[ngoodmuons-1] = mu->py();
     MuPz[ngoodmuons-1] = mu->pz(); 
     MuPt[ngoodmuons-1] = mu->pt();
     MuEta[ngoodmuons-1] = mu->eta();      
     MuPhi[ngoodmuons-1] = mu->phi(); 
     MuDxy[ngoodmuons-1] = fabs(mu->innerTrack()->dxy(bs));

   }
   
   for (int i=0;i<25;++i) {
     if (i>=ngoodjets) {
       JetE[i]=0.;
       JetPx[i]=0.;
       JetPy[i]=0.;
       JetPz[i]=0.;
       JetPt[i]=0.;
       JetEta[i]=99.;
       JetPhi[i]=99.;
       // JetEMF[i]=99.;
     }
     if (i>=ngoodelectrons) {
       EleE[i]=0.;
       ElePx[i]=0.;
       ElePy[i]=0.;
       ElePz[i]=0.;
       ElePt[i]=0.;
       EleEta[i]=99.;
       ElePhi[i]=99.;
     }
     if (i>=ngoodphotons) {
       PhE[i]=0.;
       PhPx[i]=0.;
       PhPy[i]=0.;
       PhPz[i]=0.;
       PhPt[i]=0.;
       PhEta[i]=99.;
       PhPhi[i]=99.;
       PhSwissCross[i]=99.;
     }       
     if (i>=ngoodmuons) {
       MuE[i]=0.;
       MuPx[i]=0.;
       MuPy[i]=0.;
       MuPz[i]=0.;
       MuPt[i]=0.;
       MuEta[i]=99.;
       MuPhi[i]=99.;
       MuDxy[i]=99.;
     }         
   }
   
   //Sorting
   std::sort(leadingJets.begin(), leadingJets.end());
   std::sort(leadingElectrons.begin(), leadingElectrons.end());
   std::sort(leadingPhotons.begin(), leadingPhotons.end());
   std::sort(leadingMuons.begin(), leadingMuons.end());
   std::sort(leadingLeptons.begin(), leadingLeptons.end());
   std::sort(leadingObjects.begin(), leadingObjects.end());

   std::reverse(leadingJets.begin(), leadingJets.end());
   std::reverse(leadingElectrons.begin(), leadingElectrons.end());
   std::reverse(leadingPhotons.begin(), leadingPhotons.end());
   std::reverse(leadingMuons.begin(), leadingMuons.end());
   std::reverse(leadingLeptons.begin(), leadingLeptons.end());
   std::reverse(leadingObjects.begin(), leadingObjects.end());

   for (int i=0; i<4; ++i){
     JetArr[i] = leadingJets[i];
     EleArr[i] = leadingElectrons[i];
     PhArr[i] = leadingPhotons[i];
     MuArr[i] = leadingMuons[i];
     LeadingArr[i] = leadingObjects[i];
   } 

   float trace = sumPx2 + sumPy2;
   float det = sumPx2*sumPy2 - sumPxPy*sumPxPy;
   float lambda2 = (trace - sqrt(trace*trace - 4*det))/2.0;
   Sphericity = 2*lambda2/trace;
   
   mBH = 0.;
   
   mBH = pBH.M();
   
   //h_norm -> Fill(ST);
   for (size_t i=0; i<cutNames_.size(); ++i){
     std::map<std::string,TH1*>& histo = histos_[i];       
     histo["ST"]->Fill(ST);
   }   
   
   NPV = nPVcount;
   NTracks = 0;
   NJets = ngoodjets;
   NElectrons = ngoodelectrons;
   NMuons = ngoodmuons; 	//not using muons.size() here as |dxy(bs)| < 0.2 might remove cosmics
   NPhotons = ngoodphotons; 	//not using photons.size() here as Swiss cross might remove spikes
   Multiplicity = ngoodjets + ngoodelectrons + ngoodphotons + ngoodmuons;
   
   // Attempts to reduce the size of output file and reject noise and unused events (Mult < 2)
   //if (isMCBH == false && Multiplicity < 2) return;
          
   NoScrap = noscrap; 
   isLeptonPhoton = (ngoodelectrons + ngoodmuons + ngoodphotons > 0);
   
   // Classify event by leading lepton/photon
   isEleChannel = (isLeptonPhoton && (leadingElePt > leadingMuPt) && (leadingElePt > leadingPhPt));
   isMuChannel = (isLeptonPhoton && (leadingMuPt > leadingElePt) && (leadingMuPt > leadingPhPt));
   isPhChannel = (isLeptonPhoton && (leadingPhPt > leadingElePt) && (leadingPhPt > leadingMuPt));
   
   //Filling objects "averaged" angles
   // 99. is the default for angles
   ResJetEta = 99.;
   ResJetPhi = 99.;
   ResEleEta = 99.;
   ResElePhi = 99.;   
   ResPhEta = 99.;
   ResPhPhi = 99.;     
   ResMuEta = 99.;
   ResMuPhi = 99.;
   ResLepEta = 99.;
   ResLepPhi = 99.;      
   ResObjEta = 99.;
   ResObjPhi = 99.; 
    
   if (NJets) {
     ResJetEta = pJet.Eta();
     ResJetPhi = pJet.Phi();
   } else {
     ResJetEta = 99.;
     ResJetPhi = 99.;
   }          
   ResJetM = pJet.M();
   ResJetPt = pJet.Pt();

   if (NElectrons) {
     ResEleEta = pEle.Eta();
     ResElePhi = pEle.Phi();
   } else {
     ResEleEta = 99.;
     ResElePhi = 99.;
   }   
   ResEleM = pEle.M();
   ResElePt = pEle.Pt();
   
   if (NPhotons) {
     ResPhEta = pPh.Eta();
     ResPhPhi = pPh.Phi();
   } else {
     ResPhEta = 99.;
     ResPhPhi = 99.;
   }           
   ResPhM = pPh.M();
   ResPhPt = pPh.Pt();
   
   if (NMuons) {
     ResMuEta = pMu.Eta();
     ResMuPhi = pMu.Phi();
   } else {
     ResMuEta = 99.;
     ResMuPhi = 99.;
   }           
   ResMuM = pMu.M();
   ResMuPt = pMu.Pt();
   
   if (NMuons+NElectrons) {
     ResLepEta = pLep.Eta();
     ResLepPhi = pLep.Phi();
   } else {
     ResLepEta = 99.;
     ResLepPhi = 99.;
   }       
   ResLepM = pLep.M();
   ResLepPt = pLep.Pt();
   
   if (NMuons+NElectrons+NPhotons+NJets) {
     ResObjEta = pObj.Eta();
     ResObjPhi = pObj.Phi();
   } else {
     ResObjEta = 99.;
     ResObjPhi = 99.;
   }     
   ResObjM = pObj.M();
   ResObjPt = pObj.Pt();
   
   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void BHAnalyzerTLBSM::beginJob()
{
  cutNames_.push_back("No_Cut");
  
  tree = fs_->make<TTree>("t","t");
  //h_norm = new TH1F("h_norm","",500,0,50000);

  tree->Branch("JetE",&JetE,"JetE[25]");
  tree->Branch("JetPx",&JetPx,"JetPx[25]");
  tree->Branch("JetPy",&JetPy,"JetPy[25]");
  tree->Branch("JetPz",&JetPz,"JetPz[25]");
  tree->Branch("JetPt",&JetPt,"JetPt[25]");
  tree->Branch("JetEta",&JetEta,"JetEta[25]");
  tree->Branch("JetPhi",&JetPhi,"JetPhi[25]");
  // tree->Branch("JetEMF",&JetEMF,"JetEMF[25]");  

  tree->Branch("EleE",&EleE,"EleE[25]");
  tree->Branch("ElePx",&ElePx,"ElePx[25]");
  tree->Branch("ElePy",&ElePy,"ElePy[25]");
  tree->Branch("ElePz",&ElePz,"ElePz[25]");
  tree->Branch("ElePt",&ElePt,"ElePt[25]");
  tree->Branch("EleEta",&EleEta,"EleEta[25]");
  tree->Branch("ElePhi",&ElePhi,"ElePhi[25]");
    
  tree->Branch("PhE",&PhE,"PhE[25]");
  tree->Branch("PhPx",&PhPx,"PhPx[25]");
  tree->Branch("PhPy",&PhPy,"PhPy[25]");
  tree->Branch("PhPz",&PhPz,"PhPz[25]");
  tree->Branch("PhPt",&PhPt,"PhPt[25]");
  tree->Branch("PhEta",&PhEta,"PhEta[25]");
  tree->Branch("PhPhi",&PhPhi,"PhPhi[25]");
  tree->Branch("PhSwissCross",&PhSwissCross,"PhSwissCross[25]");  

  tree->Branch("MuE",&MuE,"MuE[25]");
  tree->Branch("MuPx",&MuPx,"MuPx[25]");
  tree->Branch("MuPy",&MuPy,"MuPy[25]");
  tree->Branch("MuPz",&MuPz,"MuPz[25]");
  tree->Branch("MuPt",&MuPt,"MuPt[25]");
  tree->Branch("MuEta",&MuEta,"MuEta[25]");
  tree->Branch("MuPhi",&MuPhi,"MuPhi[25]");
  tree->Branch("MuDxy",&MuDxy,"MuDxy[25]");  
    
  tree->Branch("ST",&ST,"ST/F");
  tree->Branch("mBH",&mBH,"mBH/F");
  tree->Branch("Met",&Met,"Met/F");
  tree->Branch("MetE",&MetE,"MetE/F");
  tree->Branch("MetPx",&MetPx,"MetPx/F");
  tree->Branch("MetPy",&MetPy,"MetPy/F");
  tree->Branch("MetPz",&MetPz,"MetPz/F");
  tree->Branch("MetPt",&MetPt,"MetPt/F");
  tree->Branch("MetPhi",&MetPhi,"MetPhi/F");   
  tree->Branch("Sphericity", &Sphericity, "Sphericity/F");
  tree->Branch("Jet",JetArr,"JetArr[4]");   
  tree->Branch("Ele",EleArr,"EleArr[4]");
  tree->Branch("Mu",MuArr,"MuArr[4]");
  tree->Branch("Ph",PhArr,"PhArr[4]");   
  tree->Branch("NPV",&NPV,"NPV/I");   
  tree->Branch("NTracks",&NTracks,"NTracks/I");   
  tree->Branch("NJets",&NJets,"NJets/I");
  tree->Branch("NElectrons",&NElectrons,"NElectrons/I");
  tree->Branch("NPhotons",&NPhotons,"NPhotons/I");
  tree->Branch("NMuons",&NMuons,"NMuons/I");
  tree->Branch("NoScrap", &NoScrap, "NoScrap/B");     
  
  tree->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
  tree->Branch("isLeptonPhoton", &isLeptonPhoton, "isLeptonPhoton/B");   
  tree->Branch("isEleChannel", &isEleChannel, "isEleChannel/B");
  tree->Branch("isMuChannel", &isMuChannel, "isMuChannel/B"); 
  tree->Branch("isPhChannel", &isPhChannel, "isPhChannel/B");
  
  tree->Branch("ResJetEta",&ResJetEta,"ResJetEta/F"); 
  tree->Branch("ResJetPhi",&ResJetPhi,"ResJetPhi/F");
  tree->Branch("ResJetM",&ResJetM,"ResJetM/F");
  tree->Branch("ResJetPt",&ResJetPt,"ResJetPt/F");
  tree->Branch("ResEleEta",&ResEleEta,"ResEleEta/F"); 
  tree->Branch("ResElePhi",&ResElePhi,"ResElePhi/F");
  tree->Branch("ResEleM",&ResEleM,"ResEleM/F");
  tree->Branch("ResElePt",&ResElePt,"ResElePt/F");  
  tree->Branch("ResPhEta",&ResPhEta,"ResPhEta/F"); 
  tree->Branch("ResPhPhi",&ResPhPhi,"ResPhPhi/F");
  tree->Branch("ResPhM",&ResPhM,"ResPhM/F");
  tree->Branch("ResPhPt",&ResPhPt,"ResPhPt/F");   
  tree->Branch("ResMuEta",&ResMuEta,"ResMuEta/F"); 
  tree->Branch("ResMuPhi",&ResMuPhi,"ResMuPhi/F");
  tree->Branch("ResMuM",&ResMuM,"ResMuM/F");
  tree->Branch("ResMuPt",&ResMuPt,"ResMuPt/F");
  tree->Branch("ResLepEta",&ResLepEta,"ResLepEta/F"); 
  tree->Branch("ResLepPhi",&ResLepPhi,"ResLepPhi/F");
  tree->Branch("ResLepM",&ResLepM,"ResLepM/F");
  tree->Branch("ResLepPt",&ResLepPt,"ResLepPt/F");   
  tree->Branch("ResObjEta",&ResObjEta,"ResObjEta/F"); 
  tree->Branch("ResObjPhi",&ResObjPhi,"ResObjPhi/F");
  tree->Branch("ResObjM",&ResObjM,"ResObjM/F");
  tree->Branch("ResObjPt",&ResObjPt,"ResObjPt/F"); 
  
  tree->Branch("runno",&runno,"runno/I");
  tree->Branch("evtno",&evtno,"evtno/I"); 
  tree->Branch("lumiblock",&lumiblock,"lumiblock/I");               
  tree->Branch("isRealData",&isRealData,"isRealData/I");
  tree->Branch("muon_d0",&muon_d0,"muon_d0/F");                                                      
  
  tree->Branch("firedHLT_L1Jet6U",&firedHLT_L1Jet6U,"firedHLT_L1Jet6U/B");
  tree->Branch("firedHLT_L1Jet10U",&firedHLT_L1Jet10U,"firedHLT_L1Jet10U/B");
  tree->Branch("firedHLT_Jet15U",&firedHLT_Jet15U,"firedHLT_Jet15U/B");
  tree->Branch("firedHLT_Jet30U",&firedHLT_Jet30U,"firedHLT_Jet30U/B");
  tree->Branch("firedHLT_Jet50U",&firedHLT_Jet50U,"firedHLT_Jet50U/B");
  tree->Branch("firedHLT_DiJetAve15U_8E29",&firedHLT_DiJetAve15U_8E29,"firedHLT_DiJetAve15U_8E29/B");
  tree->Branch("firedHLT_DiJetAve30U_8E29",&firedHLT_DiJetAve30U_8E29,"firedHLT_DiJetAve30U_8E29/B");
  tree->Branch("firedHLT_MET45",&firedHLT_MET45,"firedHLT_MET45/B");
  tree->Branch("firedHLT_MET100",&firedHLT_MET100,"firedHLT_MET100/B");
  tree->Branch("firedHLT_HT100U",&firedHLT_HT100U,"firedHLT_HT100U/B");  
  tree->Branch("firedHLT_HT140U",&firedHLT_HT140U,"firedHLT_HT140U/B");
  tree->Branch("firedHLT_HT140U_Eta3_v1",&firedHLT_HT140U_Eta3_v1,"firedHLT_HT140U_Eta3_v1/B");
  tree->Branch("firedHLT_HT160U_v1",&firedHLT_HT160U_v1,"firedHLT_HT160U_v1/B");
  tree->Branch("firedHLT_HT200U",&firedHLT_HT200U,"firedHLT_HT200U/B");
  tree->Branch("firedHLT_HT200U_v1",&firedHLT_HT200U_v1,"firedHLT_HT200U_v1/B");
  tree->Branch("firedHLT_FwdJet20U",&firedHLT_FwdJet20U,"firedHLT_FwdJet20U/B");
  tree->Branch("firedHLT_QuadJet15U",&firedHLT_QuadJet15U,"firedHLT_QuadJet15U/B");
  tree->Branch("firedHLT_L1MET20",&firedHLT_L1MET20,"firedHLT_L1MET20/B");
    
  for (size_t i=0; i<cutNames_.size(); ++i)
    createHistogram(cutNames_[i]);
  
}

void BHAnalyzerTLBSM::createHistogram(const std::string& folderName){
  TFileDirectory subDir = fs_->mkdir(folderName);
  std::map<std::string, TH1*> container;
  
  container["ST"] = subDir.make<TH1F>("ST", "ST", 500, 0, 100000);
  
  histos_.push_back(container);          
}

bool BHAnalyzerTLBSM::check(std::string process, std::string pCheck) {
  bool value= false;
  
  if (process == pCheck) value = true;
  
  return value;
}

// ------------ method called once each job just after ending the event loop  ------------
void BHAnalyzerTLBSM::endJob() {
  //h_norm->Write();
}


//define this as a plug-in
DEFINE_FWK_MODULE(BHAnalyzerTLBSM);

// -*- C++ -*-
//
// Package:    BHAnalyzerTLBSM
// Class:      BHAnalyzerPATTuples
// 
/**\class BHAnalyzerPATTuplesTLBSM BHAnalyzerPATTuplesTLBSM.cc AnalysisCodeTLBSM/BHAnalyzerTLBSM/src/BHAnalyzerPATTuplesTLBSM.cc
   
Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Alexey Ferapontov,8 R-021,+41227676332,
//         Created:  Mon Apr 2 11:25:01 CEST 2010
// $Id: BHAnalyzerPATTuplesTLBSM.cc,v 1.1 2012/04/03 03:08:26 aferapon Exp $
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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

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

class BHAnalyzerPATTuplesTLBSM : public edm::EDAnalyzer {
public:
  explicit BHAnalyzerPATTuplesTLBSM(const edm::ParameterSet&);
  ~BHAnalyzerPATTuplesTLBSM();
  
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
  edm::InputTag pvSrc_;
  edm::InputTag triggerLabel_;
  bool isMCBH;
  double threshold_; 
  bool DEBUG_;
  bool DEBUGTRG_; 
  
  //BH::ObjectIDSelector selector_;
  //BH::JetResidualCorrector corrector_;
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
  float JetEMF[25];
  
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
  //float PhSwissCross[25];      
  
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
  //float JetArr[4];      
  //float EleArr[4];
  //float MuArr[4];      
  //float PhArr[4];
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
  //New method - with wildcards
  bool firedHLT_HT;  
  bool firedHLT_HT100;
  bool firedHLT_HT150;
  bool firedHLT_HT200;
  bool firedHLT_HT250;
  bool firedHLT_HT300;
  bool firedHLT_HT350;
  bool firedHLT_HT400;
  bool firedHLT_HT450;
  bool firedHLT_HT500;
  bool firedHLT_HT550;
  bool firedHLT_HT600;
  bool firedHLT_HT650;
  bool firedHLT_HT700;
  bool firedHLT_HT750;
  bool firedHLT_HT800;
  bool firedHLT_HT850;
  bool firedHLT_PFHT300;
  bool firedHLT_PFHT350;
  bool firedHLT_PFHT400;
  bool firedHLT_PFHT650;
  bool firedHLT_PFHT700;
  bool firedHLT_PFHT750;
  
  int rechits;
   
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
BHAnalyzerPATTuplesTLBSM::BHAnalyzerPATTuplesTLBSM(const edm::ParameterSet& iConfig):
  eleLabel_(iConfig.getUntrackedParameter<edm::InputTag>("electronTag")),
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  phoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("photonTag")),
  pvSrc_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertex")),
  triggerLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag")),  
  isMCBH(iConfig.getUntrackedParameter<bool>("MCLabel",false)),
  threshold_(iConfig.getUntrackedParameter<double>("PtThreshold")),
  DEBUG_(iConfig.getUntrackedParameter<bool>("DEBUG",false)),
  DEBUGTRG_(iConfig.getUntrackedParameter<bool>("DEBUGTRG",false))
{
  //now do what ever initialization is needed
  
}


BHAnalyzerPATTuplesTLBSM::~BHAnalyzerPATTuplesTLBSM()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BHAnalyzerPATTuplesTLBSM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  
  edm::Handle<View<reco::Track> >  trackHandle;
  iEvent.getByLabel("generalTracks",trackHandle);
  const edm::View<reco::Track> & tracks = *trackHandle;
  
  //Rechits
  //edm::Handle<HcalNoiseSummary> hSummary;
  //iEvent.getByType(hSummary);
  //rechits = hSummary->GetRecHitCount();

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

  float sumPx2 = 0;
  float sumPy2 = 0;
  float sumPxPy = 0; 
      
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
  int ntracks        = 0;
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
      if (!pv->isFake() && pv->ndof() > 4 && fabs(pv->z()) <= 24. && pv->position().Rho() <= 2.) ++nPVcount;
      //cout <<"vertex "<<pv->x()<<" "<<pv->y()<<" "<<pv->z()<<" "<<endl;      
    }
  }
  
  // No scraping
  bool noscrap = false;
  float fraction = 0;
  
  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks",tkRef);    
  const reco::TrackCollection* tkColl = tkRef.product();
  
  int numhighpurity=0;
  _trackQuality = reco::TrackBase::qualityByName("highPurity");
  
  if(tkColl->size()>10){ 
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      //std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
      if(itk->quality(_trackQuality)) numhighpurity++;
    }
    fraction = (float)numhighpurity/(float)tkColl->size();
    if(fraction>0.25) noscrap=true;
  }else{
    //if less than 10 Tracks accept the event anyway    
    noscrap = true;
  }
  
  // HLT results    
  //New method - with wildcards
  firedHLT_HT 	    	= false;
  firedHLT_HT100	= false;
  firedHLT_HT150	= false;
  firedHLT_HT200	= false;
  firedHLT_HT250	= false;
  firedHLT_HT300	= false;
  firedHLT_HT350	= false;
  firedHLT_HT400	= false;
  firedHLT_HT450	= false;
  firedHLT_HT500	= false;
  firedHLT_HT550	= false;
  firedHLT_HT600	= false;
  firedHLT_HT650	= false;
  firedHLT_HT700	= false;
  firedHLT_HT750	= false;
  firedHLT_HT800	= false;
  firedHLT_HT850	= false;
  firedHLT_PFHT300	= false;
  firedHLT_PFHT350	= false;
  firedHLT_PFHT400	= false;
  firedHLT_PFHT650	= false;
  firedHLT_PFHT700	= false;
  firedHLT_PFHT750	= false;
    
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
  // int ht_trig_fired = 0;
  // int rsq_trig_fired = 0;
 
  for (unsigned int i=0; i< tr.size(); i++) {
    if ( !tr[i].accept() == 1 ) continue;    
    if (DEBUGTRG_) cout<<"Trigger fired is "<<triggerList[i]<<endl;

    //New method - with wildcards    
    if (triggerList[i].find("HT") != std::string::npos) {
      firedHLT_HT = true; 
      if (DEBUGTRG_) cout<<"       HT/PFHT fired "<<triggerList[i]<<endl;
    }
    if (triggerList[i].find("HLT_HT100") != std::string::npos) { firedHLT_HT100 = true; }
    if (triggerList[i].find("HLT_HT150") != std::string::npos) { firedHLT_HT150 = true; }
    if (triggerList[i].find("HLT_HT200") != std::string::npos) { firedHLT_HT200 = true; }
    if (triggerList[i].find("HLT_HT250") != std::string::npos) { firedHLT_HT250 = true; }
    if (triggerList[i].find("HLT_HT300") != std::string::npos) { firedHLT_HT300 = true; }
    if (triggerList[i].find("HLT_HT350") != std::string::npos) { firedHLT_HT350 = true; }
    if (triggerList[i].find("HLT_HT400") != std::string::npos) { firedHLT_HT400 = true; }
    if (triggerList[i].find("HLT_HT450") != std::string::npos) { firedHLT_HT450 = true; }
    if (triggerList[i].find("HLT_HT500") != std::string::npos) { firedHLT_HT500 = true; }
    if (triggerList[i].find("HLT_HT550") != std::string::npos) { firedHLT_HT550 = true; }
    if (triggerList[i].find("HLT_HT600") != std::string::npos) { firedHLT_HT600 = true; }
    if (triggerList[i].find("HLT_HT650") != std::string::npos) { firedHLT_HT650 = true; }
    if (triggerList[i].find("HLT_HT700") != std::string::npos) { firedHLT_HT700 = true; }
    if (triggerList[i].find("HLT_HT750") != std::string::npos) { firedHLT_HT750 = true; }
    if (triggerList[i].find("HLT_HT800") != std::string::npos) { firedHLT_HT800 = true; }
    if (triggerList[i].find("HLT_HT850") != std::string::npos) { firedHLT_HT850 = true; }
    if (triggerList[i].find("HLT_PFHT300") != std::string::npos) { firedHLT_PFHT300 = true; }
    if (triggerList[i].find("HLT_PFHT350") != std::string::npos) { firedHLT_PFHT350 = true; }
    if (triggerList[i].find("HLT_PFHT400") != std::string::npos) { firedHLT_PFHT400 = true; }
    if (triggerList[i].find("HLT_PFHT650") != std::string::npos) { firedHLT_PFHT650 = true; }
    if (triggerList[i].find("HLT_PFHT700") != std::string::npos) { firedHLT_PFHT700 = true; }
    if (triggerList[i].find("HLT_PFHT750") != std::string::npos) { firedHLT_PFHT750 = true; }      
  }
    
  // get the primary vertex
  reco::Vertex primaryVertex;   
  primaryVertex = PVCollection->at(0); // No need to check if PV is valid, since it is required in the config
  
  // get the beamspot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if (!iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpotHandle)) { cout<<"No beam spot "<<endl; return;} 
  
  //Jets
  std::auto_ptr<pat::JetCollection> outputJets(new pat::JetCollection());  
  for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
    if (jet->pt() > threshold_ && fabs(jet->eta()) < 2.6) {
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
      //JetEMF[jetcnt] = jet->emEnergyFraction();
      
      ++ngoodjets;
    }   
  }
  if (DEBUG_) cout<<"NJets "<<ngoodjets<<endl;
  
  //Electrons
  std::auto_ptr<pat::ElectronCollection> outputElectrons(new pat::ElectronCollection());   
  for(edm::View<pat::Electron>::const_iterator e = electrons.begin(); e!=electrons.end(); ++e){  
    if (e->pt() > threshold_ && fabs(e->eta()) < 2.4) {
      
      const float ElectronD0Cut_ = 0.04, EtaThr_e = 2.5, ElectronETSCThr_ = 15., RemuThr = 0.1, ElectronVertexMatchThr_ = 1., RelIso_e = 0.2;
      float D0_e = 99.;
      // float Chi2_e = 99., 
      float ET_SC = -1.;
      
      bool hadId((int)e->electronID("eidTight") & 0x1);
      bool isNotConv((int)e->electronID("eidTight") & 0x4); 
      bool isGsfElectron = true;
      if (!e->gsfTrack()) isGsfElectron = false;
      int nlosthits   = e->gsfTrack()->numberOfLostHits();
      double convDist = e->convDist();
      double convDcot = e->convDcot();
      D0_e      = e->gsfTrack()->dxy(beamSpotHandle->position());
      ET_SC = e->superCluster()->energy() * sin ( e->superCluster()->position().theta() );
      double Iso = (e->chargedHadronIso()+e->neutralHadronIso()+e->photonIso())/e->pt();
      
      if (DEBUG_) {
	cout<<"isGsf "<<isGsfElectron<<", hadID "<<hadId<<", not converted? "<<isNotConv<<", ecalseed? "<<e->ecalDrivenSeed()<<endl;    
	cout<<"N lost hits "<<nlosthits<<", cut < 2"<<endl;
	cout<<"convDcot "<<fabs(convDcot)<<", convDist "<<fabs(convDist)<<", cut > 0.02"<<endl;
	cout<<"El D0 "<<fabs(D0_e)<<", cut < "<< ElectronD0Cut_<<"; alternative e->dB(pat::Electron::BS2D) "<<e->dB(pat::Electron::BS2D)<<endl;
	cout<<"ET_SC "<<ET_SC<<", cut > "<<ElectronETSCThr_<<endl;
	cout<<"El relative iso "<<Iso<<", cut < "<<RelIso_e<<endl;    
	cout<<"El vertex x/y/z "<<e->vx()<<" "<<e->vy()<<" "<<e->vz()<<", PV x/y/z "<<primaryVertex.x()<<" "<<primaryVertex.y()<<" "<<primaryVertex.z()<<", beamspot x/y/z "<<beamSpotHandle->position().x()<<" "<<beamSpotHandle->position().y()<<" "<<beamSpotHandle->position().z()<<endl;
      }
      
      if(!isGsfElectron || !hadId || !isNotConv)                  		continue; 	      
      if(nlosthits >= 2 || (fabs(convDcot) < 0.02 && fabs(convDist) < 0.02)) 	continue;	
      if(e->ecalDrivenSeed()!= 1)                     				continue;
      if(fabs(D0_e)       >= ElectronD0Cut_)         				continue; 
      if(fabs(e->eta())   >= EtaThr_e)               				continue;
      if(ET_SC            <= ElectronETSCThr_)       				continue;
      if (Iso             > RelIso_e )               				continue;
      
      float dRmue;    
      bool SharedCone = false; //DR computed between localElectrons candidates and all Global or Tracker Muon
      for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
	if ((mu->isTrackerMuon()) || (mu->isGlobalMuon())){
	  dRmue = deltaR(e->eta(),e->phi(),mu->eta(),mu->phi());
	  if (dRmue < RemuThr){
	    cout<<"Electron overlaps with global or tracker muon"<<endl;
	    SharedCone = true;
	    break;
	  }
	}
      }
      
      if ( SharedCone )                      continue;
      if ( fabs( e->vz() - primaryVertex.z() )  > ElectronVertexMatchThr_ ) continue;
      
      float dRej = 0.;
      int overlap_ej = 0;    
      for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
        dRej = deltaR(jet->eta(),jet->phi(),e->eta(),e->phi());
        if (DEBUG_) cout<<"Run "<<runno<<", Event"<<evtno<<" , Lumi "<<lumiblock<<", dR electron-jet "<<dRej<<", electron pt "<<e->pt()<<", jet pt "<<jet->pt()<<", electron eta "<<e->eta()<<", jet eta "<<jet->eta()<<", electron phi "<<e->phi()<<", jet phi "<<jet->phi()<<endl;
        if (dRej < 0.3) ++overlap_ej;
      }
      if (DEBUG_) cout << "overlap e jet "<<overlap_ej<<endl;
      //if (overlap_ej > 0) continue; //this does not apply
      
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
  }  
  
  //Photons 
  std::auto_ptr<pat::PhotonCollection> outputPhotons(new pat::PhotonCollection());
  for(edm::View<pat::Photon>::const_iterator ph = photons.begin(); ph!=photons.end(); ++ph){
    if (ph->pt() > threshold_ && fabs(ph->eta()) < 2.4) {
      
      //No need to code it here - it's being selected in the configs. If ID is removed from configs, uncomment the lines! 
      /*
	if (ph->hadronicOverEm() >= 0.05) continue; 
	if (ph->ecalRecHitSumEtConeDR04() >= 4.2 + 0.006 * ph->pt()) continue;
	if (ph->hcalTowerSumEtConeDR04() >= 2.2 + 0.0025 * ph->pt()) continue;
	if (ph->trkSumPtHollowConeDR04() >= 2.0 + 0.001 * ph->pt()) continue;
	if (ph->hasPixelSeed()) continue;
	if (ph->sigmaIetaIeta() >= (ph->isEB() ? 0.013 : 0.030)) continue;
      */
      
      float dRphj = 0.;
      int overlap_phj = 0;    
      for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
        dRphj = deltaR(jet->eta(),jet->phi(),ph->eta(),ph->phi());
        if (DEBUG_) cout<<"Run "<<runno<<", Event"<<evtno<<" , Lumi "<<lumiblock<<", dR photon-jet "<<dRphj<<", photon pt "<<ph->pt()<<", jet pt "<<jet->pt()<<", photon eta "<<ph->eta()<<", jet eta "<<jet->eta()<<", photon phi "<<ph->phi()<<", jet phi "<<jet->phi()<<endl;
        if (dRphj < 0.3) ++overlap_phj;
      }
      if (DEBUG_) cout << "overlap ph jet "<<overlap_phj<<endl;
      if (overlap_phj > 0) continue;
      
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
      //PhSwissCross[ngoodphotons-1] = 0;
    }
  } 
  
  //Muons
  std::auto_ptr<pat::MuonCollection> outputMuons(new pat::MuonCollection());
  for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
    if (mu->pt() > threshold_ && fabs(mu->eta()) < 2.1) {      
      if (!(mu->isTrackerMuon()) || !(mu->isGlobalMuon())) continue;
      
      const float MuonRelIso = 0.02, MuonD0Cut = 0.02, MuonVertexMatchThr = 1., MuonNofValidHits = 0, MuonNofValidTrHits = 10., MuonNormChi2 = 10.;
      float Chi2 = 99., D0Inner = 99., D0Standard = 99., RelIso03PF = 99.;
      int NValidHits = -1, NTrValidHits = -1;
      
      // Get the tracker track from the muon       
      const reco::TrackRef globalTrack = mu->globalTrack();
      if (globalTrack.isNonnull()) {
	// D0   = globalTrack->dxy(beamSpotHandle->position());
	Chi2 = globalTrack->normalizedChi2();
      }
      if (mu->innerTrack().isNonnull()) D0Inner = mu->innerTrack()->dxy (beamSpotHandle->position());
      if (mu->isGlobalMuon())        NValidHits = mu->globalTrack()->hitPattern().numberOfValidMuonHits();
      NTrValidHits = mu->innerTrack()->numberOfValidHits();
      D0Standard = mu->dB(pat::Muon::BS2D); //Much more elegant way to find D0. Compare it to D0Inner!
      
      float muvz = mu->vz();
      //float pvx = primaryVertex.x();
      //float pvy = primaryVertex.y();
      float pvz = primaryVertex.z();
      
      //Isolation: 
      float PATNeutralHadronIso =  mu->neutralHadronIso();
      float PATChargedHadronIso =  mu->chargedHadronIso();
      float PATPhotonIso        =  mu->photonIso();
      // float PATTrackIso         =  mu->trackIso(); 
      RelIso03PF                = ((PATNeutralHadronIso+PATChargedHadronIso+PATPhotonIso)/mu->pt());
      
      if (DEBUG_) {
	std::cout << "Chi2: " << Chi2 << " Treshold: " << MuonNormChi2 <<  std::endl;
	std::cout << "NTrValidHits: " << NTrValidHits << " Treshold: " << MuonNofValidTrHits << std::endl;
	std::cout << "NValidHits: " << NValidHits << " Treshold: " << MuonNofValidHits << std::endl;
	std::cout << "D0Inner: " << fabs(D0Inner) << " Treshold: " << MuonD0Cut << " D0 standard "<<fabs(D0Standard)<<std::endl;
	std::cout << "Eta: " << fabs( mu->eta() ) << std::endl;
	std::cout << "Pt: " << mu->pt() << " Treshold: " << threshold_ << std::endl;
	cout <<" muon z: "<<muvz<<" Primary Vertex z: "<<pvz<<endl;    
	std::cout << "Iso: " << RelIso03PF << " Threshold: " << MuonRelIso << std::endl;
      }
      
      //Selection
      if (Chi2                	>= MuonNormChi2)        continue;
      if (NTrValidHits        	<= MuonNofValidTrHits) 	continue;
      if (NValidHits          	<= MuonNofValidHits)    continue;   
      if (fabs(D0Standard)      > MuonD0Cut)	        continue;
      if (std::abs(muvz - pvz) 	> MuonVertexMatchThr)	continue;   		 
      if (RelIso03PF 		> MuonRelIso) 		continue;
      
      //Overlaps checker
      float dRmuj;
      int overlap_muj = 0;
      for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){
        dRmuj = deltaR(mu->eta(),mu->phi(),jet->eta(),jet->phi());
        if (DEBUG_) cout<<"Run "<<runno<<", Event"<<evtno<<" , Lumi "<<lumiblock<<", dR muon-jet "<<dRmuj<<", muon pt "<<mu->pt()<<", jet pt "<<jet->pt()<<", muon eta "<<mu->eta()<<", jet eta "<<jet->eta()<<", muon phi "<<mu->phi()<<", jet phi "<<jet->phi()<<" muon d0 "<<fabs(D0Standard)<<endl;
        if (dRmuj < 0.3) ++overlap_muj;
      }
      if (DEBUG_) cout<<"Overlap muon-jet "<<overlap_muj<<endl;
      //Need to think, if we want to remove such muons...
      
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
   }  
  if (DEBUG_) cout<<"Nmuons "<<muons.size()<<" good "<<ngoodmuons<<endl;
  
  //Tracks
  for(edm::View<reco::Track>::const_iterator trk = tracks.begin(); trk!=tracks.end(); ++trk){
    ++ntracks;
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
      JetEMF[i]=99.;
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
      //PhSwissCross[i]=99.;
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
  /*  
  //Sorting jets explicitly here after residual corrections!
  float tmpE, tmpPt, tmpPx, tmpPy, tmpPz, tmpEMF, tmpEta, tmpPhi;
  int g,r;
    for (r = 0; r < 24; r++) {
      for (g = r+1; g < 25; g++) {
        if (JetPt[r] < JetPt[g]) {
	  tmpPt = JetPt[r]; JetPt[r] = JetPt[g]; JetPt[g] = tmpPt;
	  tmpPx = JetPx[r]; JetPx[r] = JetPx[g]; JetPx[g] = tmpPx;
	  tmpPy = JetPy[r]; JetPy[r] = JetPy[g]; JetPy[g] = tmpPy;
	  tmpPz = JetPz[r]; JetPz[r] = JetPz[g]; JetPz[g] = tmpPz;
	  tmpE  = JetE[r]; JetE[r] = JetE[g]; JetE[g] = tmpE;
	  tmpEMF = JetEMF[r]; JetEMF[r] = JetEMF[g]; JetEMF[g] = tmpEMF;
	  tmpEta = JetEta[r]; JetEta[r] = JetEta[g]; JetEta[g] = tmpEta;
	  tmpPhi = JetPhi[r]; JetPhi[r] = JetPhi[g]; JetPhi[g] = tmpPhi;
	  }
        }   
      } 
  */  
  /*for (int i = 0; i < 4; ++i){
    JetArr[i] = leadingJets[i];
    EleArr[i] = leadingElectrons[i];
    PhArr[i] = leadingPhotons[i];
    MuArr[i] = leadingMuons[i];
    LeadingArr[i] = leadingObjects[i];
  }*/ 
  
  float trace = sumPx2 + sumPy2;
  float det = sumPx2*sumPy2 - sumPxPy*sumPxPy;
  float lambda2 = (trace - sqrt(trace*trace - 4*det))/2.0;
  Sphericity = 2*lambda2/trace;
  
  NPV = nPVcount;
  NTracks = ntracks;
  NJets = ngoodjets;
  NElectrons = ngoodelectrons;
  NMuons = ngoodmuons; 	//not using muons.size() here as |dxy(bs)| < 0.2 might remove cosmics
  NPhotons = ngoodphotons; 	//not using photons.size() here as Swiss cross might remove spikes
  Multiplicity = ngoodjets + ngoodelectrons + ngoodphotons + ngoodmuons;
  
  //MET
  Met = mets[0].pt();
  MetPhi = mets[0].phi();
  MetE = mets[0].energy();
  MetPx = mets[0].px();
  MetPy = mets[0].py();
  MetPz = mets[0].pz();
  MetPt = mets[0].pt();
  //pBH += mets[0].p4();
  
  //Sphericity
  sumPx2=mets[0].px()*mets[0].px();
  sumPy2=mets[0].py()*mets[0].py();
  sumPxPy=mets[0].px()*mets[0].py();
  
  //removing met muon corrections if no muons were found
  //if (NMuons == 0) {
  //cout<<"Muons "<<NMuons<<" MET_def "<<Met<<" MetUncorr "<< mets[0].uncorrectedPt(pat::MET::uncorrMUON)<<endl;
  //  Met = mets[0].uncorrectedPt(pat::MET::uncorrMUON);
  //  MetPhi = mets[0].uncorrectedPhi(pat::MET::uncorrMUON);
  //  MetPt = mets[0].uncorrectedPt(pat::MET::uncorrMUON);
  //} 
  
  //cout<<"ST before MET "<<ST<<endl;
  if (Met > threshold_) ST += Met;  
  //cout <<"MET "<<Met<<" ST "<<ST<<endl;       
  mBH = 0.;  
  mBH = pBH.M(); //without MET now
  
  //h_norm -> Fill(ST);
  for (size_t i=0; i<cutNames_.size(); ++i){
    std::map<std::string,TH1*>& histo = histos_[i];       
    histo["ST"]->Fill(ST);
  }
  
  // Attempts to reduce the size of output file and reject noise and unused events (Mult < 2)
  if (isMCBH == false && Multiplicity < 2) return;
  //if (ST < 1000) return;
  
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
void BHAnalyzerPATTuplesTLBSM::beginJob()
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
  tree->Branch("JetEMF",&JetEMF,"JetEMF[25]");  

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
  //tree->Branch("PhSwissCross",&PhSwissCross,"PhSwissCross[25]");  

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
  //tree->Branch("Jet",JetArr,"JetArr[4]");   
  //tree->Branch("Ele",EleArr,"EleArr[4]");
  //tree->Branch("Mu",MuArr,"MuArr[4]");
  //tree->Branch("Ph",PhArr,"PhArr[4]");   
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
  
  tree->Branch("firedHLT_HT",&firedHLT_HT,"firedHLT_HT/B");  
  tree->Branch("firedHLT_HT100",&firedHLT_HT100,"firedHLT_HT100/B");
  tree->Branch("firedHLT_HT150",&firedHLT_HT150,"firedHLT_HT150/B");
  tree->Branch("firedHLT_HT200",&firedHLT_HT200,"firedHLT_HT200/B");
  tree->Branch("firedHLT_HT250",&firedHLT_HT250,"firedHLT_HT250/B");
  tree->Branch("firedHLT_HT300",&firedHLT_HT300,"firedHLT_HT300/B");
  tree->Branch("firedHLT_HT350",&firedHLT_HT350,"firedHLT_HT350/B");
  tree->Branch("firedHLT_HT400",&firedHLT_HT400,"firedHLT_HT400/B");
  tree->Branch("firedHLT_HT450",&firedHLT_HT450,"firedHLT_HT450/B");
  tree->Branch("firedHLT_HT500",&firedHLT_HT500,"firedHLT_HT500/B");
  tree->Branch("firedHLT_HT550",&firedHLT_HT550,"firedHLT_HT550/B");
  tree->Branch("firedHLT_HT600",&firedHLT_HT600,"firedHLT_HT600/B");
  tree->Branch("firedHLT_HT650",&firedHLT_HT650,"firedHLT_HT650/B");
  tree->Branch("firedHLT_HT700",&firedHLT_HT700,"firedHLT_HT700/B");
  tree->Branch("firedHLT_HT750",&firedHLT_HT750,"firedHLT_HT750/B");
  tree->Branch("firedHLT_HT800",&firedHLT_HT800,"firedHLT_HT800/B");
  tree->Branch("firedHLT_HT850",&firedHLT_HT850,"firedHLT_HT850/B");
  tree->Branch("firedHLT_PFHT300",&firedHLT_PFHT300,"firedHLT_PFHT300/B");
  tree->Branch("firedHLT_PFHT350",&firedHLT_PFHT350,"firedHLT_PFHT350/B");
  tree->Branch("firedHLT_PFHT400",&firedHLT_PFHT400,"firedHLT_PFHT400/B");
  tree->Branch("firedHLT_PFHT650",&firedHLT_PFHT650,"firedHLT_PFHT650/B");
  tree->Branch("firedHLT_PFHT700",&firedHLT_PFHT700,"firedHLT_PFHT700/B");
  tree->Branch("firedHLT_PFHT750",&firedHLT_PFHT750,"firedHLT_PFHT750/B");

  tree->Branch("rechits",&rechits,"rechits/I");       

  for (size_t i=0; i<cutNames_.size(); ++i)
    createHistogram(cutNames_[i]);
  
}

void BHAnalyzerPATTuplesTLBSM::createHistogram(const std::string& folderName){
  TFileDirectory subDir = fs_->mkdir(folderName);
  std::map<std::string, TH1*> container;
  
  container["ST"] = subDir.make<TH1F>("ST", "ST", 500, 0, 100000);
  
  histos_.push_back(container);          
}

bool BHAnalyzerPATTuplesTLBSM::check(std::string process, std::string pCheck) {
  bool value= false;
  
  if (process == pCheck) value = true;
  
  return value;
}

// ------------ method called once each job just after ending the event loop  ------------
void BHAnalyzerPATTuplesTLBSM::endJob() {
  //h_norm->Write();
}


//define this as a plug-in
DEFINE_FWK_MODULE(BHAnalyzerPATTuplesTLBSM);

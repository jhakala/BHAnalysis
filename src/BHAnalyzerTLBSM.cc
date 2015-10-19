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


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

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

#include "DataFormats/Common/interface/ValueMap.h"





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
     virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;  

      void createHistogram(const std::string& folderName);
      bool check(std::string process, std::string pCheck);
      void init(const edm::TriggerResults &, const edm::TriggerNames & HLTNames);

      std::map<std::string, unsigned int> prescales;
      std::map<std::string, unsigned int> prescale_counter; 
      std::map<std::string, unsigned int> trigger_indices;

      edm::InputTag muoLabel_;
      edm::InputTag jetLabel_;
      edm::InputTag tauLabel_;
      edm::InputTag metLabel_;
      edm::InputTag rechitBLabel_;
      edm::InputTag rechitELabel_;
      edm::InputTag pvSrc_;
      edm::InputTag triggerLabel_;
      edm::InputTag filterLabel_;
      edm::InputTag rhoLabel_;
      

              
      bool isMCBH;
      bool DEBUG_;


      edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsToken_;
      edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsToken_;

      reco::TrackBase::TrackQuality _trackQuality;
 
      std::vector< std::map<std::string,TH1*> > histos_; 
      std::vector<std::string> cutNames_;
     
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_; 
      edm::EDGetToken eleLabelToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;
      
      //Electron Decisions
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      //Photon Decisions
      edm::EDGetToken phoLabelToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_; 
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      //TTree
      TTree* tree;
      float JetE[25];
      float JetPx[25];
      float JetPy[25];
      float JetPz[25];
      float JetPt[25];
      float JetEt[25];
      float JetEta[25];
      float JetPhi[25];      
      
      float EleE[25];
      float ElePx[25];
      float ElePy[25];
      float ElePz[25];
      float ElePt[25];
      float EleEt[25];
      float EleEta[25];
      float ElePhi[25]; 

      float PhE[25];
      float PhPx[25];
      float PhPy[25];
      float PhPz[25];
      float PhPt[25];
      float PhEt[25];
      float PhEta[25];
      float PhPhi[25];
      
      float MuE[25];
      float MuPx[25];
      float MuPy[25];
      float MuPz[25];
      float MuPt[25];
      float MuEt[25];
      float MuEta[25];
      float MuPhi[25];    
      
      float ST;
      float mBH;
      float Met;
      float MetE;
      float MetPx;
      float MetPy;
      float MetPz;
      float MetPt;      
      float MetEt;      
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
      float ResJetEt;
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
  //bool firedHLT_PFJet60_v2;
  //bool firedHLT_PFJet140_v2;
  //bool firedHLT_PFJet450_v2;
  //bool firedHLT_PFHT300_v1;
  //bool firedHLT_PFHT400_v1;
  bool firedHLT_PFHT475_v2;
  //bool firedHLT_PFHT600_v2;
  //bool firedHLT_PFHT650_v2;
  bool firedHLT_PFHT800_v2;


  // bool passed_HBHENoiseFilter;
  // bool passed_HBHENoiseIsoFilter;
  bool passed_CSCTightHaloFilter;
  // bool passed_hcalLaserEventFilter;
  bool passed_EcalDeadCellTriggerPrimitiveFilter;
  bool passed_EcalDeadCellBoundaryEnergyFilter;
  bool passed_goodVertices;
  bool passed_eeBadScFilter;
  // bool passed_ecalLaserCorrFilter;
  // bool passed_trkPOGFilters;
  // bool passed_trkPOG_manystripclus53X;
  // bool passed_trkPOG_toomanystripclus53X;
  // bool passed_trkPOG_logErrorTooManyClusters;
  bool passed_METFilters;
  //TODO 
  double Reliso_el;
  double Reliso_mu;
  int pass_eleID_medium;      
  
  Float_t rho_;

  };


  // constants, enums and typedefs

  // static data member definitions

  // constructors and destructor
  BHAnalyzerTLBSM::BHAnalyzerTLBSM(const edm::ParameterSet& iConfig):
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  pvSrc_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertex")),
  triggerLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag")),  
  filterLabel_(iConfig.getUntrackedParameter<edm::InputTag>("filterTag")),  
  rhoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rho_lable")),
  isMCBH(iConfig.getUntrackedParameter<bool>("MCLabel",false)),
  DEBUG_(iConfig.getUntrackedParameter<bool>("DEBUG",false)),
  
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"))),
  phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"))),
  phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
  {
   //now do what ever initialization is needed
   beamSpotToken_            = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
   eleLabelToken_            = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronTag"));
   vtxMiniAODToken_          = mayConsume<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("verticesMiniAOD"));
   conversionsMiniAODToken_  = mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversionsMiniAOD"));
   phoLabelToken_            = mayConsume<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"));
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
  // All PF Candidate for alternate isolation
  
   // Get the beam spot
   edm::Handle<reco::BeamSpot> theBeamSpot;
   iEvent.getByToken(beamSpotToken_,theBeamSpot);  

   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(eleLabelToken_, electrons);
  
  
   //Get the conversions collection
   edm::Handle<reco::ConversionCollection> conversions;
   iEvent.getByToken(conversionsMiniAODToken_, conversions);	

   // Get the electron ID data from the event stream.
   //   // Note: this implies that the VID ID modules have been run upstream.
   // If you need more info, check with the EGM group.
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
   
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);

   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);

   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);

  edm::Handle<edm::ValueMap<bool> > loose_id_decisions_ph;
  iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions_ph);
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions_ph;
  iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions_ph);
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions_ph;
  iEvent.getByToken(phoTightIdMapToken_ ,tight_id_decisions_ph);

   edm::Handle<double> rhoHandle;
   iEvent.getByLabel(rhoLabel_,rhoHandle);
   
   if(rhoHandle.isValid()) {
   rho_ = *(rhoHandle.product());
   }

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muoLabel_,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;   

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  const edm::View<pat::Jet> & jets = *jetHandle;
   
   edm::Handle<edm::View<pat::MET> > metHandle;
   iEvent.getByLabel(metLabel_,metHandle);
   const edm::View<pat::MET> & mets = *metHandle;
   
   edm::Handle<edm::View<reco::Photon> > photons;
   iEvent.getByToken(phoLabelToken_,photons);


   edm::Handle<edm::View<pat::Tau> > tauHandle;
   iEvent.getByLabel(tauLabel_,tauHandle);
   
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);
   // Swiss cross - ECAL cleaning
   edm::View<pat::Photon>::const_iterator photon;
   Handle<EcalRecHitCollection> Brechit;//barrel
   Handle<EcalRecHitCollection> Erechit;//endcap
   iEvent.getByLabel(rechitBLabel_,Brechit);
   iEvent.getByLabel(rechitELabel_,Erechit);

   edm::ESHandle<CaloTopology> pTopology;
   edm::ESHandle<CaloTopology> theCaloTopo_;
   iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
   
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
   MetEt = mets[0].et();
      
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
   int mucnt          = -1;
   int phcnt          = -1;
   int ngoodmuons     = 0;
   int ngoodphotons   = 0;
   int ngoodelectrons = 0;
   int ngoodjets      = 0;   
   int nPVcount       = 0;
   
   //------ Primary Vertices------- 
   edm::Handle< reco::VertexCollection > PVCollection; 
   iEvent.getByLabel(pvSrc_, PVCollection);
   const reco::Vertex & vertex_ = PVCollection->front();
   if (PVCollection->empty()) return; // skip the event if no PV found

   if (iEvent.getByLabel(pvSrc_, PVCollection )) {
     for (reco::VertexCollection::const_iterator pv = PVCollection->begin(); pv != PVCollection->end(); ++pv) {
       //--- vertex selection
       //std::cout<<" PV "<<pv->x()<<" "<<pv->y()<<" "pv->z()<<std::endl;
       //std::cout<<"PV parameters "<<pv->isFake()<<" "<<fabs(pv->z())<<" "<<pv->position().Rho()<<std::endl;
       //std::cout<<"PV position "<<pv->position()<<endl;
       if (!pv->isFake() && pv->ndof() > 4 && fabs(pv->z()) <= 24. && pv->position().Rho() <= 2.) ++nPVcount;       
     }
   }
   
   // No scraping
   bool noscrap = true;
   
   // HLT results   
   //firedHLT_PFJet60_v2  = false;
   //firedHLT_PFJet140_v2 = false;
   //firedHLT_PFJet450_v2 = false;
   //firedHLT_PFHT300_v1  = false;
   //firedHLT_PFHT400_v1  = false;
   firedHLT_PFHT475_v2  = false;
   //firedHLT_PFHT600_v2  = false;
   //firedHLT_PFHT650_v2  = false;
   firedHLT_PFHT800_v2  = false;
   //TODO
   
   TriggerResults tr;
   Handle<TriggerResults> h_trigRes;
   iEvent.getByLabel(triggerLabel_, h_trigRes);

   // MET filter results   
   // passed_HBHENoiseFilter = false;
   // passed_HBHENoiseIsoFilter = false;
   passed_CSCTightHaloFilter = false;
   // passed_hcalLaserEventFilter = false;
   passed_EcalDeadCellTriggerPrimitiveFilter = false;
   passed_EcalDeadCellBoundaryEnergyFilter = false;
   passed_goodVertices = false;
   passed_eeBadScFilter = false;
   // passed_ecalLaserCorrFilter = false;
   // passed_trkPOGFilters = false;
   // passed_trkPOG_manystripclus53X = false;
   // passed_trkPOG_toomanystripclus53X = false;
   // passed_trkPOG_logErrorTooManyClusters = false;
   passed_METFilters = false;

   TriggerResults fr;
   Handle<TriggerResults> h_filtRes;
   iEvent.getByLabel(filterLabel_, h_filtRes);
   fr = *h_filtRes;
  
   std::vector<string> triggerList;
   Service<service::TriggerNamesService> tns;
   bool foundNames = tns->getTrigPaths(tr,triggerList);
   if (!foundNames) std::cout << "Could not get trigger names!\n";
   if (tr.size()!=triggerList.size()) std::cout << "ERROR: length of names and paths not the same: " 
                                                << triggerList.size() << "," << tr.size() << endl;
   // dump trigger list at first event
   for (unsigned int i=0; i< tr.size(); i++) {
     //std::cout << "["<<i<<"] = " << triggerList[i]<<setw(40)<<
     //": Prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (tr[i].accept() ? "Event Passed" : "Event Failed") << endl;
     if ( !tr[i].accept() == 1 ) continue;
     //if( triggerList[i] == "HLT_PFJet60_v2")  { firedHLT_PFJet60_v2  = true; }
     //if( triggerList[i] == "HLT_PFJet140_v2") { firedHLT_PFJet140_v2 = true; }
     //if( triggerList[i] == "HLT_PFJet450_v2") { firedHLT_PFJet450_v2 = true; }
     //if( triggerList[i] == "HLT_PFHT300_v1")  { firedHLT_PFHT300_v1  = true; }
     //if( triggerList[i] == "HLT_PFHT400_v1")  { firedHLT_PFHT400_v1  = true; }
     if( triggerList[i] == "HLT_PFHT475_v2")  { firedHLT_PFHT475_v2  = true; }
     //if( triggerList[i] == "HLT_PFHT600_v2")  { firedHLT_PFHT600_v2  = true; }
     //if( triggerList[i] == "HLT_PFHT650_v2")  { firedHLT_PFHT650_v2  = true; }
     if( triggerList[i] == "HLT_PFHT800_v2")  { firedHLT_PFHT800_v2  = true; }

   }
   std::vector<string> filterList;
   Service<service::TriggerNamesService> fns;
   bool foundFilterNames = fns->getTrigPaths(fr,filterList);
   if (!foundFilterNames) std::cout << "Could not get filter names!\n";
   if (fr.size()!=filterList.size()) std::cout << "ERROR: length of filter names and paths not the same: " 
                                                << filterList.size() << "," << fr.size() << endl;
   // dump filter list at first event
   for (unsigned int i=0; i< fr.size(); i++) {
     //std::cout << "["<<i<<"] = " << filterList[i]<<setw(40)<<
     // ": " << (fr[i].accept() ? "Event Passed" : "Event Failed") << endl;
     if ( !fr[i].accept() == 1 ) continue;
     // if( filterList[i] == "Flag_HBHENoiseFilter")                     {  passed_HBHENoiseFilter = true; }    // needs to be re-run manually
     // if( filterList[i] == "Flag_HBHENoiseIsoFilter")                  { passed_HBHENoiseIsoFilter = true; }  // needs to be re-run manually
     if( filterList[i] == "Flag_CSCTightHaloFilter")                  { passed_CSCTightHaloFilter = true; }
     // if( filterList[i] == "Flag_hcalLaserEventFilter")                { passed_hcalLaserEventFilter = true; } // deprecated
     if( filterList[i] == "Flag_EcalDeadCellTriggerPrimitiveFilter")  { passed_EcalDeadCellTriggerPrimitiveFilter = true; } // under scrutiny
     if( filterList[i] == "Flag_EcalDeadCellBoundaryEnergyFilter")    { passed_EcalDeadCellBoundaryEnergyFilter = true; }   // under scrutiny
     if( filterList[i] == "Flag_goodVertices")                        { passed_goodVertices = true; }
     if( filterList[i] == "Flag_eeBadScFilter")                       { passed_eeBadScFilter = true; }
     // if( filterList[i] == "Flag_ecalLaserCorrFilter")                 { passed_ecalLaserCorrFilter = true; } // deprecated
     // if( filterList[i] == "Flag_trkPOGFilters")                       { passed_trkPOGFilters = true; } // deprecated
     // if( filterList[i] == "Flag_trkPOG_manystripclus53X")             { passed_trkPOG_manystripclus53X = true; } // deprecated
     // if( filterList[i] == "Flag_trkPOG_toomanystripclus53X")          { passed_trkPOG_toomanystripclus53X = true; } // deprecated
     // if( filterList[i] == "Flag_trkPOG_logErrorTooManyClusters")      { passed_trkPOG_logErrorTooManyClusters = true; } // deprecated
     if( filterList[i] == "Flag_METFilters")                          { passed_METFilters = true; } // be careful using this -- check documentation

   }
  
  for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){     
     jetcnt++;  
 
   if(DEBUG_){
      cout <<"jet["<<jetcnt<<"]" <<
      " | n_jets = "<<jets.size()<<
      " | pt = "<<jet->pt()<<
      " | eta = "<<fabs(jet->eta())<<
      " | NoD = "<<jet->numberOfDaughters()<<
      " | NHE = "<<jet->neutralHadronEnergyFraction()<<
      " | CHE = "<<jet->chargedHadronEnergyFraction()<<
      " | NEE = "<<jet->neutralEmEnergyFraction()<<
      " | CEE = "<<jet->chargedEmEnergyFraction()<<
      " | ChargedMulti = "<<jet->chargedMultiplicity()<<
      " | "<<endl;
    }
    // Loose Jet ID (equivalent to previous medium)
    if((
    jet->neutralHadronEnergyFraction()    <  0.99 &&  // 0.90 for tight
    jet->neutralEmEnergyFraction()        <  0.99 &&  // 0.90 for tight
    jet->numberOfDaughters()              >  1
    )                                             && 
    ((
    abs(jet->eta())                       <= 2.4  && 
    jet->chargedHadronEnergyFraction()    >  0    && 
    jet->chargedMultiplicity()            >  0    && 
    jet->chargedEmEnergyFraction()        <  0.99
    )                                             || 
    abs(jet->eta())                       >  2.4
    )                                             && 
    //abs(jet->eta())                       <= 2.6  &&
    jet->pt()                             >  20
    ) {
    pBH += jet->p4();
    ST  += jet->et();
    sumPx2  += jet->px()*jet->px();
    sumPy2  += jet->py()*jet->py();
    sumPxPy += jet->px()*jet->py();      
       
    leadingJets.push_back(jet->pt());
    leadingObjects.push_back(jet->pt());

    pJet += jet->p4();
    pObj += jet->p4(); 
         
    JetE[jetcnt]   = jet->energy();
    JetPx[jetcnt]  = jet->px();
    JetPy[jetcnt]  = jet->py();
    JetPz[jetcnt]  = jet->pz();     
    JetPt[jetcnt]  = jet->pt();
    JetEt[jetcnt]  = jet->et();
    JetEta[jetcnt] = jet->eta();      
    JetPhi[jetcnt] = jet->phi();
    ++ngoodjets;
   }//JetID
  }   
   
   for (size_t i = 0; i < electrons->size(); ++i){
   const auto e = electrons->ptrAt(i);
    //if(!e->gsfTrack()) continue;
    ++elecnt;
      if(DEBUG_){  
      cout <<"Ele["<<elecnt<<"]" <<
      " | n_el = "<< electrons->size()<<
      " | pt = "<<e->pt()<<
      " | dz = "<<e->gsfTrack()->dz( vertex_.position() )<<
      " | eta = "<<fabs(e->eta())<<
      " | dxy = "<<fabs(e->gsfTrack()->dxy(vertex_.position()))<< 
      " | hit = "<<e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<
      " | passConvVeto(a) = "<<!ConversionTools::hasMatchedConversion(*e,conversions,theBeamSpot->position())<<
      " | passMediumId = "<<(*medium_id_decisions)[e]<<
      " | "<<endl; 
      }
    // Electron Medium ID
    if(
    e->pt()           		  >  20.    &&
    fabs(e->eta())    		  <  2.5    &&
    (*medium_id_decisions)[e]      == 1      
    ){  
   
    pBH += e->p4();
    ST  += e->et();
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
    EleE[elecnt]   = e->energy();
    ElePx[elecnt]  = e->px();
    ElePy[elecnt]  = e->py();
    ElePz[elecnt]  = e->pz();     
    ElePt[elecnt]  = e->pt();
    EleEt[elecnt]  = e->et();
    EleEta[elecnt] = e->eta();      
    ElePhi[elecnt] = e->phi();
    //(*loose_id_decisions)[e];
    //(*tight_id_decisions)[e];
    ++ngoodelectrons;
     }//eleID
    }  
   

        
  for (size_t i = 0; i < photons->size(); ++i){
  const auto ph = photons->ptrAt(i);
  ++phcnt;
  if(DEBUG_){  
  cout <<"Pho["<<phcnt<<"]" <<
  " | n_ph = "<< photons->size()<<
  " | pt = "<< ph->pt()<<
  " | eta = "<< fabs(ph->superCluster()->eta())<<
  " | phi = " << ph->superCluster()->phi() <<
  " | pix_seed = "<< ph->hasPixelSeed()<<
    " | "<<endl;
  }
  if(
  ph->pt()                        >   20       &&
  abs(ph->eta())                  <   2.4      &&
  (*medium_id_decisions_ph)[ph]   ==  1
  )
  
  //If not a spike, increment # photons
  ++ngoodphotons;

  pBH += ph->p4();
  ST += ph->et();
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
  PhEt[ngoodphotons-1] = ph->et();
  PhEta[ngoodphotons-1] = ph->eta();
  //(*loose_id_decisions_ph)[ph];
  //(*tight_id_decisions_ph)[ph];      
  PhPhi[ngoodphotons-1] = ph->phi();
  
}
  for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
    ++mucnt;
    if(DEBUG_){
      cout <<"Mun["<<mucnt<<"]" <<
      " | n_mu = "<< muons.size()<<
      " | pt = "<<mu->pt()<<
      " | eta = "<<fabs(mu->eta())<<
      " | fabs(vtx_dxy) = "<<fabs(mu->globalTrack()->dxy(vertex_.position()))<<
      " | "<<endl;
    }
   if(
   mu->pt()                 >  20   &&
   mu->eta()                <  2.4  &&
   mu->isTightMuon(vertex_) == 1
   ){
    ++ngoodmuons;
    pBH += mu->p4();
    ST += mu->et();
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

    MuE[ngoodmuons-1]   = mu->energy();
    MuPx[ngoodmuons-1]  = mu->px();
    MuPy[ngoodmuons-1]  = mu->py();
    MuPz[ngoodmuons-1]  = mu->pz(); 
    MuPt[ngoodmuons-1]  = mu->pt();
    MuEt[ngoodmuons-1]  = mu->et();
    MuEta[ngoodmuons-1] = mu->eta();      
    MuPhi[ngoodmuons-1] = mu->phi(); 
    }//muonID
  }
   
   for (int i=0;i<25;++i) {
     if (i>=ngoodjets) {
       JetE[i]=0.;
       JetPx[i]=0.;
       JetPy[i]=0.;
       JetPz[i]=0.;
       JetPt[i]=0.;
       JetEt[i]=0.;
       JetEta[i]=99.;
       JetPhi[i]=99.;
     }
     if (i>=ngoodelectrons) {
       EleE[i]=0.;
       ElePx[i]=0.;
       ElePy[i]=0.;
       ElePz[i]=0.;
       ElePt[i]=0.;
       EleEt[i]=0.;
       EleEta[i]=99.;
       ElePhi[i]=99.;
     }
     if (i>=ngoodphotons) {
       PhE[i]=0.;
       PhPx[i]=0.;
       PhPy[i]=0.;
       PhPz[i]=0.;
       PhPt[i]=0.;
       PhEt[i]=0.;
       PhEta[i]=99.;
       PhPhi[i]=99.;
       
     }       
     if (i>=ngoodmuons) {
       MuE[i]=0.;
       MuPx[i]=0.;
       MuPy[i]=0.;
       MuPz[i]=0.;
       MuPt[i]=0.;
       MuEt[i]=0.;
       MuEta[i]=99.;
       MuPhi[i]=99.;
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
   ResJetEt = pJet.Et();

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

  tree->Branch("NJets",&NJets,"NJets/I");
  tree->Branch("JetE",&JetE,"JetE[NJets]/F");
  tree->Branch("JetPx",&JetPx,"JetPx[NJets]/F");
  tree->Branch("JetPy",&JetPy,"JetPy[NJets]/F");
  tree->Branch("JetPz",&JetPz,"JetPz[NJets]/F");
  tree->Branch("JetPt",&JetPt,"JetPt[NJets]/F");
  tree->Branch("JetEt",&JetEt,"JetEt[NJets]/F");
  tree->Branch("JetEta",&JetEta,"JetEta[NJets]/F");
  tree->Branch("JetPhi",&JetPhi,"JetPhi[NJets]/F");

  tree->Branch("EleE" ,&EleE, "EleE[25]/F");
  tree->Branch("ElePx",&ElePx,"ElePx[25]/F");
  tree->Branch("ElePy",&ElePy,"ElePy[25]/F");
  tree->Branch("ElePz",&ElePz,"ElePz[25]/F");
  tree->Branch("ElePt",&ElePt,"ElePt[25]/F");
  tree->Branch("EleEt",&EleEt,"EleEt[25]/F");
  tree->Branch("EleEta",&EleEta,"EleEta[25]/F");
  tree->Branch("ElePhi",&ElePhi,"ElePhi[25]/F");
 
  tree->Branch("PhE",&PhE,"PhE[25]/F");
  tree->Branch("PhPx",&PhPx,"PhPx[25]/F");
  tree->Branch("PhPy",&PhPy,"PhPy[25]/F");
  tree->Branch("PhPz",&PhPz,"PhPz[25]/F");
  tree->Branch("PhPt",&PhPt,"PhPt[25]/F");
  tree->Branch("PhEt",&PhEt,"PhEt[25]/F");
  tree->Branch("PhEta",&PhEta,"PhEta[25]/F");
  tree->Branch("PhPhi",&PhPhi,"PhPhi[25]/F");

  tree->Branch("MuE",&MuE,"MuE[25]");
  tree->Branch("MuPx",&MuPx,"MuPx[25]");
  tree->Branch("MuPy",&MuPy,"MuPy[25]");
  tree->Branch("MuPz",&MuPz,"MuPz[25]");
  tree->Branch("MuPt",&MuPt,"MuPt[25]");
  tree->Branch("MuEt",&MuEt,"MuEt[25]");
  tree->Branch("MuEta",&MuEta,"MuEta[25]");
  tree->Branch("MuPhi",&MuPhi,"MuPhi[25]");
    
  tree->Branch("ST",&ST,"ST/F");
  tree->Branch("mBH",&mBH,"mBH/F");
  tree->Branch("Met",&Met,"Met/F");
  tree->Branch("MetE",&MetE,"MetE/F");
  tree->Branch("MetPx",&MetPx,"MetPx/F");
  tree->Branch("MetPy",&MetPy,"MetPy/F");
  tree->Branch("MetPz",&MetPz,"MetPz/F");
  tree->Branch("MetPt",&MetPt,"MetPt/F");
  tree->Branch("MetEt",&MetEt,"MetEt/F");
  tree->Branch("MetPhi",&MetPhi,"MetPhi/F");   
  tree->Branch("Sphericity", &Sphericity, "Sphericity/F");
  tree->Branch("Jet",JetArr,"JetArr[4]");   
  tree->Branch("Ele",EleArr,"EleArr[4]");
  tree->Branch("Mu",MuArr,"MuArr[4]");
  tree->Branch("Ph",PhArr,"PhArr[4]");   
  tree->Branch("NPV",&NPV,"NPV/I");   
  tree->Branch("NTracks",&NTracks,"NTracks/I");   
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
  tree->Branch("ResJetEt",&ResJetEt,"ResJetEt/F");
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
  
  //tree->Branch("firedHLT_PFJet60_v2",&firedHLT_PFJet60_v2,"firedHLT_PFJet60_v2/B");
  //tree->Branch("firedHLT_PFJet140_v2",&firedHLT_PFJet140_v2,"firedHLT_PFJet140_v2/B");
  //tree->Branch("firedHLT_PFJet450_v2",&firedHLT_PFJet450_v2,"firedHLT_PFJet450_v2/B");
  //tree->Branch("firedHLT_PFHT300_v1",&firedHLT_PFHT300_v1,"firedHLT_PFHT300_v1/B");
  //tree->Branch("firedHLT_PFHT400_v1",&firedHLT_PFHT400_v1,"firedHLT_PFHT400_v1/B");
  tree->Branch("firedHLT_PFHT475_v2",&firedHLT_PFHT475_v2,"firedHLT_PFHT475_v2/B");
  //tree->Branch("firedHLT_PFHT600_v2",&firedHLT_PFHT600_v2,"firedHLT_PFHT600_v2/B");
  //tree->Branch("firedHLT_PFHT650_v2",&firedHLT_PFHT650_v2,"firedHLT_PFHT650_v2/B");
  tree->Branch("firedHLT_PFHT800_v2",&firedHLT_PFHT800_v2,"firedHLT_PFHT800_v2/B");

  //tree->Branch("passed_HBHENoiseFilter", &passed_HBHENoiseFilter, "passed_HBHENoiseFilter/B"); 
  //tree->Branch("passed_HBHENoiseIsoFilter", &passed_HBHENoiseIsoFilter, "passed_HBHENoiseIsoFilter/B"); 
  tree->Branch("passed_CSCTightHaloFilter",&passed_CSCTightHaloFilter,"passed_CSCTightHaloFilter/B");
  //tree->Branch("passed_hcalLaserEventFilter", &passed_hcalLaserEventFilter, "passed_hcalLaserEventFilter/B"); 
  tree->Branch("passed_EcalDeadCellTriggerPrimitiveFilter", &passed_EcalDeadCellTriggerPrimitiveFilter, "passed_EcalDeadCellTriggerPrimitiveFilter/B"); 
  tree->Branch("passed_EcalDeadCellBoundaryEnergyFilter", &passed_EcalDeadCellBoundaryEnergyFilter, "passed_EcalDeadCellBoundaryEnergyFilter/B"); 
  tree->Branch("passed_goodVertices", &passed_goodVertices, "passed_goodVertices/B"); 
  tree->Branch("passed_eeBadScFilter", &passed_eeBadScFilter, "passed_eeBadScFilter/B"); 
  //tree->Branch("passed_ecalLaserCorrFilter", &passed_ecalLaserCorrFilter, "passed_ecalLaserCorrFilter/B"); 
  //tree->Branch("passed_trkPOGFilters", &passed_trkPOGFilters, "passed_trkPOGFilters/B"); 
  //tree->Branch("passed_trkPOG_manystripclus53X", &passed_trkPOG_manystripclus53X, "passed_trkPOG_manystripclus53X/B"); 
  //tree->Branch("passed_trkPOG_toomanystripclus53X", &passed_trkPOG_toomanystripclus53X, "passed_trkPOG_toomanystripclus53X/B"); 
  //tree->Branch("passed_trkPOG_logErrorTooManyClusters", &passed_trkPOG_logErrorTooManyClusters, "passed_trkPOG_logErrorTooManyClusters/B"); 
  tree->Branch("passed_METFilters", &passed_METFilters, "passed_METFilters/B"); 
    
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

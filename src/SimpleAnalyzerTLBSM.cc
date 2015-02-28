#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TH1D.h"
#include "TH2D.h"
#include <boost/foreach.hpp>

//
// class declaration
//
using namespace edm;
using namespace std;

class SimpleAnalyzerTLBSM : public edm::EDAnalyzer {
public:
  explicit SimpleAnalyzerTLBSM(const edm::ParameterSet&);
  ~SimpleAnalyzerTLBSM();
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  edm::InputTag photonSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag metSrc_;
  double metThreshold_;
  bool DEBUG_;
  
  TH1D* hN_;
  TH1D* hST_;
  TH2D* hN_vs_ST_;
  TH2D* hNjet_vs_ST_; //for multijet events, no lepton/photon
};

SimpleAnalyzerTLBSM::SimpleAnalyzerTLBSM(const edm::ParameterSet& iConfig):
  muonSrc_(iConfig.getParameter<edm::InputTag>("muonSrc")),
  electronSrc_(iConfig.getParameter<edm::InputTag>("electronSrc")),
  photonSrc_(iConfig.getParameter<edm::InputTag>("photonSrc")),
  jetSrc_(iConfig.getParameter<edm::InputTag>("jetSrc")),
  metSrc_(iConfig.getParameter<edm::InputTag>("metSrc")),
  metThreshold_(iConfig.getParameter<double>("metThreshold")),
  DEBUG_(iConfig.getParameter<bool>("DEBUG"))
  
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  hN_ = fs->make<TH1D>("multiplicity", "multiplicity", 50, 0, 50);
  hST_ = fs->make<TH1D>("ST", "ST", 1000, 0, 10000);
  hN_vs_ST_ = fs->make<TH2D>("N_vs_ST", "N_vs_ST",
			     1000, 0, 10000,
			     50, 0, 50);
  hNjet_vs_ST_ = fs->make<TH2D>("Njet_vs_ST", "Njet_vs_ST",
				1000, 0, 10000,
				50, 0, 50);
  
}

SimpleAnalyzerTLBSM::~SimpleAnalyzerTLBSM()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SimpleAnalyzerTLBSM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  double ST = 0.0;
  int N = 0;
  
  //------ Primary Vertices     
  edm::Handle< reco::VertexCollection > PVCollection; 
  iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), PVCollection );  
  reco::Vertex primaryVertex;   
  primaryVertex = PVCollection->at(0); // No need to check if PV is valid, since it is required in the config
  
  // get the beamspot
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  if (!iEvent.getByLabel(edm::InputTag("offlineBeamSpot"), beamSpotHandle)) { cout<<"No beam spot "<<endl; return;} 
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_, jets);
  
  //Muons 
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_, muons);
  
  //Jets
  int njets = 0;
  BOOST_FOREACH(const pat::Jet& jet, *jets){
    if (jet.pt() > metThreshold_ && fabs(jet.eta()) < 2.6) { 
      ST += jet.pt();
      ++njets;
    }
  }
  N += njets;
  
  //Electrons
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);
  
  int nelectrons = 0;
  BOOST_FOREACH(const pat::Electron& electron, *electrons){
    if (electron.pt() > metThreshold_ && fabs(electron.eta()) < 2.4) {
      
      const float ElectronD0Cut_ = 0.04, EtaThr_e = 2.5, ElectronETSCThr_ = 15., RemuThr = 0.1, ElectronVertexMatchThr_ = 1., RelIso_e = 0.2;
      float D0_e = 99., Chi2_e = 99., ET_SC = -1.;
      
      bool hadId((int)electron.electronID("eidTight") & 0x1);
      bool isNotConv((int)electron.electronID("eidTight") & 0x4); 
      bool isGsfElectron = true;
      if (!electron.gsfTrack()) isGsfElectron = false;
      int nlosthits   = electron.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      double convDist = electron.convDist();
      double convDcot = electron.convDcot();
      D0_e      = electron.gsfTrack()->dxy(beamSpotHandle->position());
      ET_SC = electron.superCluster()->energy() * TMath::Sin ( electron.superCluster()->position().theta() );
      double Iso = (electron.chargedHadronIso()+electron.neutralHadronIso()+electron.photonIso())/electron.pt();
      
      if (DEBUG_) {
	cout<<"isGsf "<<isGsfElectron<<", hadID "<<hadId<<", not converted? "<<isNotConv<<", ecalseed? "<<electron.ecalDrivenSeed()<<endl;    
	cout<<"N lost hits "<<nlosthits<<", cut < 2"<<endl;
	cout<<"convDcot "<<fabs(convDcot)<<", convDist "<<fabs(convDist)<<", cut > 0.02"<<endl;
	cout<<"El D0 "<<fabs(D0_e)<<", cut < "<< ElectronD0Cut_<<endl;
	cout<<"ET_SC "<<ET_SC<<", cut > "<<ElectronETSCThr_<<endl;
	cout<<"El relative iso "<<Iso<<", cut < "<<RelIso_e<<endl;    
	cout<<"El vertex x/y/z "<<electron.vx()<<" "<<electron.vy()<<" "<<electron.vz()<<", PV x/y/z "<<primaryVertex.x()<<" "<<primaryVertex.y()<<" "<<primaryVertex.z()<<", beamspot x/y/z "<<beamSpotHandle->position().x()<<" "<<beamSpotHandle->position().y()<<" "<<beamSpotHandle->position().z()<<endl;
      }
      
      if(!isGsfElectron || !hadId || !isNotConv)                  		continue; 	      
      if(nlosthits >= 2 || (fabs(convDcot) < 0.02 && fabs(convDist) < 0.02)) 	continue;	
      if(electron.ecalDrivenSeed() != 1)                     			continue;
      if(fabs(D0_e)       >= ElectronD0Cut_)         				continue; 
      if(fabs(electron.eta())   >= EtaThr_e)               				continue;
      if(ET_SC            <= ElectronETSCThr_)       				continue;
      if (Iso             > RelIso_e )               				continue;
      
      float dRmue = 0.;    
      bool SharedCone = false; //DR computed between localElectrons candidates and all Global or Tracker Muon
      BOOST_FOREACH(const pat::Muon& muon, *muons){
	if ((muon.isTrackerMuon()) || (muon.isGlobalMuon())){
	  dRmue = deltaR(electron.eta(),electron.phi(),muon.eta(),muon.phi());
	  if (dRmue < RemuThr){
	    cout<<"Electron overlaps with global or tracker muon"<<endl;
	    SharedCone = true;
	    break;
	  }
	}
      }
      
      if ( SharedCone )                      continue;
      if ( fabs( electron.vz() - primaryVertex.z() )  > ElectronVertexMatchThr_ ) continue;
      
      float dRej = 0.;
      int overlap_ej = 0;    
      BOOST_FOREACH(const pat::Jet& jet, *jets){
        dRej = deltaR(jet.eta(),jet.phi(),electron.eta(),electron.phi());
        if (DEBUG_) cout<<"Run "<<iEvent.id().run()<<", Event"<<iEvent.id().event()<<" , Lumi "<<iEvent.id().luminosityBlock()<<", dR electron-jet "<<dRej<<", electron pt "<<electron.pt()<<", jet pt "<<jet.pt()<<", electron eta "<<electron.eta()<<", jet eta "<<jet.eta()<<", electron phi "<<electron.phi()<<", jet phi "<<jet.phi()<<endl;
        if (dRej < 0.3) ++overlap_ej;
      }
      if (DEBUG_) cout << "Overlap e-jet "<<overlap_ej<<endl;
      //if (overlap_ej > 0) continue; //this does not apply
      
      ST += electron.pt();
      ++nelectrons;
    }
    
  }
  N += nelectrons;
  
  int nphot = 0.;
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByLabel(photonSrc_, photons);
  
  //cout <<"Photons size "<<photons.size()<<endl;
  BOOST_FOREACH(const pat::Photon& photon, *photons){
    if (photon.pt() < metThreshold_ && fabs(photon.eta()) < 2.4) continue;
    
    float dRphj = 0.;
    int overlap_phj = 0;    
    BOOST_FOREACH(const pat::Jet& jet, *jets){
      dRphj = deltaR(jet.eta(),jet.phi(),photon.eta(),photon.phi());
      if (DEBUG_) cout<<"Run "<<iEvent.id().run()<<", Event"<<iEvent.id().event()<<" , Lumi "<<iEvent.id().luminosityBlock()<<", dR photon-jet "<<dRphj<<", photon pt "<<photon.pt()<<", jet pt "<<jet.pt()<<", photon eta "<<photon.eta()<<", jet eta "<<jet.eta()<<", photon phi "<<photon.phi()<<", jet phi "<<jet.phi()<<endl;
      if (dRphj < 0.3) ++overlap_phj;
    }
    if (DEBUG_) cout << "overlap ph jet "<<overlap_phj<<endl;
    if (overlap_phj > 0) continue;
    
    ST += photon.pt();
    //cout<<"Photons pt "<<photon.pt()<<endl;
    ++nphot;      
  }
  N += nphot; //Do not use photons.size() for Nphotons! 
  if (DEBUG_) cout <<"NPhot "<<nphot<<endl;
  
  //Muons  
  int nmuons = 0;
  BOOST_FOREACH(const pat::Muon& muon, *muons){
    if (muon.pt() > metThreshold_ && fabs(muon.eta()) < 2.1) { 
      
      if (!(muon.isTrackerMuon()) || !(muon.isGlobalMuon())) continue;
      
      const float MuonRelIso = 0.02, MuonD0Cut = 0.02, MuonVertexMatchThr = 1., MuonNofValidHits = 0, MuonNofValidTrHits = 10., MuonNormChi2 = 10.;
      float D0 = 99., Chi2 = 99., D0Inner = 99., D0Standard = 99., RelIso03PF = 99.;
      int NValidHits = -1, NTrValidHits = -1;
      
      // Get the tracker track from the muon       
      const reco::TrackRef globalTrack = muon.globalTrack();
      if (globalTrack.isNonnull()) {
	D0   = globalTrack->dxy(beamSpotHandle->position());
	Chi2 = globalTrack->normalizedChi2();
      }
      if (muon.innerTrack().isNonnull()) D0Inner = muon.innerTrack()->dxy (beamSpotHandle->position());
      if (muon.isGlobalMuon())        NValidHits = muon.globalTrack()->hitPattern().numberOfValidMuonHits();
      NTrValidHits = muon.innerTrack()->numberOfValidHits();
      D0Standard = muon.dB(pat::Muon::BS2D); //Much more elegant way to find D0. Compare it to D0Inner!
      
      float muvz = muon.vz();
      float pvx = primaryVertex.x();
      float pvy = primaryVertex.y();
      float pvz = primaryVertex.z();
      
      //Isolation: 
      float PATNeutralHadronIso =  muon.neutralHadronIso();
      float PATChargedHadronIso =  muon.chargedHadronIso();
      float PATPhotonIso        =  muon.photonIso();
      float PATTrackIso         =  muon.trackIso(); 
      RelIso03PF                = ((PATNeutralHadronIso+PATChargedHadronIso+PATPhotonIso)/muon.pt());
      
      if (DEBUG_) {
	std::cout << "Chi2: " << Chi2 << " Treshold: " << MuonNormChi2 <<  std::endl;
	std::cout << "NTrValidHits: " << NTrValidHits << " Treshold: " << MuonNofValidTrHits << std::endl;
	std::cout << "NValidHits: " << NValidHits << " Treshold: " << MuonNofValidHits << std::endl;
	std::cout << "D0Inner: " << fabs(D0Inner) << " Treshold: " << MuonD0Cut << " D0 standard "<<fabs(D0Standard)<<std::endl;
	std::cout << "Eta: " << fabs( muon.eta() ) << std::endl;
	std::cout << "Pt: " << muon.pt() << " Treshold: " << metThreshold_ << std::endl;
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
      BOOST_FOREACH(const pat::Jet& jet, *jets){
        dRmuj = deltaR(muon.eta(),muon.phi(),jet.eta(),jet.phi());
        if (DEBUG_) cout<<"Run "<<iEvent.id().run()<<", Event"<<iEvent.id().event()<<" , Lumi "<<iEvent.id().luminosityBlock()<<", dR muon-jet "<<dRmuj<<", muon pt "<<muon.pt()<<", jet pt "<<jet.pt()<<", muon eta "<<muon.eta()<<", jet eta "<<jet.eta()<<", muon phi "<<muon.phi()<<", jet phi "<<jet.phi()<<" muon d0 "<<fabs(D0Standard)<<endl;
        if (dRmuj < 0.3) ++overlap_muj;
      }
      if (DEBUG_) cout<<"Overlap muon-jet "<<overlap_muj<<endl;
      //Need to think, if we want to remove such muons...     
      
      ST += muon.pt();
      ++nmuons;
    }
  }
  N += nmuons;
   
  //cout<<"ST before MET "<<ST<<endl;
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_, mets);
  const pat::MET& met = mets->at(0);
  
  //double met_pt = muons->size() > 0
  //   ? met.pt() : met.uncorrectedPt(pat::MET::uncorrMUON);
  if (met.pt() > metThreshold_)
    ST += met.pt();


  //MET tests
  //cout <<"MET Before "<<met.pt()<<" muons size "<<muons->size()<<endl;
  //double met_pt = muons->size() > 0
  //    ? met.pt() : met.uncorrectedPt(pat::MET::uncorrMUON);
  //cout <<"MET After  "<<met_pt<<endl;

  
  //cout <<"Njets (Simple) "<<njets<<" , ST "<<ST<<" , mult "<<N<<" , jets->size() "<<jets->size()<<endl;

  //cout<<"MET "<<met.pt()<<" ST "<<ST<<endl;
  hN_->Fill(N);
  hST_->Fill(ST);
  hN_vs_ST_->Fill(ST, N);
  if (N == njets)
    hNjet_vs_ST_->Fill(ST, N);
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
SimpleAnalyzerTLBSM::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleAnalyzerTLBSM::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(SimpleAnalyzerTLBSM);

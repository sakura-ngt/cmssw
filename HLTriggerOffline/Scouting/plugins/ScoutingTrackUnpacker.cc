// system includes
#include <memory>

// user includes
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "fastjet/contrib/SoftKiller.hh"

/*
 * HLTScoutingUnpackProducer unpacks Run3Scouting data formats to their reco format counterparts
 * and cross link collections (currently PFCandidate - Track)
*/

class HLTScoutingUnpackProducer : public edm::stream::EDProducer<> {
public:
  using Run3ScoutingPFJetCollection = std::vector<Run3ScoutingPFJet>;
  using Run3ScoutingParticleCollection = std::vector<Run3ScoutingParticle>;
  using Run3ScoutingTrackCollection = std::vector<Run3ScoutingTrack>;
  using Run3ScoutingVertexCollection = std::vector<Run3ScoutingVertex>;
  using PFCandCollection = std::vector<pat::PackedCandidate>;

  template <typename T>
  using RefCollection = std::vector<edm::Ref<std::vector<T>>>;
  template <typename T>
  using RefMap = edm::ValueMap<edm::Ref<std::vector<T>>>;

  explicit HLTScoutingUnpackProducer(edm::ParameterSet const& params);
  ~HLTScoutingUnpackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup) override;

  // private helper functions
  reco::PFJet createPFJet(Run3ScoutingPFJet const& scoutingPFJet);
  reco::Vertex createVertex(Run3ScoutingVertex const& scoutingVertex);
  reco::Track createTrack(Run3ScoutingTrack const& scoutingTrack);
  reco::PFCandidate createPFCandidate(Run3ScoutingParticle const& scoutingPFCandidate);
  bool hasPFTrackDetails(Run3ScoutingParticle const& scoutingPFCandidate);
  int findCompatibleScoutingTrack(edm::Handle<Run3ScoutingParticleCollection> const& scoutingParticleCollection,
                                  Run3ScoutingTrack const& scoutingTrack);
  void buildHitPattern(Run3ScoutingParticle const& scoutingPFCandidate,
                       Run3ScoutingTrack const& scoutingTrack,
                       reco::Track& recoTrack);
  std::unique_ptr<reco::Track> createPFTrack(
      Run3ScoutingParticle const& scoutingPFCandidate,
      edm::Handle<Run3ScoutingVertexCollection> const& scoutingPrimaryVertex_collection_handle);

  // register products for reco object and corresponding Value to Ref to original scouting objects
  template <typename RecoObjectType, typename ScoutingObjectType>
  void produceWithRef(std::string const& name);

  // put reco object together with ValueMap to Ref to original scouting object
  template <typename RecoObjectType, typename ScoutingObjectType>
  void putWithRef(edm::Event& iEvent,
                  std::string const& name,
                  std::unique_ptr<std::vector<RecoObjectType>>& recoObject_collection_ptr,
                  std::unique_ptr<RefCollection<ScoutingObjectType>>& scoutingObjectRef_collection_ptr);

  edm::EDGetTokenT<Run3ScoutingTrackCollection> scoutingTrack_collection_token_;
  edm::EDGetTokenT<Run3ScoutingVertexCollection> scoutingPrimaryVertex_collection_token_;
  edm::EDGetTokenT<Run3ScoutingParticleCollection> scoutingParticle_collection_token_;
  edm::EDGetTokenT<PFCandCollection> pfCand_collection_token_;
  edm::EDGetTokenT<PFCandCollection> lostTrack_collection_token_;
  bool isScouting_;
  bool produce_PFCandidate_;
  bool produce_PFCandidateMatchTrack_;
  bool produce_PFCHSCandidate_;  // CHS = charged hadron subtraction
  bool produce_PFSKCandidate_;   // SK = soft killer

  inline static const std::string REF_TO_SCOUTING_LABEL_SUFFIX_ = "-RefToOriginal";
};

// constructor
HLTScoutingUnpackProducer::HLTScoutingUnpackProducer(edm::ParameterSet const& params)
    : scoutingTrack_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingTrack"))),
      scoutingPrimaryVertex_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingPrimaryVertex"))),
      scoutingParticle_collection_token_(consumes(params.getParameter<edm::InputTag>("scoutingParticle"))),
      pfCand_collection_token_(consumes(params.getParameter<edm::InputTag>("pfCand"))),
      lostTrack_collection_token_(consumes(params.getParameter<edm::InputTag>("lostTrack"))),
      isScouting_(params.getParameter<bool>("isScouting")),
      produce_PFCHSCandidate_(params.getParameter<bool>("producePFCHSCandidate")) {
  //produce_PFSKCandidate_(params.getParameter<bool>("producePFSKCandidate"))

  if (isScouting_) {
    produceWithRef<reco::Track, Run3ScoutingTrack>("Track");
    produceWithRef<reco::Vertex, Run3ScoutingVertex>("PrimaryVertex");
  } else {
    produces<std::vector<reco::Track>>("Track");
  }

  if (produce_PFCandidate_) {
    produceWithRef<reco::PFCandidate, Run3ScoutingParticle>("PFCandidate");
  }
  if (produce_PFCandidateMatchTrack_) {
    produceWithRef<reco::PFCandidate, Run3ScoutingParticle>("PFCandidateMatchTrack");
  }
  if (produce_PFCHSCandidate_) {
    produceWithRef<reco::PFCandidate, Run3ScoutingParticle>("PFCHSCandidate");
  }
  /*if (produce_PFSKCandidate_){
        produces<reco::PFCandidateCollection>("PFSKCandidate");
    }*/
  if (produce_PFCandidate_ || produce_PFCandidateMatchTrack_ || produce_PFCHSCandidate_ || produce_PFSKCandidate_) {
    produces<reco::TrackCollection>("PFTrack");  // reco::track collection from best track parameters of PF
  }
}

template <typename RecoObjectType, typename ScoutingObjectType>
void HLTScoutingUnpackProducer::produceWithRef(std::string const& name) {
  produces<std::vector<RecoObjectType>>(name);
  produces<RefMap<ScoutingObjectType>>(name + REF_TO_SCOUTING_LABEL_SUFFIX_);
}

void HLTScoutingUnpackProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup) {
  // produce reco::Vertex
  edm::Handle<Run3ScoutingVertexCollection> scoutingPrimaryVertex_collection_handle =
      iEvent.getHandle(scoutingPrimaryVertex_collection_token_);
  auto recoPrimaryVertex_collection_ptr = std::make_unique<reco::VertexCollection>();
  auto scoutingPrimaryVertexRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingVertex>>();
  if (scoutingPrimaryVertex_collection_handle.isValid()) {
    //std::cout<<"number of scouting vertices: "<<scoutingPrimaryVertex_collection_handle->size()<<std::endl;
    for (size_t scoutingPrimaryVertex_index = 0;
         scoutingPrimaryVertex_index < scoutingPrimaryVertex_collection_handle->size();
         scoutingPrimaryVertex_index++) {
      auto& scoutingPrimaryVertex = scoutingPrimaryVertex_collection_handle->at(scoutingPrimaryVertex_index);
      recoPrimaryVertex_collection_ptr->push_back(createVertex(scoutingPrimaryVertex));
      scoutingPrimaryVertexRef_collection_ptr->push_back(
          edm::Ref<Run3ScoutingVertexCollection>(scoutingPrimaryVertex_collection_handle, scoutingPrimaryVertex_index));
    }
  }

  edm::Handle<Run3ScoutingParticleCollection> scoutingParticle_collection_handle =
      iEvent.getHandle(scoutingParticle_collection_token_);
  if (scoutingParticle_collection_handle.isValid()) {
    //std::cout<<"scouting particle size: "<<scoutingParticle_collection_handle->size()<<std::endl;
  }

  // produce reco::Track
  edm::Handle<Run3ScoutingTrackCollection> scoutingTrack_collection_handle =
      iEvent.getHandle(scoutingTrack_collection_token_);
  auto recoTrack_collection_ptr = std::make_unique<reco::TrackCollection>();
  auto scoutingTrackRef_collection_ptr = std::make_unique<RefCollection<Run3ScoutingTrack>>();
  if (scoutingTrack_collection_handle.isValid()) {
    //std::cout<<"scouting track size: "<<scoutingTrack_collection_handle->size()<<std::endl;
    for (size_t scoutingTrack_index = 0; scoutingTrack_index < scoutingTrack_collection_handle->size();
         scoutingTrack_index++) {
      auto& scoutingTrack = scoutingTrack_collection_handle->at(scoutingTrack_index);
      int match_index = findCompatibleScoutingTrack(scoutingParticle_collection_handle, scoutingTrack);
      //std::cout<<"matched particle index: "<<match_index<<std::endl;
      reco::Track track = createTrack(scoutingTrack);
      if (match_index != -1) {
        buildHitPattern(scoutingParticle_collection_handle->at(match_index), scoutingTrack, track);
        //std::cout<<"pf lost inner hits: "<<(int) scoutingParticle_collection_handle->at(match_index).lostInnerHits()<<" missing inner hits: "<<track.missingInnerHits()<<" pf pdgid: "<<scoutingParticle_collection_handle->at(match_index).pdgId()<<" number of lost hits: "<<track.missingInnerHits()<<std::endl;
        //std::cout<<"pf track quality: " << (int) scoutingParticle_collection_handle->at(match_index).quality()<<std::endl;
      }
      recoTrack_collection_ptr->push_back(track);
      scoutingTrackRef_collection_ptr->push_back(
          edm::Ref<Run3ScoutingTrackCollection>(scoutingTrack_collection_handle, scoutingTrack_index));
    }
  }

  edm::Handle<PFCandCollection> pfCand_collection_handle = iEvent.getHandle(pfCand_collection_token_);
  auto pfRecoTrack_collection_ptr = std::make_unique<reco::TrackCollection>();
  if (pfCand_collection_handle.isValid()) {
    for (size_t pfCand_index = 0; pfCand_index < pfCand_collection_handle->size(); pfCand_index++) {
      auto& pfCand = pfCand_collection_handle->at(pfCand_index);
      const reco::Track* pfTrack = pfCand.bestTrack();
      if (pfTrack != nullptr) {
        pfRecoTrack_collection_ptr->push_back(*pfTrack);
      }
    }
  }

  edm::Handle<PFCandCollection> lostTrack_collection_handle = iEvent.getHandle(lostTrack_collection_token_);
  if (lostTrack_collection_handle.isValid()) {
    for (size_t lostTrack_index = 0; lostTrack_index < lostTrack_collection_handle->size(); lostTrack_index++) {
      auto& lostTrackPF = lostTrack_collection_handle->at(lostTrack_index);
      const reco::Track* lostTrack = lostTrackPF.bestTrack();
      if (lostTrack != nullptr) {
        pfRecoTrack_collection_ptr->push_back(*lostTrack);
        //std::cout<<"lost chi2: "<<lostTrack->chi2()<<" ndof: "<<lostTrack->ndof()<<" reduced chi2: "<<float(lostTrack->chi2())/lostTrack->ndof()<<std::endl;
      }
    }
  }

  // build fake HitPattern for ScoutingTracks which are not matched with ScoutingPFCandidate
  // this might be needed if perform vertexing with all tracks

  // put products in Event
  if (isScouting_) {
    putWithRef<reco::Vertex, Run3ScoutingVertex>(
        iEvent, "PrimaryVertex", recoPrimaryVertex_collection_ptr, scoutingPrimaryVertexRef_collection_ptr);
    putWithRef<reco::Track, Run3ScoutingTrack>(
        iEvent, "Track", recoTrack_collection_ptr, scoutingTrackRef_collection_ptr);
  } else {
    iEvent.put(std::move(pfRecoTrack_collection_ptr), "Track");
  }
}

reco::Vertex HLTScoutingUnpackProducer::createVertex(Run3ScoutingVertex const& scoutingVertex) {
  // fill point coordinate
  //std::cout<<"scouting vertex position:"<<scoutingVertex.x()<<" "<<scoutingVertex.y()<<" "<<scoutingVertex.z()<<std::endl;
  reco::Vertex::Point point(scoutingVertex.x(), scoutingVertex.y(), scoutingVertex.z());

  // fill error
  std::vector<float> error_vec(6);
  error_vec[0] = scoutingVertex.xError() * scoutingVertex.xError();  // cov(0, 0)
  error_vec[1] = 0;                                                  // cov(0, 1)
  error_vec[2] = 0;                                                  // cov(0, 2)
  error_vec[3] = scoutingVertex.yError() * scoutingVertex.yError();  // cov(1, 1)
  error_vec[4] = 0;                                                  // cov(1, 2)
  error_vec[5] = scoutingVertex.zError() * scoutingVertex.zError();  // cov(2, 2)

  // off-diagonal errors are added in the begining of 2024
  // see https://github.com/cms-sw/cmssw/pull/43758
  /*try { 
      error_vec[1] = scoutingVertex.xyCov();
      error_vec[2] = scoutingVertex.xzCov();
      error_vec[4] = scoutingVertex.yzCov();
      }
    catch (...) { 
    }*/
  //Commented out the above to use for 2023

  reco::Vertex::Error error(error_vec.begin(), error_vec.end());

  return scoutingVertex.isValidVtx()
             ? reco::Vertex(point, error, scoutingVertex.chi2(), scoutingVertex.ndof(), scoutingVertex.tracksSize())
             : reco::Vertex(point, error);
}

reco::Track HLTScoutingUnpackProducer::createTrack(Run3ScoutingTrack const& scoutingTrack) {
  float chi2 = scoutingTrack.tk_chi2();
  //std::cout<<"scout chi2: "<<chi2<<std::endl;
  float ndof = scoutingTrack.tk_ndof();
  reco::TrackBase::Point referencePoint(scoutingTrack.tk_vx(), scoutingTrack.tk_vy(), scoutingTrack.tk_vz());

  float px = scoutingTrack.tk_pt() * cos(scoutingTrack.tk_phi());
  float py = scoutingTrack.tk_pt() * sin(scoutingTrack.tk_phi());
  float pz = scoutingTrack.tk_pt() * sinh(scoutingTrack.tk_eta());
  reco::TrackBase::Vector momentum(px, py, pz);

  int charge = scoutingTrack.tk_charge();

  std::vector<float> cov_vec(15);                                                  // 5*(5+1)/2 = 15
  cov_vec[0] = scoutingTrack.tk_qoverp_Error() * scoutingTrack.tk_qoverp_Error();  // cov(0, 0)
  cov_vec[1] = scoutingTrack.tk_qoverp_lambda_cov();                               // cov(0, 1)
  cov_vec[3] = scoutingTrack.tk_qoverp_phi_cov();                                  // cov(0, 2)
  cov_vec[6] = scoutingTrack.tk_qoverp_dxy_cov();                                  // cov(0, 3)
  cov_vec[10] = scoutingTrack.tk_qoverp_dsz_cov();                                 // cov(0, 4)
  cov_vec[2] = scoutingTrack.tk_lambda_Error() * scoutingTrack.tk_lambda_Error();  // cov(1, 1)
  cov_vec[4] = scoutingTrack.tk_lambda_phi_cov();                                  // cov(1, 2)
  cov_vec[7] = scoutingTrack.tk_lambda_dxy_cov();                                  // cov(1, 3)
  cov_vec[11] = scoutingTrack.tk_lambda_dsz_cov();                                 // cov(1, 4)
  cov_vec[5] = scoutingTrack.tk_phi_Error() * scoutingTrack.tk_phi_Error();        // cov(2, 2)
  cov_vec[8] = scoutingTrack.tk_phi_dxy_cov();                                     // cov(2, 3)
  cov_vec[12] = scoutingTrack.tk_phi_dsz_cov();                                    // cov(2, 4)
  cov_vec[9] = scoutingTrack.tk_dxy_Error() * scoutingTrack.tk_dxy_Error();        // cov(3, 3)
  cov_vec[13] = scoutingTrack.tk_dxy_dsz_cov();                                    // cov(3, 4)
  cov_vec[14] = scoutingTrack.tk_dsz_Error() * scoutingTrack.tk_dsz_Error();       // cov(4, 4)
  reco::TrackBase::CovarianceMatrix cov(cov_vec.begin(), cov_vec.end());

  //std::cout<<"scouting track vertex index: "<<scoutingTrack.tk_vtxInd()<<" track dz: "<<scoutingTrack.tk_dz()<<std::endl;
  //std::cout<<"track reference point: "<<scoutingTrack.tk_vx()<<" "<<scoutingTrack.tk_vy()<<" "<<scoutingTrack.tk_vz()<<std::endl;
  reco::TrackBase::TrackAlgorithm algo(reco::TrackBase::undefAlgorithm);  // undefined
  reco::TrackBase::TrackQuality quality(reco::TrackBase::confirmed);      // confirmed

  // the rests are default: t0 = 0, beta = 0, covt0t0 = -1, covbetabeta = -1

  reco::Track recoTrack(chi2, ndof, referencePoint, momentum, charge, cov, algo, quality);

  return recoTrack;
}

int HLTScoutingUnpackProducer::findCompatibleScoutingTrack(
    edm::Handle<Run3ScoutingParticleCollection> const& scoutingParticleCollection,
    Run3ScoutingTrack const& scoutingTrack) {
  int index = 0;

  auto is_close = [](float a, float b, float relative_tolerance) -> bool {
    return fabs(a - b) <= relative_tolerance * fmax(fabs(a), fabs(b));
  };

  for (size_t pf_index = 0; pf_index < scoutingParticleCollection->size(); pf_index++) {
    auto& scoutingPFCandidate = scoutingParticleCollection->at(pf_index);
    // retrieve track parameters from PFCandidate
    float pt = scoutingPFCandidate.trk_pt();
    float eta = scoutingPFCandidate.trk_eta();
    float phi = scoutingPFCandidate.trk_phi();
    if (scoutingPFCandidate.relative_trk_vars()) {
      pt += scoutingPFCandidate.pt();
      eta += scoutingPFCandidate.eta();
      phi += scoutingPFCandidate.phi();
    }
    float normchi2 = scoutingPFCandidate.normchi2();

    if (is_close(normchi2, scoutingTrack.tk_chi2() / scoutingTrack.tk_ndof(), 0.001) &&
        is_close(pt, scoutingTrack.tk_pt(), 0.000001) && is_close(eta, scoutingTrack.tk_eta(), 0.000001) &&
        deltaPhi(phi, scoutingTrack.tk_phi()) <= 0.00001) {
      //printf("matched track!\n");
      return index;
    }
    index++;
  }
  //printf("did not match track!\n");
  return -1;
}

// example from https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/src/PackedCandidate.cc#L219
void HLTScoutingUnpackProducer::buildHitPattern(Run3ScoutingParticle const& scoutingPFCandidate,
                                                Run3ScoutingTrack const& scoutingTrack,
                                                reco::Track& recoTrack) {
  // retrieve information from scoutingPFCandidate
  auto lost_inner_hits =
      static_cast<pat::PackedCandidate::LostInnerHits>(static_cast<int8_t>(scoutingPFCandidate.lostInnerHits()));

  // retrieve information from scoutingTrack
  int number_of_pixel_hits = scoutingTrack.tk_nValidPixelHits();
  int number_of_strip_hits = scoutingTrack.tk_nValidStripHits();
  int number_of_tracker_layers_with_measurement =
      scoutingTrack.tk_nTrackerLayersWithMeasurement();  // pixel + strip number of layers?
  int number_of_hits = number_of_pixel_hits + number_of_strip_hits;

  // pixel has 7 layers: 4 in BPIX + 3 FPIX
  int number_of_pxb_layers = 4;
  int number_of_pxf_layers = 3;
  int number_of_pixel_layers = number_of_pxb_layers + number_of_pxf_layers;
  // we assume the number of pixel layers with measurement is min(number_of_pixel_layers, number_of_layers_with_measurement)
  // in fact, even if we have number of layers with measurement more than number of pixel layers, we can have pixel layers with measurement less than number of pixel layers
  // but we don't have access to this information, so we will make this assumption here
  int number_of_pixel_layers_with_measurement = number_of_tracker_layers_with_measurement > number_of_pixel_layers
                                                    ? number_of_pixel_layers
                                                    : number_of_tracker_layers_with_measurement;
  // number of pixel layers with measurement cannot exceed number of pixel hits
  number_of_pixel_layers_with_measurement = number_of_pixel_layers_with_measurement > number_of_pixel_hits
                                                ? number_of_pixel_hits
                                                : number_of_pixel_layers_with_measurement;
  if ((number_of_pixel_layers_with_measurement == number_of_tracker_layers_with_measurement) &&
      (number_of_strip_hits > 0)) {
    number_of_pixel_layers_with_measurement -= 1;
  }

  // strip has 22 layers: 4 TIB + 3 TID + 6 TOB + 9 TEC
  int number_of_tib_layers = 4;
  int number_of_tid_layers = 3;
  int number_of_tob_layers = 6;
  int number_of_tec_layers = 9;
  int number_of_strip_layers =
      number_of_tib_layers + number_of_tob_layers + number_of_tid_layers + number_of_tec_layers;

  // assume the remaining layers with measuremnet are in strip
  int number_of_strip_layers_with_measurement =
      number_of_tracker_layers_with_measurement - number_of_pixel_layers_with_measurement;
  // number of strip layers with measurement cannot exceed number of strip layers
  // in fact, this it exceeds, we should throw error
  number_of_strip_layers_with_measurement = number_of_strip_layers_with_measurement > number_of_strip_layers
                                                ? number_of_strip_layers
                                                : number_of_strip_layers_with_measurement;
  // number of strip layers with measurement cannot exceed number of strip hits
  number_of_strip_layers_with_measurement = number_of_strip_layers_with_measurement > number_of_strip_hits
                                                ? number_of_strip_hits
                                                : number_of_strip_layers_with_measurement;

  // now, proceed similar to PackedCandidate
  uint16_t first_hit = 0;  // assume 0 since ScoutingTrack does not save HitPattern

  int i = 0;
  if (number_of_pixel_layers_with_measurement > 0) {
    if (first_hit == 0) {
      if (lost_inner_hits == pat::PackedCandidate::validHitInFirstPixelBarrelLayer) {
        recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, 1, 0, TrackingRecHit::valid);
        i = 1;
      }
    } else {
      recoTrack.appendHitPattern(first_hit, TrackingRecHit::valid);
      if (reco::HitPattern::pixelHitFilter(first_hit)) {
        i = 1;
      }
    }
  }
  // add hits to match the number of layers and validHitInFirstPixelBarrelLayer
  if (lost_inner_hits == pat::PackedCandidate::validHitInFirstPixelBarrelLayer) {
    for (; i < number_of_pixel_layers_with_measurement; i++) {
      if (i <= 3) {
        recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, i + 1, 0, TrackingRecHit::valid);
      } else {
        recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelEndcap, i - 3, 0, TrackingRecHit::valid);
      }
    }
  } else {
    int i_offset = 0;
    if (first_hit != 0 && reco::HitPattern::pixelHitFilter(first_hit)) {
      i_offset = reco::HitPattern::getLayer(first_hit);
      if (reco::HitPattern::getSubStructure(first_hit) == PixelSubdetector::PixelEndcap) {
        i_offset += 3;
      }
    } else {
      i_offset = 1;
    }
    for (; i < number_of_pixel_layers_with_measurement; i++) {
      if (i + i_offset <= 2) {
        recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, i + i_offset + 1, 0, TrackingRecHit::valid);
      } else {
        recoTrack.appendTrackerHitPattern(
            PixelSubdetector::PixelEndcap, i + i_offset + 1 - 3, 0, TrackingRecHit::valid);
      }
    }
  }

  // add extra hits (overlaps, etc), all on the first layer with a hit
  // to avoid increasing the layer count
  for (; i < number_of_pixel_hits; i++) {
    if (first_hit != 0 && reco::HitPattern::pixelHitFilter(first_hit)) {
      recoTrack.appendTrackerHitPattern(reco::HitPattern::getSubStructure(first_hit),
                                        reco::HitPattern::getLayer(first_hit),
                                        0,
                                        TrackingRecHit::valid);
    } else {
      recoTrack.appendTrackerHitPattern(
          PixelSubdetector::PixelBarrel,
          lost_inner_hits == pat::PackedCandidate::validHitInFirstPixelBarrelLayer ? 1 : 2,
          0,
          TrackingRecHit::valid);
    }
  }

  // now start adding strip layers, putting one hit on each layer
  // we don't know what the layers where, so we just start with 4 TIB, then 3 TID, then 6 TOB, and then 9 TEC
  // note the comment in PackedCandidate is incorrect
  if (first_hit != 0 && reco::HitPattern::stripHitFilter(first_hit)) {
    i += 1;
  }
  int sl_offset = 0;
  if (first_hit != 0 && reco::HitPattern::stripHitFilter(first_hit)) {
    sl_offset = reco::HitPattern::stripHitFilter(first_hit) - 1;
    if (reco::HitPattern::getSubStructure(first_hit) == StripSubdetector::TID)
      sl_offset += 4;
    if (reco::HitPattern::getSubStructure(first_hit) == StripSubdetector::TOB)
      sl_offset += 7;
    if (reco::HitPattern::getSubStructure(first_hit) == StripSubdetector::TEC)
      sl_offset += 13;
  }

  for (int sl = sl_offset; sl < number_of_strip_layers_with_measurement + sl_offset; sl++, i++) {
    if (sl < 4) {
      recoTrack.appendTrackerHitPattern(StripSubdetector::TIB, sl + 1, 1, TrackingRecHit::valid);
    } else if (sl < 4 + 3) {
      recoTrack.appendTrackerHitPattern(StripSubdetector::TID, (sl - 4) + 1, 1, TrackingRecHit::valid);
    } else if (sl < 7 + 6) {
      recoTrack.appendTrackerHitPattern(StripSubdetector::TOB, (sl - 7) + 1, 1, TrackingRecHit::valid);
    } else if (sl < 13 + 9) {
      recoTrack.appendTrackerHitPattern(StripSubdetector::TEC, (sl - 13) + 1, 1, TrackingRecHit::valid);
    } else {
      break;
    }
  }

  for (; i < number_of_hits; i++) {
    if (reco::HitPattern::stripHitFilter(first_hit)) {
      recoTrack.appendTrackerHitPattern(reco::HitPattern::getSubStructure(first_hit),
                                        reco::HitPattern::getLayer(first_hit),
                                        1,
                                        TrackingRecHit::valid);
    } else {
      recoTrack.appendTrackerHitPattern(StripSubdetector::TIB, 1, 1, TrackingRecHit::valid);
    }
  }

  switch (lost_inner_hits) {
    case pat::PackedCandidate::validHitInFirstPixelBarrelLayer:
      break;
    case pat::PackedCandidate::noLostInnerHits:
      break;
    case pat::PackedCandidate::oneLostInnerHit:
      recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, 1, 0, TrackingRecHit::missing_inner);
      break;
    case pat::PackedCandidate::moreLostInnerHits:
      recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, 1, 0, TrackingRecHit::missing_inner);
      recoTrack.appendTrackerHitPattern(PixelSubdetector::PixelBarrel, 2, 0, TrackingRecHit::missing_inner);
      break;
  };
  //if((number_of_pixel_hits!=recoTrack.hitPattern().numberOfValidPixelHits()) || (number_of_strip_hits!=recoTrack.hitPattern().numberOfValidStripHits()) || (number_of_tracker_layers_with_measurement!=recoTrack.hitPattern().trackerLayersWithMeasurement())){
  //std::cout << "pixel hits expected: " << number_of_pixel_hits << " get: " << recoTrack.hitPattern().numberOfValidPixelHits() << " strip hits expected: " << number_of_strip_hits << " get: " << recoTrack.hitPattern().numberOfValidStripHits();
  //std::cout << " layers with measurement expect: " << number_of_tracker_layers_with_measurement << " get: " << recoTrack.hitPattern().trackerLayersWithMeasurement() << std::endl;
  //std::cout<<"missing inner hits: "<<(int) lost_inner_hits<<" number of pixel layers with measurement: "<<number_of_pixel_layers_with_measurement<<" number of strip layers with measurement: "<<number_of_strip_layers_with_measurement<<std::endl;
  //}
}

template <typename RecoObjectType, typename ScoutingObjectType>
void HLTScoutingUnpackProducer::putWithRef(
    edm::Event& iEvent,
    std::string const& name,
    std::unique_ptr<std::vector<RecoObjectType>>& recoObject_collection_ptr,
    std::unique_ptr<RefCollection<ScoutingObjectType>>& scoutingObjectRef_collection_ptr) {
  auto recoObject_collection_handle = iEvent.put(std::move(recoObject_collection_ptr), name);

  std::unique_ptr<RefMap<ScoutingObjectType>> refmap_to_scouting(new RefMap<ScoutingObjectType>());
  typename RefMap<ScoutingObjectType>::Filler filler_refmap_to_scouting(*refmap_to_scouting);
  filler_refmap_to_scouting.insert(
      recoObject_collection_handle, scoutingObjectRef_collection_ptr->begin(), scoutingObjectRef_collection_ptr->end());
  filler_refmap_to_scouting.fill();
  iEvent.put(std::move(refmap_to_scouting), name + REF_TO_SCOUTING_LABEL_SUFFIX_);
}

void HLTScoutingUnpackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("scoutingTrack", edm::InputTag("hltScoutingTrack"));
  desc.add<edm::InputTag>("scoutingPrimaryVertex", edm::InputTag("hltScoutingPrimaryVertex"));
  desc.add<edm::InputTag>("scoutingParticle", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("pfCand", edm::InputTag(""));
  desc.add<edm::InputTag>("lostTrack", edm::InputTag(""));
  desc.add<bool>("isScouting", true);

  desc.add<bool>("producePFCHSCandidate", false);
  //desc.add<bool>("producePFSKCandidate", false);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HLTScoutingUnpackProducer);

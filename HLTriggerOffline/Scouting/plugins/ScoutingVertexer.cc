// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingVertexer
// Class:      ScoutingVertexer
//
/**\class ScoutingVertexer ScoutingVertexer.cc Run3ScoutingAnalysisTools/ScoutingVertexer/plugins/ScoutingVertexer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bruno Lopes
//         Created:  Wed, 04 Sep 2024 12:37:22 GMT
//
//

// system include files
#include <memory>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/ESGetToken.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//Scouting data formats
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

//Vertex tools
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TH1.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace edm;

//
// class declaration
//

class ScoutingVertexer : public edm::stream::EDProducer<> {
public:
  ~ScoutingVertexer() override;
  explicit ScoutingVertexer(edm::ParameterSet const& params);

private:
  typedef std::set<reco::TrackRef> track_set;
  typedef std::vector<reco::TrackRef> track_vec;

  void beginStream(edm::StreamID) override;
  void endStream() override;
  void produce(edm::Event&, const edm::EventSetup&) override;

  TH1D* h_dxyErr_weighted_sum_barrel;
  TH1D* h_dszErr_weighted_sum_barrel;
  TH1D* h_dszdxyCov_weighted_sum_barrel;
  TH1D* h_dxyErr_weighted_sq_sum_barrel;
  TH1D* h_dszErr_weighted_sq_sum_barrel;
  TH1D* h_dszdxyCov_weighted_sq_sum_barrel;
  TH1D* h_weight_sum_barrel;
  TH1D* h_weight_sq_sum_barrel;

  TH1D* h_dxyErr_weighted_sum_disk;
  TH1D* h_dszErr_weighted_sum_disk;
  TH1D* h_dszdxyCov_weighted_sum_disk;
  TH1D* h_dxyErr_weighted_sq_sum_disk;
  TH1D* h_dszErr_weighted_sq_sum_disk;
  TH1D* h_dszdxyCov_weighted_sq_sum_disk;
  TH1D* h_weight_sum_disk;
  TH1D* h_weight_sq_sum_disk;

  const double pt_min_cut;
  const double dxySig_min_cut;
  const double dxySig_max_cut;
  const int npixelHits_min_cut;
  const int nstripHits_min_cut;
  const int ntrackerLayers_min_cut;
  const int n_tracks_per_seed_vertex;
  const double max_seed_vertex_chi2;
  const bool use_2d_vertex_dist;
  const bool use_2d_track_dist;
  const bool remove_one_track_at_a_time;
  const double merge_shared_dist;
  const double merge_shared_sig;
  const double max_track_vertex_dist;
  const double max_track_vertex_sig;
  const double min_track_vertex_sig_to_remove;
  const bool resolve_split_vertices_loose;
  const bool resolve_split_vertices_tight;
  const double merge_anyway_sig;
  const double merge_anyway_dist;
  const double max_nm1_refit_dist3;
  const double max_nm1_refit_distz;
  const int max_nm1_refit_count;
  const bool investigate_merged_vertices;
  const bool verbose;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_token;
  const edm::EDGetTokenT<std::vector<reco::Track>> seed_tracks_token_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> token_builder;

  edm::EDPutTokenT<reco::VertexCollection> putToken_;
  edm::EDPutTokenT<edm::ValueMap<std::pair<float, float>>> vertexShiftZToken_;
  edm::EDPutTokenT<edm::ValueMap<std::pair<float, float>>> vertexShift3DToken_;

  // ----------member data ---------------------------

  VertexDistanceXY vertex_dist_2d;
  VertexDistance3D vertex_dist_3d;

  bool is_track_subset(const track_set& a, const track_set& b) const {
    bool is_subset = true;
    const track_set& smaller = a.size() <= b.size() ? a : b;
    const track_set& bigger = a.size() <= b.size() ? b : a;

    for (auto t : smaller)
      if (bigger.count(t) < 1) {
        is_subset = false;
        break;
      }

    return is_subset;
  }

  Measurement1D vertex_dist(const reco::Vertex& v0, const reco::Vertex& v1) {
    if (use_2d_vertex_dist)
      return vertex_dist_2d.distance(v0, v1);
    else
      return vertex_dist_3d.distance(v0, v1);
  }

  track_set vertex_track_set(const reco::Vertex& v, const double min_weight = 0.5) const {
    track_set result;

    for (auto it = v.tracks_begin(), ite = v.tracks_end(); it != ite; ++it) {
      const double w = v.trackWeight(*it);
      const bool use = w >= min_weight;
      assert(use);
      if (use)
        result.insert(it->castTo<reco::TrackRef>());
    }

    return result;
  }

  std::pair<bool, Measurement1D> track_dist(const reco::TransientTrack& t, const reco::Vertex& v) const {
    if (use_2d_track_dist)
      return IPTools::absoluteTransverseImpactParameter(t, v);
    else
      return IPTools::absoluteImpactParameter3D(t, v);
  }

  track_vec vertex_track_vec(const reco::Vertex& v, const double min_weight = 0.5) const {
    track_set s = vertex_track_set(v, min_weight);
    return track_vec(s.begin(), s.end());
  }

  template <typename T>
  void print_track_set(const T& ts) const {
    for (auto r : ts)
      printf(" %u", r.key());
  }

  template <typename T>
  void print_track_set(const T& ts, const reco::Vertex& v) const {
    for (auto r : ts)
      printf(" %u%s", r.key(), (v.trackWeight(r) < 0.5 ? "!" : ""));
  }

  void print_track_set(const reco::Vertex& v) const {
    for (auto r = v.tracks_begin(), re = v.tracks_end(); r != re; ++r)
      printf(" %lu%s", r->key(), (v.trackWeight(*r) < 0.5 ? "!" : ""));
  }
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

ScoutingVertexer::ScoutingVertexer(edm::ParameterSet const& params)
    : pt_min_cut(params.getParameter<double>("pt_min_cut")),
      dxySig_min_cut(params.getParameter<double>("dxySig_min_cut")),
      dxySig_max_cut(params.getParameter<double>("dxySig_max_cut")),
      npixelHits_min_cut(params.getParameter<int>("npixelHits_min_cut")),
      nstripHits_min_cut(params.getParameter<int>("nstripHits_min_cut")),
      ntrackerLayers_min_cut(params.getParameter<int>("ntrackerLayers_min_cut")),
      n_tracks_per_seed_vertex(params.getParameter<int>("n_tracks_per_seed_vertex")),
      max_seed_vertex_chi2(params.getParameter<double>("max_seed_vertex_chi2")),
      use_2d_vertex_dist(params.getParameter<bool>("use_2d_vertex_dist")),
      use_2d_track_dist(params.getParameter<bool>("use_2d_track_dist")),
      remove_one_track_at_a_time(params.getParameter<bool>("remove_one_track_at_a_time")),
      merge_shared_dist(params.getParameter<double>("merge_shared_dist")),
      merge_shared_sig(params.getParameter<double>("merge_shared_sig")),
      max_track_vertex_dist(params.getParameter<double>("max_track_vertex_dist")),
      max_track_vertex_sig(params.getParameter<double>("max_track_vertex_sig")),
      min_track_vertex_sig_to_remove(params.getParameter<double>("min_track_vertex_sig_to_remove")),
      resolve_split_vertices_loose(params.getParameter<bool>("resolve_split_vertices_loose")),
      resolve_split_vertices_tight(params.getParameter<bool>("resolve_split_vertices_tight")),
      merge_anyway_sig(params.getParameter<double>("merge_anyway_sig")),
      merge_anyway_dist(params.getParameter<double>("merge_anyway_dist")),
      max_nm1_refit_dist3(params.getParameter<double>("max_nm1_refit_dist3")),
      max_nm1_refit_distz(params.getParameter<double>("max_nm1_refit_distz")),
      max_nm1_refit_count(params.getParameter<int>("max_nm1_refit_count")),
      investigate_merged_vertices(params.getParameter<bool>("investigate_merged_vertices")),
      verbose(params.getParameter<bool>("verbose")),

      beamspot_token(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamspot_src"))),
      seed_tracks_token_(consumes(params.getParameter<edm::InputTag>("seed_tracks_src"))),
      token_builder(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),

      putToken_{produces()} {
  vertexShiftZToken_ = produces<edm::ValueMap<std::pair<float, float>>>("vtxZShift");
  vertexShift3DToken_ = produces<edm::ValueMap<std::pair<float, float>>>("vtx3DShift");
}

ScoutingVertexer::~ScoutingVertexer() {}

//
// member functions
//

KalmanVertexFitter kv_reco;
std::vector<TransientVertex> kv_reco_dropin(std::vector<reco::TransientTrack>& ttks) {
  if (ttks.size() < 2)
    return std::vector<TransientVertex>();
  std::vector<TransientVertex> v(1, kv_reco.vertex(ttks));
  if (v[0].normalisedChiSquared() > 5)
    return std::vector<TransientVertex>();
  return v;
}

// ------------ method called to produce the data  ------------
void ScoutingVertexer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //////////////////////////////////////////////////////////////////////
  // DataFormats setup and track preselection
  //////////////////////////////////////////////////////////////////////

  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(beamspot_token, beamspot);
  const double bsx = beamspot->position().x();
  const double bsy = beamspot->position().y();
  const double bsz = beamspot->position().z();
  const reco::Vertex fake_bs_vtx(beamspot->position(), beamspot->covariance3D());

  //Get the Transient Track Builder
  auto const& tt_builder = iSetup.getData(token_builder);

  //Get the reco tracks from the events
  edm::Handle<std::vector<reco::Track>> seed_track_handle;
  iEvent.getByToken(seed_tracks_token_, seed_track_handle);

  // Build the references to the tracks
  std::vector<reco::TrackRef> seed_track_refs;
  std::map<reco::TrackRef, size_t> seed_track_index_map;

  for (size_t i_tk = 0; i_tk < seed_track_handle->size(); i_tk++) {
    const edm::Ref<reco::TrackCollection> tk_ref(seed_track_handle, i_tk);
    reco::TransientTrack ttk = tt_builder.build(tk_ref);
    std::pair<bool, Measurement1D> ttk_dist = IPTools::absoluteTransverseImpactParameter(ttk, fake_bs_vtx);
    //std::pair<bool, Measurement1D> ttk_dist = track_dist(ttk, fake_bs_vtx);
    float IP_sig = ttk_dist.second.significance();
    if ((tk_ref->pt() > pt_min_cut) && (tk_ref->hitPattern().numberOfValidPixelHits() > npixelHits_min_cut) &&
        (tk_ref->hitPattern().numberOfValidStripHits() > nstripHits_min_cut) &&
        (tk_ref->hitPattern().trackerLayersWithMeasurement() > ntrackerLayers_min_cut) && (fabs(tk_ref->eta()) < 2.4)) {
      //if ((tk_ref->pt()>0.9) && (fabs(tk_ref->eta())<2.4)){
      if (fabs(tk_ref->eta()) < 1.5) {
        h_dxyErr_weighted_sum_barrel->Fill(tk_ref->pt(), tk_ref->dxyError());
        h_dszErr_weighted_sum_barrel->Fill(tk_ref->pt(), tk_ref->dszError());
        h_dszdxyCov_weighted_sum_barrel->Fill(tk_ref->pt(), fabs(tk_ref->covariance(3, 4)));
        h_dxyErr_weighted_sq_sum_barrel->Fill(tk_ref->pt(), pow(tk_ref->dxyError(), 2));
        h_dszErr_weighted_sq_sum_barrel->Fill(tk_ref->pt(), pow(tk_ref->dszError(), 2));
        h_dszdxyCov_weighted_sq_sum_barrel->Fill(tk_ref->pt(), pow(tk_ref->covariance(3, 4), 2));
        h_weight_sum_barrel->Fill(tk_ref->pt());
        h_weight_sq_sum_barrel->Fill(tk_ref->pt());
      } else {
        h_dxyErr_weighted_sum_disk->Fill(tk_ref->pt(), tk_ref->dxyError());
        h_dszErr_weighted_sum_disk->Fill(tk_ref->pt(), tk_ref->dszError());
        h_dszdxyCov_weighted_sum_disk->Fill(tk_ref->pt(), fabs(tk_ref->covariance(3, 4)));
        h_dxyErr_weighted_sq_sum_disk->Fill(tk_ref->pt(), pow(tk_ref->dxyError(), 2));
        h_dszErr_weighted_sq_sum_disk->Fill(tk_ref->pt(), pow(tk_ref->dszError(), 2));
        h_dszdxyCov_weighted_sq_sum_disk->Fill(tk_ref->pt(), pow(tk_ref->covariance(3, 4), 2));
        h_weight_sum_disk->Fill(tk_ref->pt());
        h_weight_sq_sum_disk->Fill(tk_ref->pt());
      }

      if ((((dxySig_max_cut > 0) && (IP_sig < dxySig_max_cut)) || (dxySig_max_cut <= 0)) && (IP_sig > dxySig_min_cut)) {
        //if(IP_sig > 2){
        seed_track_refs.push_back(tk_ref);
        seed_track_index_map[tk_ref] = i_tk;
      }
    }
    //if ((IP_sig > 4) && (tk_ref->pt()>0.9)) seed_track_refs.push_back(tk_ref);
    if (verbose)
      printf("Build track references. IP_sig = %f\n", IP_sig);
  }

  //Build transient tracks from reco tracks
  std::vector<reco::TransientTrack> seed_tracks;

  std::map<reco::TrackRef, size_t> seed_track_ref_map;
  for (const reco::TrackRef& tk : seed_track_refs) {
    seed_tracks.push_back(tt_builder.build(tk));
    seed_track_ref_map[tk] = seed_tracks.size() - 1;
  }

  std::vector<std::pair<float, float>> shiftZVec(seed_track_handle->size(), std::make_pair(-999.0, -999.0));
  std::vector<std::pair<float, float>> shift3DVec(seed_track_handle->size(), std::make_pair(-999.0, -999.0));

  //////////////////////////////////////////////////////////////////////
  // Form seed vertices from all pairs of tracks whose vertex fit
  // passes cuts.
  //////////////////////////////////////////////////////////////////////

  const size_t ntk = seed_tracks.size();

  std::unique_ptr<reco::VertexCollection> vertices(new reco::VertexCollection);
  auto outMapZ = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler fillerZ(*outMapZ);

  auto outMap3D = std::make_unique<edm::ValueMap<std::pair<float, float>>>();
  edm::ValueMap<std::pair<float, float>>::Filler filler3D(*outMap3D);

  if (ntk == 0) {
    iEvent.emplace(putToken_, std::move(*vertices));
    fillerZ.insert(seed_track_handle, shiftZVec.begin(), shiftZVec.end());
    fillerZ.fill();
    iEvent.put(std::move(outMapZ), "vtxZShift");
    filler3D.insert(seed_track_handle, shift3DVec.begin(), shift3DVec.end());
    filler3D.fill();
    iEvent.put(std::move(outMap3D), "vtx3DShift");
    return;
  }

  std::vector<size_t> itks(n_tracks_per_seed_vertex, 0);

  auto try_seed_vertex = [&]() {
    std::vector<reco::TransientTrack> ttks(n_tracks_per_seed_vertex);
    for (int i = 0; i < n_tracks_per_seed_vertex; ++i)
      ttks[i] = seed_tracks[itks[i]];

    TransientVertex seed_vertex = kv_reco.vertex(ttks);
    if (seed_vertex.isValid() && seed_vertex.normalisedChiSquared() < max_seed_vertex_chi2) {
      vertices->push_back(reco::Vertex(seed_vertex));

      if (verbose) {
        const reco::Vertex& v = vertices->back();
        const double vchi2 = v.normalizedChi2();
        const double vndof = v.ndof();
        const double vx = v.position().x() - bsx;
        const double vy = v.position().y() - bsy;
        const double vz = v.position().z() - bsz;
        const double phi = atan2(vy, vx);
        const double rho = sqrt(vx * vx + vy * vy);
        const double r = sqrt(vx * vx + vy * vy + vz * vz);

        printf("from tracks");
        for (auto itk : itks)
          printf(" %lu", itk);
        printf(
            ": vertex #%3lu: chi2/dof: %7.3f dof: %7.3f pos: <%7.3f, %7.3f, %7.3f>  rho: %7.3f  phi: %7.3f  r: %7.3f\n",
            vertices->size() - 1,
            vchi2,
            vndof,
            vx,
            vy,
            vz,
            rho,
            phi,
            r);
      }
    }
  };

  // ha
  for (size_t itk = 0; itk < ntk; ++itk) {
    itks[0] = itk;
    for (size_t jtk = itk + 1; jtk < ntk; ++jtk) {
      itks[1] = jtk;
      if (n_tracks_per_seed_vertex == 2) {
        try_seed_vertex();
        continue;
      }
      for (size_t ktk = jtk + 1; ktk < ntk; ++ktk) {
        itks[2] = ktk;
        if (n_tracks_per_seed_vertex == 3) {
          try_seed_vertex();
          continue;
        }
        for (size_t ltk = ktk + 1; ltk < ntk; ++ltk) {
          itks[3] = ltk;
          if (n_tracks_per_seed_vertex == 4) {
            try_seed_vertex();
            continue;
          }
          for (size_t mtk = ltk + 1; mtk < ntk; ++mtk) {
            itks[4] = mtk;
            try_seed_vertex();
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Take care of track sharing. If a track is in two vertices, and
  // the vertices are "close", refit the tracks from the two together
  // as one vertex. If the vertices are not close, keep the track in
  // the vertex to which it is "closer".
  //////////////////////////////////////////////////////////////////////

  //printf("entering the track sharing part\n");

  track_set discarded_tracks;
  int n_resets = 0;
  int n_onetracks = 0;
  std::vector<reco::Vertex>::iterator v[2];

  size_t ivtx[2];

  for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
    track_set tracks[2];
    ivtx[0] = v[0] - vertices->begin();
    tracks[0] = vertex_track_set(*v[0]);

    if (tracks[0].size() < 2) {
      if (verbose)
        printf("track-sharing: vertex-0 #%lu is down to one track, junking it\n", ivtx[0]);
      v[0] = vertices->erase(v[0]) - 1;
      ++n_onetracks;
      continue;
    }

    bool duplicate = false;
    bool merge = false;
    bool refit = false;
    track_set tracks_to_remove_in_refit[2];

    for (v[1] = v[0] + 1; v[1] != vertices->end(); ++v[1]) {
      ivtx[1] = v[1] - vertices->begin();
      tracks[1] = vertex_track_set(*v[1]);

      if (tracks[1].size() < 2) {
        if (verbose)
          printf("track-sharing: vertex-1 #%lu is down to one track, junking it\n", ivtx[1]);
        v[1] = vertices->erase(v[1]) - 1;
        ++n_onetracks;
        continue;
      }

      if (verbose) {
        printf("track-sharing: # vertices = %lu. considering vertices #%lu (chi2/dof %.3f, track set",
               vertices->size(),
               ivtx[0],
               v[0]->chi2() / v[0]->ndof());
        print_track_set(tracks[0], *v[0]);
        printf(") and #%lu (chi2/dof %.3f, track set", ivtx[1], v[1]->chi2() / v[1]->ndof());
        print_track_set(tracks[1], *v[1]);
        printf("):\n");
      }

      if (is_track_subset(tracks[0], tracks[1])) {
        if (verbose)
          printf("   subset/duplicate vertices %lu and %lu, erasing second and starting over\n", ivtx[0], ivtx[1]);
        duplicate = true;
        break;
      }

      std::vector<reco::TrackRef> shared_tracks;
      for (auto tk : tracks[0])
        if (tracks[1].count(tk) > 0)
          shared_tracks.push_back(tk);

      if (verbose) {
        if (shared_tracks.size()) {
          printf("   shared tracks are: ");
          print_track_set(shared_tracks);
          printf("\n");
        } else
          printf("   no shared tracks\n");
      }

      if (shared_tracks.size() > 0) {
        Measurement1D v_dist = vertex_dist(*v[0], *v[1]);

        if (verbose)
          printf(
              "   vertex dist (2d? %i) %7.3f  sig %7.3f\n", use_2d_vertex_dist, v_dist.value(), v_dist.significance());

        if (v_dist.value() < merge_shared_dist || v_dist.significance() < merge_shared_sig) {
          if (verbose)
            printf("          dist < %7.3f || sig < %7.3f, will try using merge result first before arbitration\n",
                   merge_shared_dist,
                   merge_shared_sig);
          merge = true;
        } else
          refit = true;

        for (auto tk : shared_tracks) {
          const reco::TransientTrack& ttk = seed_tracks[seed_track_ref_map[tk]];
          std::pair<bool, Measurement1D> t_dist_0 = track_dist(ttk, *v[0]);
          std::pair<bool, Measurement1D> t_dist_1 = track_dist(ttk, *v[1]);

          t_dist_0.first = t_dist_0.first && (t_dist_0.second.value() < max_track_vertex_dist ||
                                              t_dist_0.second.significance() < max_track_vertex_sig);
          t_dist_1.first = t_dist_1.first && (t_dist_1.second.value() < max_track_vertex_dist ||
                                              t_dist_1.second.significance() < max_track_vertex_sig);
          bool remove_from_0 = !t_dist_0.first;
          bool remove_from_1 = !t_dist_1.first;
          if (t_dist_0.second.significance() < min_track_vertex_sig_to_remove &&
              t_dist_1.second.significance() < min_track_vertex_sig_to_remove) {
            if (tracks[0].size() > tracks[1].size())
              remove_from_1 = true;
            else
              remove_from_0 = true;
          } else if (t_dist_0.second.significance() < t_dist_1.second.significance())
            remove_from_1 = true;
          else
            remove_from_0 = true;

          if (remove_from_0)
            tracks_to_remove_in_refit[0].insert(tk);
          if (remove_from_1)
            tracks_to_remove_in_refit[1].insert(tk);

          if (remove_one_track_at_a_time)
            break;
        }

        break;
      }
    }

    if (duplicate) {
      vertices->erase(v[1]);
    }

    else if (merge) {
      track_set tracks_to_fit;
      for (int i = 0; i < 2; ++i)
        for (auto tk : tracks[i])
          tracks_to_fit.insert(tk);

      std::vector<reco::TransientTrack> ttks;
      for (auto tk : tracks_to_fit)
        ttks.push_back(seed_tracks[seed_track_ref_map[tk]]);

      reco::VertexCollection new_vertices;

      for (const TransientVertex& tv : kv_reco_dropin(ttks)) {
        new_vertices.push_back(reco::Vertex(tv));
      }
      // If we got two new vertices, maybe it took A B and A C D and made a better one from B C D, and left a broken one A B! C! D!.
      // If we get one that is truly the merger of the track lists, great. If it is just something like A B , A C . A B C!, or we get nothing, then default to arbitration.
      if (new_vertices.size() > 1) {
        assert(new_vertices.size() == 2);
        *v[1] = reco::Vertex(new_vertices[1]);
        *v[0] = reco::Vertex(new_vertices[0]);

      } else if (new_vertices.size() == 1 && vertex_track_set(new_vertices[0], 0) == tracks_to_fit) {
        vertices->erase(v[1]);

        *v[0] = reco::Vertex(
            new_vertices[0]);  // ok to use v[0] after the erase(v[1]) because v[0] is by construction before v[1]

      } else
        refit = true;
    }
    if (refit) {
      bool erase[2] = {false};
      reco::Vertex vsave[2] = {*v[0], *v[1]};

      for (int i = 0; i < 2; ++i) {
        if (tracks_to_remove_in_refit[i].empty())
          continue;

        std::vector<reco::TransientTrack> ttks;
        for (auto tk : tracks[i])
          if (tracks_to_remove_in_refit[i].count(tk) == 0)
            ttks.push_back(seed_tracks[seed_track_ref_map[tk]]);

        reco::VertexCollection new_vertices;
        for (const TransientVertex& tv : kv_reco_dropin(ttks)) {
          new_vertices.push_back(reco::Vertex(tv));
        }
        if (new_vertices.size() == 1)
          *v[i] = new_vertices[0];
        else
          erase[i] = true;
      }

      if (erase[1])
        vertices->erase(v[1]);
      if (erase[0])
        vertices->erase(v[0]);
    }

    // If we changed the vertices at all, start loop over completely.
    if (duplicate || merge || refit) {
      //printf("duplicate = %d, merge = %d, refit = %d\n", duplicate, merge, refit);

      v[0] = vertices->begin() - 1;  // -1 because about to ++sv
      ++n_resets;

      //if (n_resets == 3000)
      //  throw "I'm dumb";
    }
  }

  //checkpoint

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Merge vertices that are still "close" in 2D, aka "loose" merging (typically off by default)
  //////////////////////////////////////////////////////////////////////////////////////////////

  if (resolve_split_vertices_loose) {
    if (merge_anyway_sig > 0 || merge_anyway_dist > 0) {
      //double v0x;
      //double v0y;
      //double phi0;

      for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
        //ivtx[0] = v[0] - vertices->begin();

        //double v1x;
        //double v1y;
        //double phi1;

        for (v[1] = v[0] + 1; v[1] != vertices->end(); ++v[1]) {
          //ivtx[1] = v[1] - vertices->begin();

          Measurement1D v_dist = vertex_dist(*v[0], *v[1]);

          //v0x = v[0]->x() - bsx;
          //v0y = v[0]->y() - bsy;
          //phi0 = atan2(v0y, v0x);
          //v1x = v[1]->x() - bsx;
          //v1y = v[1]->y() - bsy;
          //phi1 = atan2(v1y, v1x);

          if (v_dist.value() < merge_anyway_dist || v_dist.significance() < merge_anyway_sig) {
            std::vector<reco::TransientTrack> ttks;

            for (int i = 0; i < 2; ++i) {
              for (auto tk : vertex_track_set(*v[i])) {
                ttks.push_back(tt_builder.build(tk));
              }
            }

            reco::VertexCollection merged_vertices;
            for (const TransientVertex& tv : kv_reco_dropin(ttks)) {
              merged_vertices.push_back(reco::Vertex(tv));

              for (auto it = merged_vertices[0].tracks_begin(), ite = merged_vertices[0].tracks_end(); it != ite;
                   ++it) {
                reco::TransientTrack seed_track;
                seed_track = tt_builder.build(*it.operator*());
                std::pair<bool, Measurement1D> tk_vtx_dist = track_dist(seed_track, merged_vertices[0]);
              }
            }

            if (merged_vertices.size() == 1) {
              //std::cout << "check no mem out of ranges (before) : " << v[1] - vertices->begin() << std::endl;
              *v[0] = merged_vertices[0];
              //std::cout << "check no mem out of ranges (after) : " << v[1] - vertices->begin() << std::endl;

              v[1] = vertices->erase(v[1]) - 1;
            }
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Drop tracks that "move" the vertex too much by refitting without each track.
  //////////////////////////////////////////////////////////////////////

  if (max_nm1_refit_dist3 > 0 || max_nm1_refit_distz > 0) {
    for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
      const track_vec tks = vertex_track_vec(*v[0]);
      const size_t ntks = tks.size();
      if (ntks < 3)
        continue;

      std::vector<reco::TransientTrack> ttks(ntks - 1);
      for (size_t i = 0; i < ntks; ++i) {
        for (size_t j = 0; j < ntks; ++j)
          if (j != i)
            ttks[j - (j >= i)] = tt_builder.build(tks[j]);
        reco::Vertex vnm1(TransientVertex(kv_reco.vertex(ttks)));
        const double dist3_2 = (vnm1.x() - v[0]->x()) * (vnm1.x() - v[0]->x()) +
                               (vnm1.y() - v[0]->y()) * (vnm1.y() - v[0]->y()) +
                               (vnm1.z() - v[0]->z()) * (vnm1.z() - v[0]->z());
        const double distz = sqrt((vnm1.z() - v[0]->z()) * (vnm1.z() - v[0]->z()));
        shiftZVec[seed_track_index_map[tks[i]]] =
            std::make_pair(distz, sqrt(fabs(vnm1.covariance(2, 2) - v[0]->covariance(2, 2))));
        AlgebraicVector3 vDiff;
        vDiff[0] = (vnm1.x() - v[0]->x()) / sqrt(dist3_2);
        vDiff[1] = (vnm1.y() - v[0]->y()) / sqrt(dist3_2);
        vDiff[2] = (vnm1.z() - v[0]->z()) / sqrt(dist3_2);
        //AlgebraicSymMatrix33 error = fabs(vnm1.covariance()-v[0]->covariance());
        AlgebraicSymMatrix33 error1 = vnm1.covariance();
        AlgebraicSymMatrix33 error2 = v[0]->covariance();
        //double err2 = ROOT::Math::Similarity(error, vDiff);
        double err1_2 = ROOT::Math::Similarity(error1, vDiff);
        double err2_2 = ROOT::Math::Similarity(error2, vDiff);
        double err3D = sqrt(fabs(err1_2 - err2_2));
        //shift3DVec[seed_track_index_map[tks[i]]] = std::make_pair(sqrt(dist3_2),sqrt(err2));
        shift3DVec[seed_track_index_map[tks[i]]] = std::make_pair(sqrt(dist3_2), err3D);

        /*
        if (vnm1.chi2() < 0 ||
            (max_nm1_refit_dist3 > 0 && err3D > max_nm1_refit_dist3) ||
            (max_nm1_refit_distz > 0 && distz > max_nm1_refit_distz)) {

          *v[0] = vnm1;
          ++refit_count[iv];
          --v[0], --iv;
          break;
	  }*/
      }  //loop over tracks

      std::vector<reco::TransientTrack> ttks_shift;
      for (size_t i = 0; i < ntks; ++i) {
        //double shift3DErr = shift3DVec[seed_track_index_map[tks[i]]].second;
        double shift3DValue = shift3DVec[seed_track_index_map[tks[i]]].first;
        if (max_nm1_refit_dist3 > 0 && shift3DValue < max_nm1_refit_dist3) {
          ttks_shift.push_back(tt_builder.build(tks[i]));
        }
      }
      if (ttks_shift.size() != ntks) {
        if (ttks_shift.size() < 2) {
          v[0] = vertices->erase(v[0]) - 1;
          continue;
        }
        reco::Vertex vnm1_shift(TransientVertex(kv_reco.vertex(ttks_shift)));
        *v[0] = vnm1_shift;
      }
    }
    //some vertices after dz refiting have normalized chi2 > 5

    for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
      if ((*v[0]).normalizedChi2() > 5) {
        v[0] = vertices->erase(v[0]) - 1;
        continue;
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Merge every pair of output vertices that satisfy the following criteria to resolve split-vertices:
  //   - >=2trk/vtx
  //   - dBV > 100 um
  //   - |dPhi(vtx0,vtx1)| < 0.5
  //   - svdist2d < 300 um
  // Note that the merged vertex must pass chi2/dof < 5
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  if (resolve_split_vertices_tight) {
    reco::VertexCollection potential_merged_vertices;

    for (v[0] = vertices->begin(); v[0] != vertices->end(); ++v[0]) {
      track_set tracks[2];
      tracks[0] = vertex_track_set(*v[0]);

      bool merge = false;
      for (v[1] = v[0] + 1; v[1] != vertices->end(); ++v[1]) {
        if (vertices->size() >= 2 && v[0]->nTracks() >= 2 && v[1]->nTracks() >= 2) {
          tracks[1] = vertex_track_set(*v[1]);

          Measurement1D v_dist = vertex_dist_2d.distance(*v[0], *v[1]);

          Measurement1D dBV0_Meas1D = vertex_dist_2d.distance(*v[0], fake_bs_vtx);
          double dBV0 = dBV0_Meas1D.value();

          Measurement1D dBV1_Meas1D = vertex_dist_2d.distance(*v[1], fake_bs_vtx);
          double dBV1 = dBV1_Meas1D.value();

          double v0x = v[0]->x() - bsx;
          double v0y = v[0]->y() - bsy;

          double phi0 = atan2(v0y, v0x);

          double v1x = v[1]->x() - bsx;
          double v1y = v[1]->y() - bsy;

          double phi1 = atan2(v1y, v1x);

          if (fabs(reco::deltaPhi(phi0, phi1)) < 0.5 && v_dist.value() < 0.0500 && dBV0 > 0.0100 && dBV1 > 0.0100) {
            track_set tracks_to_fit;
            for (int i = 0; i < 2; ++i)
              for (auto tk : tracks[i])
                tracks_to_fit.insert(tk);
            std::vector<reco::TransientTrack> ttks;
            for (auto tk : tracks_to_fit)
              ttks.push_back(seed_tracks[seed_track_ref_map[tk]]);

            if (investigate_merged_vertices) {
              std::vector<TransientVertex> tv(1, kv_reco.vertex(ttks));
              potential_merged_vertices.push_back(reco::Vertex(tv[0]));
              //std::cout << "ntrack in potental merged: " << potential_merged_vertices.back().nTracks() << std::endl;
            }

            reco::VertexCollection merged_vertices;
            for (const TransientVertex& tv : kv_reco_dropin(ttks)) {
              merged_vertices.push_back(reco::Vertex(tv));
            }

            if (merged_vertices.size() == 1 && vertex_track_set(merged_vertices[0], 0) == tracks_to_fit) {
              merge = true;

              v[1] = vertices->erase(v[1]) - 1;  // (1) erase and point the iterator at the previous entry
              *v[0] = reco::Vertex(
                  merged_vertices
                      [0]);  // (2) updated v[0] (ok to use v[0] after the erase(v[1]) because v[0] is by construction before v[1])
            }
          }
        }
      }
      // going through all the pairs of of v[1] and a fixed v[0] for merging, if merge happens (1) each v[1] is erased (2) v[0] is updated (recurring until exit loop) (3) reset the combination again
      if (merge)
        v[0] = vertices->begin() - 1;  // (3) reset the combination if a valid merge happens
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Put the output.
  //////////////////////////////////////////////////////////////////////

  //Save the vertices
  iEvent.emplace(putToken_, std::move(*vertices));
  fillerZ.insert(seed_track_handle, shiftZVec.begin(), shiftZVec.end());
  fillerZ.fill();
  iEvent.put(std::move(outMapZ), "vtxZShift");
  filler3D.insert(seed_track_handle, shift3DVec.begin(), shift3DVec.end());
  filler3D.fill();
  iEvent.put(std::move(outMap3D), "vtx3DShift");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void ScoutingVertexer::beginStream(edm::StreamID) {
  // please remove this method if not needed
  edm::Service<TFileService> fs;
  h_dxyErr_weighted_sum_barrel =
      fs->make<TH1D>("dxyErr_weighted_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszErr_weighted_sum_barrel =
      fs->make<TH1D>("dszErr_weighted_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszdxyCov_weighted_sum_barrel =
      fs->make<TH1D>("dszdxyCov_weighted_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dxyErr_weighted_sq_sum_barrel =
      fs->make<TH1D>("dxyErr_weighted_sq_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszErr_weighted_sq_sum_barrel =
      fs->make<TH1D>("dszErr_weighted_sq_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszdxyCov_weighted_sq_sum_barrel =
      fs->make<TH1D>("dszdxyCov_weighted_sq_sum_barrel", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_weight_sum_barrel = fs->make<TH1D>("weight_sum_barrel", ";Track p_{T}; Weight Sum", 200, 0, 200);
  h_weight_sq_sum_barrel = fs->make<TH1D>("weight_sq_sum_barrel", ";Track p_{T}; Weight Squared Sum", 200, 0, 200);

  h_dxyErr_weighted_sum_disk = fs->make<TH1D>("dxyErr_weighted_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszErr_weighted_sum_disk = fs->make<TH1D>("dszErr_weighted_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszdxyCov_weighted_sum_disk =
      fs->make<TH1D>("dszdxyCov_weighted_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dxyErr_weighted_sq_sum_disk =
      fs->make<TH1D>("dxyErr_weighted_sq_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszErr_weighted_sq_sum_disk =
      fs->make<TH1D>("dszErr_weighted_sq_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_dszdxyCov_weighted_sq_sum_disk =
      fs->make<TH1D>("dszdxyCov_weighted_sq_sum_disk", ";Track p_{T}; Weighted Sum", 200, 0, 200);
  h_weight_sum_disk = fs->make<TH1D>("weight_sum_disk", ";Track p_{T}; Weight Sum", 200, 0, 200);
  h_weight_sq_sum_disk = fs->make<TH1D>("weight_sq_sum_disk", ";Track p_{T}; Weight Squared Sum", 200, 0, 200);
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void ScoutingVertexer::endStream() {
  // please remove this method if not needed
  h_dxyErr_weighted_sum_barrel->Draw();
  h_dxyErr_weighted_sum_barrel->Write();

  h_dszErr_weighted_sum_barrel->Draw();
  h_dszErr_weighted_sum_barrel->Write();

  h_dszdxyCov_weighted_sum_barrel->Draw();
  h_dszdxyCov_weighted_sum_barrel->Write();

  h_dxyErr_weighted_sq_sum_barrel->Draw();
  h_dxyErr_weighted_sq_sum_barrel->Write();

  h_dszErr_weighted_sq_sum_barrel->Draw();
  h_dszErr_weighted_sq_sum_barrel->Write();

  h_dszdxyCov_weighted_sq_sum_barrel->Draw();
  h_dszdxyCov_weighted_sq_sum_barrel->Write();

  h_weight_sum_barrel->Draw();
  h_weight_sum_barrel->Write();

  h_weight_sq_sum_barrel->Draw();
  h_weight_sq_sum_barrel->Write();

  h_dxyErr_weighted_sum_disk->Draw();
  h_dxyErr_weighted_sum_disk->Write();

  h_dszErr_weighted_sum_disk->Draw();
  h_dszErr_weighted_sum_disk->Write();

  h_dszdxyCov_weighted_sum_disk->Draw();
  h_dszdxyCov_weighted_sum_disk->Write();

  h_dxyErr_weighted_sq_sum_disk->Draw();
  h_dxyErr_weighted_sq_sum_disk->Write();

  h_dszErr_weighted_sq_sum_disk->Draw();
  h_dszErr_weighted_sq_sum_disk->Write();

  h_dszdxyCov_weighted_sq_sum_disk->Draw();
  h_dszdxyCov_weighted_sq_sum_disk->Write();

  h_weight_sum_disk->Draw();
  h_weight_sum_disk->Write();

  h_weight_sq_sum_disk->Draw();
  h_weight_sq_sum_disk->Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
ScoutingVertexer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
ScoutingVertexer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ScoutingVertexer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ScoutingVertexer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

/* ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingVertexer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingVertexer);

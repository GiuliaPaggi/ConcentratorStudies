#define Analiser_cxx
#include "Analiser.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <TVector.h>

#include <iostream>
#include <vector>

#include "Cluster.h"
#include "TriggerPrimitive.h"
#include "Segment.h"
#include "Digi.h"

// CB not optimal, but readable
const std::vector<int> WHEELS{-2, -1, 0, 1, 2};
const std::vector<int> SECTORS{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
const std::vector<int> STATIONS{1, 2, 3, 4};

double T_MIN{-9800};
double T_MAX{-9200};
double BX_MIN{-392};
double BX_MAX{-368};

const std::array<double, 4> MB{402.2, 490.5, 597.5, 700.0};
const int LOW_QUAL_CUT{0};
const int HIGH_QUAL_CUT{5};
const double PSI_CUT{TMath::Pi() / 6.0};
double PHI_CUT{0.02};
const double PHI_CUT_2{0.01};
const double T0_CUT{12.5};
const double X_CUT{5.0};
const double DIGI_CUT{10.0};
const int CORRECT_BX{-380};

std::vector<Cluster> buildClusters(std::vector<TriggerPrimitive> &tps, std::vector<Segment> &seg, std::vector<Digi> &d, double x_cut, double digi_cut) {
  std::vector<Cluster> clusters;
  std::vector<TriggerPrimitive> tpsToCluster = tps;
  std::vector<Segment> segToCluster = seg;
  std::vector<Digi> digiToCluster = d;

  for (const auto wh : WHEELS) {
    for (const auto sec : SECTORS) {
      for (const auto st : STATIONS) {
          while (true){   
          Cluster cluster{tpsToCluster, segToCluster, digiToCluster, x_cut, digi_cut, wh, sec, st};   
          if (cluster.bestTPQuality() > -1 || cluster.bestSegPhiHits() > -1 || cluster.WhichSL()) {
            clusters.push_back(cluster);  // CB can be improved 
          }
          else break;
        }
      }
    }
  }

  return clusters;
};

template<typename T> T getXY(TClonesArray * arr, int x, int y) 
{ 
  return static_cast<T>((*((TVectorT<float> *)(arr->At(x))))[y]); 
};

void Analiser::DefinePlot(){

  m_plots["N_DigiPerCluster"] = new TH1D("N_DigiPerCluster", "N_DigiPerCluster", 40, 0, 40); 
  m_plots["DigiSL"] = new TH1D("DigiSL", "DigiSL", 6, -0.5, 5.5); 
  m_plots["t0_LowQuality"] = new TH1D("t0_LowQuality", "t0_LowQuality; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["t0_Selected"] = new TH1D("t0_Selected", "t0_Selected; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["LowQ_matched"] = new TH1D("LowQ_matched", "LowQ_matched; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["BX_LowQuality"] = new TH1D("BX_LowQuality", "BX_LowQuality; BX; Entries", 24, BX_MIN, BX_MAX);
  m_plots["BX_LowQ_more1HQ"] =  new TH1D("BX_LowQ_more1HQ", "BX_LowQ_more1HQ,BX;Entries", 24, BX_MIN, BX_MAX);
  m_plots["LowQ_more1HQ_Phi"] = new TH1D("LowQ_more1HQ_Phi", "LowQ_more1HQ_Phi", 100, T_MIN, T_MAX);
  m_plots["BX_LowQ_matched"] =  new TH1D("BX_LowQ_matched", "BX_LowQ_matched,BX;Entries", 24, BX_MIN, BX_MAX);

  for (auto CLtype : tags){
    m_plots[Form("%s_N_Ghost", CLtype.c_str())] = new TH1I(Form("%s_N_Ghost", CLtype.c_str()), Form("%s_N_Ghost", CLtype.c_str()), 20, 0, 20);
    m_plots[Form("%s_Q_Best", CLtype.c_str())] = new TH1I(Form("%s_Q_Best", CLtype.c_str()), Form("%s_Q_Best", CLtype.c_str()), 10, 0, 10);
    m_plots[Form("%s_Q_Ghost", CLtype.c_str())] = new TH1I(Form("%s_Q_Ghost", CLtype.c_str()), Form("%s_Q_Ghost", CLtype.c_str()), 10, 0, 10);
    m_plots[Form("%s_OoTGhosts", CLtype.c_str())] = new TH1D(Form("%s_OoTGhosts", CLtype.c_str()), Form("%s_OoTGhosts", CLtype.c_str()), 20, 0, 20);
    m_plots[Form("%s_ITGhosts", CLtype.c_str())] = new TH1D(Form("%s_ITGhosts", CLtype.c_str()), Form("%s_ITGhosts", CLtype.c_str()), 20, 0, 20);
    m_plots[Form("%s_BX_ITGhosts", CLtype.c_str())] = new TH1I(Form("%s_BX_ITGhosts", CLtype.c_str()), Form("%s_BX_ITGhosts", CLtype.c_str()), 24, BX_MIN, BX_MAX);
    m_plots[Form("%s_BX_OoTGhosts", CLtype.c_str())] = new TH1I(Form("%s_BX_OoTGhosts", CLtype.c_str()), Form("%s_BX_OoTGhosts", CLtype.c_str()), 25, BX_MIN, BX_MAX);
    m_plots[Form("%s_Res_ITGhosts", CLtype.c_str())] = new TH1D(Form("%s_Res_ITGhosts", CLtype.c_str()), Form("%s_Res_ITGhosts", CLtype.c_str()), 111, -5.5, 5.5);
    m_plots[Form("%s_Res_OoTGhosts", CLtype.c_str())] = new TH1D(Form("%s_Res_OoTGhosts", CLtype.c_str()), Form("%s_Res_OoTGhosts", CLtype.c_str()), 111, -5.5, 5.5);
    m_plots[Form("%s_Res_SegMatched", CLtype.c_str())]= new TH1D(Form("%s_Res_SegMatched", CLtype.c_str()), Form("%s_Res_SegMatched; cluster.bestTP().xLoc - seg.xLoc; Entries", CLtype.c_str()), 101, -10, 10);
    
    m_2Dplots[Form("%s_Q_OoTGhosts", CLtype.c_str())] = new TH2D(Form("%s_Q_OoTGhosts", CLtype.c_str()), Form("%s_Q_OoTGhosts;High Quality;Out of time Ghost Quality", CLtype.c_str()), 10, 0, 10, 10, 0, 10);
    m_2Dplots[Form("%s_Q_ITGhosts", CLtype.c_str())] = new TH2D(Form("%s_Q_ITGhosts", CLtype.c_str()), Form("%s_Q_ITGhosts;High Quality;In time Ghost Quality", CLtype.c_str()), 10, 0, 10, 10, 0, 10);
    

    for (const auto st : STATIONS) {
      m_plots[Form("%s_DigiTime_st%d", CLtype.c_str(), st)]   = new TH1D( Form("%s_DigiTime_st%d", CLtype.c_str(), st), Form("%s_DigiTime_st%d; mean digi cluster time (ns); entries", CLtype.c_str(), st), 300, 400, 1000 );
      m_plots[Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st)] = new TH1D(Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st), Form("%s_x_LowBestQ_st%d; xLoc; Entries",CLtype.c_str(), st ), 100, -220, 220);
      
      m_2Dplots[Form("%s_N_Cluster_st%d", CLtype.c_str(), st)] = new TH2I(Form("%s_N_Cluster_st%d", CLtype.c_str(), st), Form("%s_N_Cluster_st%d", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5); 
      m_2Dplots[Form("%s_N_Digi_st%d", CLtype.c_str(), st)] = new TH2I(Form("%s_N_Digi_st%d", CLtype.c_str(), st), Form("%s_N_Digi_st%d", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", CLtype.c_str(), st)] = new TH2D (Form("%s_TPs_LocalDirectionvsPosition_st%d", CLtype.c_str(), st), Form("AllTPs_LocalDirectionvsPosition_st%d", st ), 200, -200, 200, 11, -1, 1 );
      
      m_effs[Form("%s_ClusterEfficiency_st%d", CLtype.c_str(), st)] = new TEfficiency(Form("%s_ClusterEfficiency_st%d", CLtype.c_str(), st), Form("%s_ClusterEfficiency_st%d; sector; wheel", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_effs[Form("%s_Eff_SegMatch_st%d", CLtype.c_str(), st)] = new TEfficiency(Form("%s_Eff_SegMatch_st%d", CLtype.c_str(), st), Form("%s_Eff_SegMatch_st%d; sector; wheel", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_effs[Form("%s_Eff_DigiMatch_st%d", CLtype.c_str(), st)] = new TEfficiency(Form("%s_Eff_DigiMatch_st%d", CLtype.c_str(), st), Form("%s_Eff_DigiMatch_st%d; sector; wheel", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      
      for (const auto wh : WHEELS){
        m_plots[Form("%s_Digi_residual_st%d_wh%d", CLtype.c_str(), st, wh)] = new TH1D( Form("%s_Digi_residual_st%d_wh%d", CLtype.c_str(), st, wh), Form("%s_Digi_residual_st%d_wh%d; cluster.bestSeg.xLoc - digi.xLoc; entries", CLtype.c_str(), st, wh), 51, -10, 10 );
      }

    }
  }

  for (const auto st : STATIONS) {
    m_plots[Form("Res_MuMatched_st%d", st)] = new TH1D(Form("Res_MuMatched_st%d", st), Form("Res_MuMatched; cluster.bestTP().xLoc - mu.xLoc; Entries"), 101, -10, 10);
    
    m_2Dplots[Form("N_MuMatch_st%d", st)] = new TH2I(Form("N_MuMatch_st%d", st), Form("N_MuMatch_st%d; sector ; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);

    m_effs[Form("Eff_MuMatch_st%d", st)] = new TEfficiency(Form("Eff_MuMatch_st%d", st), Form("Eff_MuMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
  }
}

void Analiser::ClusterAnalisis(std::vector<Cluster> CLtoAnalize, string CLtype,  std::vector<Segment> Segments){
    

    for (auto const &cluster : CLtoAnalize) {
      auto wh{cluster.wheel};
      auto sec{cluster.sector};
      if (sec == 13 ) sec = 4;
      if (sec == 14 ) sec = 10;
      auto st{cluster.station};

      m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", CLtype.c_str(), st)]->Fill(cluster.bestTP().xLoc, cluster.bestTP().psi);

      // ########## Study TP ghost distribution #############
      ++nClusters;
      ooTHQCount += cluster.ootCountIf([=](TriggerPrimitive const & tp) { return tp.quality > HIGH_QUAL_CUT; });

      int bestQ = cluster.bestTPQuality();
      int ootSize{cluster.ootSize()};
      int itSize{cluster.itSize()};

      m_plots[Form("%s_Q_Best", CLtype.c_str())]->Fill(bestQ);
      m_plots[Form("%s_ITGhosts", CLtype.c_str())]->Fill(itSize);
      m_plots[Form("%s_OoTGhosts", CLtype.c_str())]->Fill(ootSize);
      m_plots[Form("%s_N_Ghost", CLtype.c_str())]->Fill(ootSize + itSize);

      if (bestQ == 1) m_plots[Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st)]->Fill(cluster.bestTP().xLoc);
      
      m_2Dplots[Form("%s_N_Cluster_st%d", CLtype.c_str(), st)]->Fill(sec, wh);
      
      if (cluster.hasGhosts()) {
        ++nClustersGhosts;

        for (const auto &ghost : cluster.ootGhosts()) {
          m_plots[Form("%s_BX_OoTGhosts", CLtype.c_str())]  ->Fill(ghost.BX);
          m_plots[Form("%s_Res_OoTGhosts", CLtype.c_str())] ->Fill(ghost.xLoc - cluster.bestTP().xLoc);
          m_2Dplots[Form("%s_Q_OoTGhosts", CLtype.c_str())] ->Fill(bestQ, ghost.quality);
          m_plots[Form("%s_Q_Ghost", CLtype.c_str())]  ->Fill(ghost.quality);
        }

        for (const auto &ghost : cluster.itGhosts()) {
          m_plots[Form("%s_BX_ITGhosts", CLtype.c_str())]->Fill(ghost.BX);
          m_plots[Form("%s_Res_ITGhosts", CLtype.c_str())]->Fill(ghost.xLoc - cluster.bestTP().xLoc);
          m_2Dplots[Form("%s_Q_ITGhosts", CLtype.c_str())]->Fill(bestQ, ghost.quality);
          m_plots[Form("%s_Q_Ghost", CLtype.c_str())]->Fill(ghost.quality);
        }
      }

      if (cluster.muMatched && cluster.bestSeg().nPhiHits >= 4 ) {
        bool efficient = std::abs(cluster.bestTP().BX - CORRECT_BX) < 1;
        m_effs[Form("%s_ClusterEfficiency_st%d", CLtype.c_str(), st)]->Fill(efficient , sec, wh); 
      }

      // ########## Study segment matching #############
      m_effs[Form("%s_Eff_SegMatch_st%d", CLtype.c_str(), st)]->Fill(cluster.segMatched, sec, wh);
      for (const auto s : Segments)
      if (cluster.foundTP && s.wheel == wh && s.sector == sec && s.station == st){
        m_plots[Form("%s_Res_SegMatched", CLtype.c_str())]->Fill( cluster.bestTP().xLoc - s.xLoc );
      }
      
      // ########## Study digi clusters #############
      m_effs[Form("%s_Eff_DigiMatch_st%d", CLtype.c_str(), st)]->Fill(cluster.digiMatched, sec, wh);
      
      if (cluster.digiMatched){
        if (cluster.MeanDigiTime() > 0) m_plots[Form("%s_DigiTime_st%d", CLtype.c_str(), st)]->Fill(cluster.MeanDigiTime());
        for (auto &digi : cluster.matchedDigi()){
          m_plots[Form("%s_Digi_residual_st%d_wh%d", CLtype.c_str(), st, wh)]->Fill( cluster.bestSeg().xLoc - digi.xLoc );
        }
      }
     
      m_plots["N_DigiPerCluster"]->Fill(cluster.GetNDigi());
      m_plots["DigiSL"]->Fill(cluster.WhichSL());

      if (cluster.foundDigi) m_2Dplots[Form("%s_N_Digi_st%d", CLtype.c_str(), st)]->Fill(sec, wh);
    }

}


void Analiser::Loop() {

  TFile *file = new TFile("./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root");   
  //./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root
  //  ./VtxSmeared/DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25_VtxSmeared.root

  TFile outputFile("prompt/outputFile.root","RECREATE");
  outputFile.cd();

  tags.push_back("PreFilter");
  tags.push_back("PostFilter");
 
// 4+2 3+2 qualities, not there if we use slice-test configuration for emulator
 
  DefinePlot();


  if (fChain == 0) return;

  Long64_t n_entries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < n_entries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // ########## LOOP ON EVENTS #############
    if (jentry % 100 == 0) cout << "Processing event: " << jentry << '\r' << flush;

    if (std::abs(gen_pdgId->at(0)) != 13 || std::abs(gen_eta->at(0)) > 0.8) continue;

    // ########## CREATE TPs std::vector #############
    std::vector<TriggerPrimitive> tps;

    for (std::size_t j = 0; j < ph2TpgPhiEmuAm_nTrigs; ++j) {
      tps.emplace_back(TriggerPrimitive{j, ph2TpgPhiEmuAm_wheel->at(j),
                                        ph2TpgPhiEmuAm_sector->at(j), ph2TpgPhiEmuAm_station->at(j),
                                        ph2TpgPhiEmuAm_quality->at(j), ph2TpgPhiEmuAm_phi->at(j),
                                        ph2TpgPhiEmuAm_phiB->at(j), ph2TpgPhiEmuAm_BX->at(j),
                                        ph2TpgPhiEmuAm_t0->at(j), ph2TpgPhiEmuAm_posLoc_x->at(j)});
    }

    // ########## CREATE Digis std::vector #############
    std::vector<Digi> digis;
    for (std::size_t i = 0; i < digi_nDigis; ++i ){
      digis.emplace_back(Digi(i, digi_wheel->at(i), digi_sector->at(i), 
                              digi_station->at(i), digi_superLayer->at(i), 
                              digi_layer->at(i), digi_wire->at(i), 
                              digi_time->at(i)));
    }

    //########## BUILD segments std::vector #############
    std::vector<Segment> segments;
    for (std::size_t j = 0; j < seg_nSegments; ++j) {
      segments.emplace_back(Segment(j, seg_station->at(j), seg_wheel->at(j), 
                                    seg_sector->at(j), seg_phi_nHits->at(j), 
                                    seg_posLoc_x->at(j)));
    }

    // ########## BUILD clusters std::vector #############
    auto clusters = buildClusters(tps, segments, digis, X_CUT, DIGI_CUT);

    // ########## BUILD cluster with phi matchin information #############
    // tag TPs that match using the extrapolation on a straight line 
        for (TriggerPrimitive &tp : tps) {
      if (tp.quality == 1) {
        m_plots["t0_LowQuality"]->Fill(tps.back().t0);
        m_plots["BX_LowQuality"]->Fill(tps.back().BX);
      }

      // select HQ tp to try and match in previous/following station with the expected value of phi from straight line extrapolation
      if (tp.quality > HIGH_QUAL_CUT && !tp.hasMatched) {
        // select HQ TPs which are not matched 
        m_plots["t0_Selected"]->Fill(tp.t0);
        for (TriggerPrimitive &other_tp : tps) {
          if (tp.index != other_tp.index && tp.Match(other_tp, PHI_CUT, T0_CUT)) {    
            if (other_tp.quality == 1) {
              m_plots["t0_Selected"]->Fill(other_tp.t0);
              m_plots["LowQ_matched"]->Fill(other_tp.t0);
              m_plots["BX_LowQ_matched"]->Fill(other_tp.BX);
            } else {
              m_plots["t0_Selected"]->Fill(other_tp.t0);
            }
          }
        }
      }
    }

      // repeat clustering with phi match information in tps
    auto PostMatch_clusters = buildClusters(tps, segments, digis, X_CUT, DIGI_CUT);


    // ########## ATTEMPT cluster - muon extrapolation matching #############
    
    for (int iMu = 0; iMu < mu_nMuons; ++iMu){
      for (int i = 0; i < mu_nMatches->at(iMu); ++i){
        int muTrkWheel = getXY<float>(mu_matches_wheel, iMu, i);
        int muTrkStation = getXY<float>(mu_matches_station, iMu, i);
        int muTrkSector = getXY<float>(mu_matches_sector, iMu, i);
        if (muTrkSector == 13 ) muTrkSector = 4;
        if (muTrkSector == 14 ) muTrkSector = 10; 
        double muTrkX = getXY<float>(mu_matches_x, iMu, i);
        double edgeX = getXY<float>(mu_matches_edgeX, iMu, i);
        double edgeY = getXY<float>(mu_matches_edgeY, iMu, i);

        for (auto &cluster: clusters){
          auto wh{cluster.wheel};
          auto sec{cluster.sector};
          if (sec == 13 ) sec = 4;
          if (sec == 14 ) sec = 10;
          auto st{cluster.station};
          cluster.MatchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
          cluster.MatchDigi(digis, DIGI_CUT);
          m_effs[Form("Eff_MuMatch_st%d", st)]->Fill(cluster.muMatched, sec, wh);
        if (cluster.muMatched) {
          m_plots[Form("Res_MuMatched_st%d", st)]->Fill(cluster.bestTP().xLoc - getXY<float>(mu_matches_x, cluster.muMatchedIndex[0], cluster.muMatchedIndex[1]));
          m_2Dplots[Form("N_MuMatch_st%d", st)]->Fill(sec, wh);
      }
        }

        for (auto &PMcluster : PostMatch_clusters) {
          PMcluster.MatchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
          PMcluster.MatchDigi(digis, DIGI_CUT);
        }
      }
    }

    


    // ########## RUN SOME ANALYSIS #############
    // ########## CLUSTER ANALYSIS #############
    ClusterAnalisis(clusters, tags[0], segments);
    ClusterAnalisis(PostMatch_clusters, tags[1], segments);

    // ########## PHI MATCHING TPs ANALYSIS #############
    for (TriggerPrimitive tp : tps) {
      if (tp.quality == 1 && tp.Matches.size() > 0) {
        m_plots["LowQ_more1HQ_Phi"]->Fill(tp.t0);
        m_plots["BX_LowQ_more1HQ"] ->Fill(tp.BX);
      }
    }


  }

  double ghostFraction = nClustersGhosts / nClusters;

  cout << " Ratio LQ/selected with phi=  " << m_plots["LowQ_matched"]->GetEntries() << "/ "
       << m_plots["t0_LowQuality"]->GetEntries() << " = "
       << m_plots["LowQ_matched"]->GetEntries() / m_plots["t0_LowQuality"]->GetEntries() << endl;
  cout << " Ratio LQ/matched with phi=  " << m_plots["LowQ_more1HQ_Phi"]->GetEntries() << "/ "
       << m_plots["t0_LowQuality"]->GetEntries() << " = "
       << m_plots["LowQ_more1HQ_Phi"]->GetEntries() / m_plots["t0_LowQuality"]->GetEntries() << endl;
  cout << " Fraction of clusters with ghost (" << nClustersGhosts << ") on total (" << nClusters
       << ") = " << ghostFraction << endl;
  cout << " HQ out of time clusters: " << ooTHQCount << endl;

  outputFile.Write();
  outputFile.Close();
}


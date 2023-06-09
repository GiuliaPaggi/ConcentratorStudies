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

  m_plots["DigiSL"] = new TH1D("DigiSL", "DigiSL", 6, -0.5, 5.5); 
  m_plots["t0_LowQuality"] = new TH1D("t0_LowQuality", "t0_LowQuality; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["BX_LowQuality"] = new TH1D("BX_LowQuality", "BX_LowQuality; BX; Entries", 24, BX_MIN, BX_MAX);
  m_plots["BX_LowQ_more1HQ"] =  new TH1D("BX_LowQ_more1HQ", "BX_LowQ_more1HQ; BX; Entries", 24, BX_MIN, BX_MAX);

  
  for (auto CLtype : tags){
    m_plots[Form("%s_N_DigiPerCluster", CLtype.c_str())] = new TH1D(Form("%s_N_DigiPerCluster", CLtype.c_str()), Form("%s_N_DigiPerCluster; # digi in cluster; Entries", CLtype.c_str()), 50, 0, 50); 
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
    //m_plots[Form("%s_LowQ_more1HQ_Phi", CLtype.c_str())] = new TH1D(Form("%s_LowQ_more1HQ_Phi", CLtype.c_str()), Form("%s_LowQ_more1HQ_Phi", CLtype.c_str()), 100, T_MIN, T_MAX);
    m_plots[Form("%s_LowQ_matched", CLtype.c_str())] = new TH1D(Form("%s_LowQ_matched", CLtype.c_str()), Form("%s_LowQ_matched; t0 (ns); Entries", CLtype.c_str()), 100, T_MIN, T_MAX); 
    m_plots[Form("%s_BX_LowQ_matched", CLtype.c_str())] =  new TH1D(Form("%s_BX_LowQ_matched", CLtype.c_str()), Form("%s_BX_LowQ_matched,BX;Entries", CLtype.c_str()), 24, BX_MIN, BX_MAX);
    m_plots[Form("%s_t0_Selected", CLtype.c_str())] = new TH1D(Form("%s_t0_Selected", CLtype.c_str()), Form("%s_t0_Selected; t0 (ns); Entries", CLtype.c_str()), 100, T_MIN, T_MAX);
    m_plots[Form("%s_RemovedQuality", CLtype.c_str())] = new TH1D(Form("%s_RemovedQuality", CLtype.c_str()), Form("%s_RemovedQuality; Quality; Entries", CLtype.c_str()), 10, -0.5, 9.5);
    m_plots[Form("%s_BestQ_RemovedCluster", CLtype.c_str())] = new TH1D(Form("%s_BestQ_RemovedCluster", CLtype.c_str()), Form("%s_BestQ_RemovedCluster; Quality; Entries", CLtype.c_str()), 10, -0.5, 9.5);
    m_plots[Form("%s_ResFilteredCluster", CLtype.c_str())] = new TH1D(Form("%s_ResFilteredCluster", CLtype.c_str()), Form("%s_ResFilteredCluster; FilteredCl - OriginalCl; Entries ", CLtype.c_str()), 101, -50, 50);
    m_plots[Form("%s_ClusterSize", CLtype.c_str())] = new TH1I(Form("%s_ClusterSize", CLtype.c_str()), Form("%s_ClusterSize; # TPs in cluster; Entries", CLtype.c_str()), 21, -.5, 20.5 );

    m_2Dplots[Form("%s_Q_OoTGhosts", CLtype.c_str())] = new TH2D(Form("%s_Q_OoTGhosts", CLtype.c_str()), Form("%s_Q_OoTGhosts;High Quality;Out of time Ghost Quality", CLtype.c_str()), 10, 0, 10, 10, 0, 10);
    m_2Dplots[Form("%s_Q_ITGhosts", CLtype.c_str())] = new TH2D(Form("%s_Q_ITGhosts", CLtype.c_str()), Form("%s_Q_ITGhosts;High Quality;In time Ghost Quality", CLtype.c_str()), 10, 0, 10, 10, 0, 10);
    
    m_counters[Form("%s_nClustersGhosts", CLtype.c_str())] = 0;
    m_counters[Form("%s_ooTHQCount", CLtype.c_str())] = 0;
    m_counters[Form("%s_nClusters", CLtype.c_str())] = 0;
    m_counters[Form("%s_nTPs", CLtype.c_str())] = 0;

    for (const auto st : STATIONS) {
      m_plots[Form("%s_DigiTime_st%d", CLtype.c_str(), st)]   = new TH1D( Form("%s_DigiTime_st%d", CLtype.c_str(), st), Form("%s_DigiTime_st%d; mean digi cluster time (ns); entries", CLtype.c_str(), st), 300, 400, 1000 );
      m_plots[Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st)] = new TH1D(Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st), Form("%s_x_LowBestQ_st%d; xLoc; Entries",CLtype.c_str(), st ), 100, -220, 220);
      m_plots[Form("%s_meanXLoc_RemovedClusters_st%d",CLtype.c_str(), st)] = new TH1D(Form("%s_meanXLoc_RemovedClusters_st%d",CLtype.c_str(), st), Form("%s_meanXLoc_RemovedClusters_st%d;xLoc;Entries",CLtype.c_str(), st), 101, -220, 220);
      m_2Dplots[Form("%s_N_Cluster_st%d", CLtype.c_str(), st)] = new TH2I(Form("%s_N_Cluster_st%d", CLtype.c_str(), st), Form("%s_N_Cluster_st%d", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5); 
      m_2Dplots[Form("%s_N_Digi_st%d", CLtype.c_str(), st)] = new TH2I(Form("%s_N_Digi_st%d", CLtype.c_str(), st), Form("%s_N_Digi_st%d", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", CLtype.c_str(), st)] = new TH2D(Form("%s_TPs_LocalDirectionvsPosition_st%d", CLtype.c_str(), st), Form("AllTPs_LocalDirectionvsPosition_st%d", st ), 200, -200, 200, 11, -1, 1 );
      m_2Dplots[Form("%s_NCLusterRemoved_st%d", CLtype.c_str(), st)] = new TH2I(Form("%s_NCLusterRemoved_st%d", CLtype.c_str(), st), Form("%s_NCLusterRemoved_st%d; sector; wheel", CLtype.c_str(), st), 14, -0.5, 13.5, 7, -3.5, 3.5);

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
    m_2Dplots[Form("LQnotinHQ_st%d", st)] = new TH2I(Form("LQnotinHQ_st%d", st), Form("LQnotinHQ_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    m_2Dplots[Form("HQnotinLQ_st%d", st)] = new TH2I(Form("HQnotinLQ_st%d", st), Form("HQnotinLQ_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);

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
      m_plots[Form("%s_ClusterSize", CLtype.c_str())]->Fill(cluster.tpClusterSize()); 
      // ########## Study TP ghost distribution #############
      //++nClusters;
      ++m_counters[Form("%s_nClusters", CLtype.c_str())];
      //ooTHQCount += cluster.ootCountIf([=](TriggerPrimitive const & tp) { return tp.quality > HIGH_QUAL_CUT; });
      m_counters[Form("%s_ooTHQCount", CLtype.c_str())] += cluster.ootCountIf([=](TriggerPrimitive const & tp) { return tp.quality > HIGH_QUAL_CUT; });

      int bestQ = cluster.bestTPQuality();
      int ootSize{cluster.ootSize()};
      int itSize{cluster.itSize()};

      m_plots[Form("%s_Q_Best", CLtype.c_str())]->Fill(bestQ);
      m_plots[Form("%s_ITGhosts", CLtype.c_str())]->Fill(itSize);
      m_plots[Form("%s_OoTGhosts", CLtype.c_str())]->Fill(ootSize);
      m_plots[Form("%s_N_Ghost", CLtype.c_str())]->Fill(ootSize + itSize);
      m_counters[Form("%s_nTPs", CLtype.c_str())] += cluster.tpClusterSize();


      if (bestQ == 1) m_plots[Form("%s_x_LowBestQ_st%d", CLtype.c_str(), st)]->Fill(cluster.bestTP().xLoc);
      
      m_2Dplots[Form("%s_N_Cluster_st%d", CLtype.c_str(), st)]->Fill(sec, wh);
      
      if (cluster.hasGhosts()) {
        //++nClustersGhosts;
        ++m_counters[Form("%s_nClustersGhosts", CLtype.c_str())];

        for (const auto &ghost : cluster.ootGhosts()) {
          m_plots[Form("%s_BX_OoTGhosts", CLtype.c_str())]  ->Fill(ghost.BX);
          m_plots[Form("%s_Res_OoTGhosts", CLtype.c_str())] ->Fill(ghost.xLoc - cluster.bestTP().xLoc);
          m_2Dplots[Form("%s_Q_OoTGhosts", CLtype.c_str())] ->Fill(bestQ, ghost.quality);
          m_plots[Form("%s_Q_Ghost", CLtype.c_str())]  ->Fill(ghost.quality);

          if (ghost.quality == 1) {
            m_plots[Form("%s_t0_Selected", CLtype.c_str())]->Fill(ghost.t0);
            m_plots[Form("%s_LowQ_matched", CLtype.c_str())]->Fill(ghost.t0);
            m_plots[Form("%s_BX_LowQ_matched", CLtype.c_str())]->Fill(ghost.BX);
          }
        }

        for (const auto &ghost : cluster.itGhosts()) {
          m_plots[Form("%s_BX_ITGhosts", CLtype.c_str())]->Fill(ghost.BX);
          m_plots[Form("%s_Res_ITGhosts", CLtype.c_str())]->Fill(ghost.xLoc - cluster.bestTP().xLoc);
          m_2Dplots[Form("%s_Q_ITGhosts", CLtype.c_str())]->Fill(bestQ, ghost.quality);
          m_plots[Form("%s_Q_Ghost", CLtype.c_str())]->Fill(ghost.quality);

          if (ghost.quality == 1) {
            m_plots[Form("%s_t0_Selected", CLtype.c_str())]->Fill(ghost.t0);
            m_plots[Form("%s_LowQ_matched", CLtype.c_str())]->Fill(ghost.t0);
            m_plots[Form("%s_BX_LowQ_matched", CLtype.c_str())]->Fill(ghost.BX);
          }
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
      //if (cluster.GetNDigi() > 0) std::cout << " N digi in cluster " << cluster.GetNDigi() << std::endl;
      m_effs[Form("%s_Eff_DigiMatch_st%d", CLtype.c_str(), st)]->Fill(cluster.digiMatched, sec, wh);
      
      if (cluster.digiMatched){
        if (cluster.MeanDigiTime() > 0) m_plots[Form("%s_DigiTime_st%d", CLtype.c_str(), st)]->Fill(cluster.MeanDigiTime());
        for (auto &digi : cluster.matchedDigi()){
          m_plots[Form("%s_Digi_residual_st%d_wh%d", CLtype.c_str(), st, wh)]->Fill( cluster.bestSeg().xLoc - digi.xLoc );
        }
      }
     
      m_plots[Form("%s_N_DigiPerCluster", CLtype.c_str())]->Fill(cluster.GetNDigi());
      m_plots["DigiSL"]->Fill(cluster.WhichSL());

      if (cluster.foundDigi) m_2Dplots[Form("%s_N_Digi_st%d", CLtype.c_str(), st)]->Fill(sec, wh);
    }


}

std::vector<Cluster> MissingClusters(std::vector<Cluster> FirstCLVector, std::vector<Cluster> SecondCLvector) {
  // FirstCL is the smaller set (filtered) , SecondCL is the bigger one (all) 
  std::vector<Cluster> missingClusters;

  for (auto FirstCluster : FirstCLVector) {
    if (!FirstCluster.foundTP) continue;
    bool found = false;
    auto wh{FirstCluster.wheel};
    auto sec{FirstCluster.sector};
    auto st{FirstCluster.station};

    for (auto SecondCluster : SecondCLvector) {
      if (!SecondCluster.foundTP) continue;
      if ( std::abs(FirstCluster.bestTP().xLoc - SecondCluster .bestTP().xLoc) < 3.5 ) {
        found = true;
      }
    }

    if (found == false) missingClusters.push_back(FirstCluster);

  }
  return missingClusters;
}


void Analiser::Loop() {

  TFile *file = new TFile("./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root");   
  //./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root
  //  ./VtxSmeared/DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25_VtxSmeared.root

  TFile outputFile("prompt/outputFile.root","RECREATE");
  outputFile.cd();

  tags.push_back("PreFilter");
  tags.push_back("LQFilter");
  tags.push_back("HQFilter");
 
// 4+2 3+2 qualities, not there if we use slice-test configuration for emulator
 
  DefinePlot();

  double removedHQFilter{0};
  double removedLQFilter{0}; 

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
    //std::cout << " n Digis: " << digis.size() << std::endl;
    //########## BUILD segments std::vector #############
    std::vector<Segment> segments;
    for (std::size_t j = 0; j < seg_nSegments; ++j) {
      segments.emplace_back(Segment(j, seg_station->at(j), seg_wheel->at(j), 
                                    seg_sector->at(j), seg_phi_nHits->at(j), 
                                    seg_posLoc_x->at(j)));
    }

    std::vector<Cluster> MatchFromLQ_clusters;
    std::vector<Cluster> MatchFromHQ_clusters;

    // ########## BUILD clusters std::vector #############
    auto clusters = buildClusters(tps, segments, digis, X_CUT, DIGI_CUT);


    // ########## BUILD cluster with phi matchin information #############
    // tag TPs that match using the extrapolation on a straight line 
  
    std::vector<TriggerPrimitive> MatchFromLQ_tps = tps;
    std::vector<TriggerPrimitive> MatchFromHQ_tps = tps;

    //extrapolate from LQ TP
    for (TriggerPrimitive &tp : MatchFromLQ_tps) {

      if (!tp.hasMatched) {

        for (TriggerPrimitive &other_tp : MatchFromLQ_tps) {
          if (tp.index != other_tp.index && tp.MatchFromLQ(other_tp, PHI_CUT, T0_CUT)) {    

          }
        }
        // If after the loop it did not match any HQ i remove it from the sample
        //if (!tp.hasMatched) MatchFromLQ_tps.erase(std::remove(MatchFromLQ_tps.begin(), MatchFromLQ_tps.end(), tp), MatchFromLQ_tps.end());
      }
    }
    MatchFromLQ_tps.erase(std::remove_if(MatchFromLQ_tps.begin(), MatchFromLQ_tps.end(), [](auto &tp){ return !tp.hasMatched && tp.quality == 1;}), MatchFromLQ_tps.end());
    // repeat clustering after filtering
    MatchFromLQ_clusters = buildClusters(MatchFromLQ_tps, segments, digis, X_CUT, DIGI_CUT);



    //extrapolate from HQ TP
    for (TriggerPrimitive &tp : MatchFromHQ_tps) {
      if (!tp.hasMatched) {

        for (TriggerPrimitive &other_tp : MatchFromHQ_tps) {
          
          if (tp.index != other_tp.index && tp.MatchFromHQ(other_tp, PHI_CUT, T0_CUT)) {    

          }
        }
      }
    }
    MatchFromHQ_tps.erase(std::remove_if(MatchFromHQ_tps.begin(), MatchFromHQ_tps.end(), [](auto &tp){ return !tp.hasMatched && tp.quality == 1;}), MatchFromHQ_tps.end());
    // repeat clustering after filtering
    MatchFromHQ_clusters = buildClusters(MatchFromHQ_tps, segments, digis, X_CUT, DIGI_CUT);

    std::vector<TriggerPrimitive> TPRemovedfromLQ = tps;
    for (auto tp : MatchFromLQ_tps) TPRemovedfromLQ.erase(std::remove(TPRemovedfromLQ.begin(), TPRemovedfromLQ.end(), tp), TPRemovedfromLQ.end());
    std::vector<TriggerPrimitive> TPRemovedfromHQ = tps;
    for (auto tp : MatchFromHQ_tps) TPRemovedfromHQ.erase(std::remove(TPRemovedfromHQ.begin(), TPRemovedfromHQ.end(), tp), TPRemovedfromHQ.end());

    for (auto tp : TPRemovedfromLQ) m_plots[Form("%s_RemovedQuality", "LQFilter")]->Fill(tp.quality);
    for (auto tp : TPRemovedfromHQ) m_plots[Form("%s_RemovedQuality", "HQFilter")]->Fill(tp.quality);

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

        for (auto &PMcluster : MatchFromLQ_clusters) {
          PMcluster.MatchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
          PMcluster.MatchDigi(digis, DIGI_CUT);
        }
        for (auto &PMcluster : MatchFromHQ_clusters) {
          PMcluster.MatchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
          PMcluster.MatchDigi(digis, DIGI_CUT);
        }
      }
    }

    // ########## RUN SOME ANALYSIS #############
    // ########## CLUSTER ANALYSIS #############
    ClusterAnalisis(clusters, tags[0], segments);
    ClusterAnalisis(MatchFromLQ_clusters, tags[1], segments);
    ClusterAnalisis(MatchFromHQ_clusters, tags[2], segments);

    for (auto FirstCluster : MatchFromLQ_clusters){
      for (auto SecondCluster : clusters){
        if (SecondCluster.sector != FirstCluster.sector || SecondCluster.wheel != FirstCluster.wheel || SecondCluster.station != FirstCluster.station) continue;
        double res = FirstCluster.bestTP().xLoc - SecondCluster .bestTP().xLoc;
        if (res != 0 ) m_plots[Form("%s_ResFilteredCluster", "LQFilter")]->Fill(res);
      }
    }

    for (auto FirstCluster : MatchFromHQ_clusters){
      for (auto SecondCluster : clusters){
        if (SecondCluster.sector != FirstCluster.sector || SecondCluster.wheel != FirstCluster.wheel || SecondCluster.station != FirstCluster.station) continue;
        double res = FirstCluster.bestTP().xLoc  - SecondCluster .bestTP().xLoc;
        if (res != 0) m_plots[Form("%s_ResFilteredCluster", "HQFilter")]->Fill(res);
      }
    }
    
    std::vector<Cluster> ClusterCut_LQFilter = MissingClusters(MatchFromLQ_clusters, clusters);
    removedLQFilter += ClusterCut_LQFilter.size();
    if (ClusterCut_LQFilter.size() > 0) {
      double meanxLoc{0};
      for (auto Cl: ClusterCut_LQFilter) {
        meanxLoc+=Cl.bestTP().xLoc;
        m_2Dplots[Form("%s_NCLusterRemoved_st%d", "LQFilter", Cl.station)]->Fill(Cl.sector, Cl.wheel);
        m_plots[Form("%s_BestQ_RemovedCluster", "LQFilter")]->Fill(Cl.bestTP().quality);
      }
      m_plots[Form("%s_meanXLoc_RemovedClusters_st%d", "LQFilter" , ClusterCut_LQFilter[0].station)]->Fill(meanxLoc/ClusterCut_LQFilter.size());
    }

    std::vector<Cluster> ClusterCut_HQFilter = MissingClusters(MatchFromHQ_clusters, clusters);
    removedHQFilter += ClusterCut_HQFilter.size();

    if (ClusterCut_HQFilter.size() > 0) {
      double meanxLoc{0};
      for (auto Cl: ClusterCut_HQFilter) {
        meanxLoc+=Cl.bestTP().xLoc;
        m_2Dplots[Form("%s_NCLusterRemoved_st%d", "HQFilter", Cl.station)]->Fill(Cl.sector, Cl.wheel);
        m_plots[Form("%s_BestQ_RemovedCluster", "HQFilter")]->Fill(Cl.bestTP().quality);
        
        if (Cl.bestTPQuality() > 1 ) {
          auto missingHQ = Cl.bestTPIndex();
          for (auto OriginalCluster : clusters) {
            if (missingHQ == OriginalCluster.bestTPIndex()) {
              std::cout << " è la best TP di un cluster " << std::endl; 
              std::cout << " Original cluster in " << OriginalCluster.wheel << " " << OriginalCluster.sector << " " << OriginalCluster.station << " xLoc:" << OriginalCluster.bestTP().xLoc << std::endl;
              std::cout << " HQFilter cluster in " << Cl.wheel << " " << Cl.sector << " " << Cl.station << " xLoc:" << Cl.bestTP().xLoc << std::endl;
              break;
            }
            else {
              bool wasghost = false;
              for (auto itGhost : OriginalCluster.itGhosts()) {
                if (missingHQ == itGhost.index ){
                  std::cout << " era un ghost in time BX = " << itGhost.BX << " quality = " << itGhost.quality << std::endl;
                  std::cout << " ora è la best TP in un cluster di dimensione " << Cl.tpClusterSize() << " che ha " << Cl.segClusterSize() << " segmenti e " << Cl.GetNDigi() << " digi " <<  endl;
                  wasghost = true;
                  break;
                } 
              }
              if (!wasghost ) {
                for (auto ootGhost : OriginalCluster.ootGhosts()){
                  if (missingHQ == ootGhost.index) std::cout << " era un ghost out of time BX = " << ootGhost.BX << std::endl; 
                }
              }
            }
          }
        }
      }

      m_plots[Form("%s_meanXLoc_RemovedClusters_st%d", "HQFilter" , ClusterCut_HQFilter[0].station)]->Fill(meanxLoc/ClusterCut_HQFilter.size());
    }

    for (auto &LQ : ClusterCut_LQFilter){
      for (auto &HQ :ClusterCut_HQFilter){
        if (! (LQ == HQ)) m_2Dplots[Form("LQnotinHQ_st%d", LQ.station)]->Fill(LQ.sector, LQ.wheel);
      }
    }


    for (auto &HQ : ClusterCut_HQFilter){
      for (auto &LQ :ClusterCut_LQFilter){
        if (!(HQ == LQ)) m_2Dplots[Form("HQnotinLQ_st%d", HQ.station)]->Fill(HQ.sector, HQ.wheel);
      }
    }


    // ########## PHI MATCHING TPs ANALYSIS #############
    for (TriggerPrimitive tp : tps) {
      if (tp.quality == 1) {
        m_plots["t0_LowQuality"]->Fill(tps.back().t0);
        m_plots["BX_LowQuality"]->Fill(tps.back().BX);
      }
    }
  }


  for (auto CLtype : tags){
  double ghostFraction = m_counters[Form("%s_nClustersGhosts", CLtype.c_str())] / m_counters[Form("%s_nClusters", CLtype.c_str())];
  cout << " \nFor " << CLtype << endl;
  cout << " Ratio LQ/selected with phi=  " << m_plots[Form("%s_LowQ_matched", CLtype.c_str())] ->GetEntries() << "/ "
       << m_plots["t0_LowQuality"]->GetEntries() << " = "
       << m_plots[Form("%s_LowQ_matched", CLtype.c_str())] ->GetEntries() / m_plots["t0_LowQuality"]->GetEntries() << endl;
  //cout << " Ratio LQ/matched with phi=  " << m_plots["LowQ_more1HQ_Phi"]->GetEntries() << "/ "
  //     << m_plots["t0_LowQuality"]->GetEntries() << " = "
  //     << m_plots["LowQ_more1HQ_Phi"]->GetEntries() / m_plots["t0_LowQuality"]->GetEntries() << endl;
  
  cout << " Fraction of clusters with ghost (" << m_counters[Form("%s_nClustersGhosts", CLtype.c_str())] << ") on total (" << m_counters[Form("%s_nClusters", CLtype.c_str())]
       << ") = " << ghostFraction << " made with "  << m_counters[Form("%s_nTPs", CLtype.c_str())]<< " TPs" << endl; 
  cout << " HQ out of time clusters: " << m_counters[Form("%s_ooTHQCount", CLtype.c_str())] << endl;
  }

  std::cout << "LQFilter removed " << removedLQFilter << " clusters, HQFilter removed " << removedHQFilter << " clusters " << std::endl;

  outputFile.Write();
  outputFile.Close();

}


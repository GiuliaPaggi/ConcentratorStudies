#define AnalyserBase_cxx
#include "include/Analyser.h"

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TLorentzVector.h>

#include <cassert>
#include <iostream>
#include <vector>

#include "include/Geometry.h"

// CB separate with comments what is about ghost cleaning,
//    clustering and analysis ...

// const std::array<double, 4> MB{402.2, 490.5, 597.5, 700.0}; CB not used
// const int LOW_QUAL_CUT{0};
// const double PSI_CUT{TMath::Pi() / 6.0};
// const double PHI_CUT_2{0.01};
const double PHI_CUT{0.02};
const int HIGH_QUAL_CUT{5};
const double T0_CUT{12.5};
const double X_CUT{10.0};
const double DIGI_CUT{15.0};
const int CORRECT_BX{20};

double BX_MIN{CORRECT_BX - 12};
double BX_MAX{CORRECT_BX + 12};

double T_MIN{BX_MIN * 25};
double T_MAX{BX_MAX * 25};

namespace MU {
const double MAX_ETA{0.8};
const double MIN_PT{5.0};
const double MAX_GEN_DR{0.25};
const double MASS{105.7};
const double X_CUT{5.0}; //cm
}  // namespace MU

const std::array<int, 3> QUAL_PLOT{0, 3, 6};

void Analyser::FillEfficiency(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster){
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  double BXTP = cluster.bestTP().BX;
  double qualTP = cluster.bestTPQuality();
  bool eff = (std::abs(BXTP - CORRECT_BX) < 1) && qualTP >= qual;

  m_effs[Form("%s_ClusterEfficiencyVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(eff, valueToFill);
}

void Analyser::FillEfficiency2D(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFillx, const double &valueToFilly, const Cluster &cluster){
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  double BXTP = cluster.bestTP().BX;
  double qualTP = cluster.bestTPQuality();
  bool eff = (std::abs(BXTP - CORRECT_BX) < 1) && qualTP >= qual;

  m_effs[Form("%s_ClusterEfficiencyVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(eff, valueToFillx, valueToFilly);
}

void Analyser::FillGhostRatio(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster){
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  double qualTP = cluster.bestTPQuality();
  bool eff = cluster.hasGhosts();
  bool ITeff = cluster.itSize();
  bool OOTeff = cluster.ootSize();
  if ( qualTP >= qual ) {
    m_effs[Form("%s_GhostFractionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(eff, valueToFill);
    m_effs[Form("%s_ITGhostFractionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(ITeff, valueToFill);
    m_effs[Form("%s_OOTGhostFractionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(OOTeff, valueToFill);
  }
}

void Analyser::FillGhostProfile(const std::string &typeStr, const std::string &varStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster){
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  const double ITghost = cluster.itSize();
  const double OOTghost = cluster.ootSize();

  double qualTP = cluster.bestTPQuality();
  if ( cluster.hasGhosts() && qualTP >= qual ) {
    m_plots[Form("%s_GhostDistributionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(valueToFill, ITghost + OOTghost);
    m_plots[Form("%s_ITGhostDistributionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(valueToFill, ITghost);
    m_plots[Form("%s_OOTGhostDistributionVS%s_st%d_minqual%d", type, var, st, qual)]->Fill(valueToFill, OOTghost);
  }
}

void Analyser::FillBackground(const std::string &typeStr, const int &st, const int &qual, const double &valueToFill, const Cluster &cluster){
  const auto type{typeStr.c_str()};
  
  const double ITtps = cluster.itSize();
  const double OOTtps = cluster.ootSize();
  if (cluster.bestTPQuality() >= qual){
    m_plots[Form("%s_ITBackground_st%d_minqual%d", type, st, qual)]->Fill(valueToFill, ITtps+1);
    m_plots[Form("%s_OOTBackground_st%d_minqual%d", type, st, qual)]->Fill(valueToFill, OOTtps);
  }
  if (ITtps+OOTtps > 0){
    const int bin = st*5+(cluster.wheel-2); //station*5+(wheel-2)
    m_plots[Form("%s_ITBackgroundDistribution", type)]->Fill(bin, ITtps);
    m_plots[Form("%s_OOTBackgroundDistribution", type)]->Fill(bin, OOTtps);
    if (cluster.bestTPQuality() < 0) {
      m_plots[Form("%s_NoBestTPBackgroundDistribution", type)]->Fill(bin, OOTtps);
    }    
  }

}

void Analyser::DefinePlot() {
  Geometry geom{};

  m_plots["DigiSL"] = new TH1D("DigiSL", "DigiSL", 6, -0.5, 5.5);
  m_plots["t0_LowQuality"] = new TH1D("t0_LowQuality", "t0_LowQuality; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["BX_LowQuality"] = new TH1D("BX_LowQuality", "BX_LowQuality; BX; Entries", 24, BX_MIN, BX_MAX);

  m_plots["ClusterPerEvent"] = new TH1I("ClusterPerEvent", "ClusterPerEvent; # clusters; Entries", 30, 0, 30);
  m_plots["GoodMuClusterPerEvent"] = new TH1I("GoodMuClusterPerEvent", "GoodMuClusterPerEvent; # clusters; Entries", 30, 0, 30);
  // Muon matching control plots
  m_plots["DeltaR"] = new TH1D("DeltaR", "DeltaR; #DeltaR(reco,gen); Entries", 100, -0.0, 10.0);
  m_plots["NGoodMu"] = new TH1D("NGoodMu", "NGoodMu; # of \"good\" muons; Entries", 11, -.5, 10.5);
  m_plots["NClusterMu"] = new TH1D("NClusterMu", "NClusterMu; # of cluster matchers per muon; Entries", 21, -.5, 20.5);
  m_plots["PtGoodMu"] = new TH1D("PtGoodMu", "PtGoodMu; pT (GeV); Entries", 200, .5, 100.5);
  m_plots["PhiGoodMu"] = new TH1D("PhiGoodMu", "PhiGoodMu; phi(rad); Entries", 50, -TMath::Pi(), TMath::Pi());
  m_plots["EtaGoodMu"] = new TH1D("EtaGoodMu", "EtaGoodMu; eta; Entries", 50, -1, 1);

  for (auto tag : tags) {
    const auto type{tag.c_str()};
    m_plots[Form("%s_N_DigiPerCluster", type)] = new TH1D(
        Form("%s_N_DigiPerCluster", type), Form("%s_N_DigiPerCluster; # digi in cluster; Entries", type), 50, 0, 50);
    m_plots[Form("%s_N_Ghost", type)] = new TH1I(Form("%s_N_Ghost", type), Form("%s_N_Ghost", type), 20, 0, 20);
    m_plots[Form("%s_Q_Best", type)] = new TH1I(Form("%s_Q_Best", type), Form("%s_Q_Best", type), 10, 0, 10);
    m_plots[Form("%s_Q_Ghost", type)] = new TH1I(Form("%s_Q_Ghost", type), Form("%s_Q_Ghost", type), 10, 0, 10);
    m_plots[Form("%s_OoTGhosts", type)] = new TH1D(Form("%s_OoTGhosts", type), Form("%s_OoTGhosts", type), 20, 0, 20);
    m_plots[Form("%s_ITGhosts", type)] = new TH1D(Form("%s_ITGhosts", type), Form("%s_ITGhosts", type), 20, 0, 20);
    m_plots[Form("%s_BX_ITGhosts", type)] =
        new TH1I(Form("%s_BX_ITGhosts", type), Form("%s_BX_ITGhosts", type), 24, BX_MIN, BX_MAX);
    m_plots[Form("%s_BX_OoTGhosts", type)] =
        new TH1I(Form("%s_BX_OoTGhosts", type), Form("%s_BX_OoTGhosts", type), 25, BX_MIN, BX_MAX);
    m_plots[Form("%s_Res_ITGhosts", type)] =
        new TH1D(Form("%s_Res_ITGhosts", type), Form("%s_Res_ITGhosts", type), 111, -5.5, 5.5);
    m_plots[Form("%s_Res_OoTGhosts", type)] =
        new TH1D(Form("%s_Res_OoTGhosts", type), Form("%s_Res_OoTGhosts", type), 111, -5.5, 5.5);
    m_plots[Form("%s_Res_SegMatched", type)] =
        new TH1D(Form("%s_Res_SegMatched", type),
                 Form("%s_Res_SegMatched; cluster.bestTP().xLoc - seg.xLoc; Entries", type), 101, -10, 10);
    m_plots[Form("%s_LowQ_matched", type)] =
        new TH1D(Form("%s_LowQ_matched", type), Form("%s_LowQ_matched; t0 (ns); Entries", type), 100, T_MIN, T_MAX);
    m_plots[Form("%s_BX_LowQ_matched", type)] =
        new TH1D(Form("%s_BX_LowQ_matched", type), Form("%s_BX_LowQ_matched,BX;Entries", type), 24, BX_MIN, BX_MAX);
    m_plots[Form("%s_t0_Selected", type)] =
        new TH1D(Form("%s_t0_Selected", type), Form("%s_t0_Selected; t0 (ns); Entries", type), 100, T_MIN, T_MAX);
    m_plots[Form("%s_ClusterSize", type)] =
        new TH1I(Form("%s_ClusterSize", type), Form("%s_ClusterSize; # TPs in cluster; Entries", type), 21, -.5, 20.5);

    m_plots[Form("%s_BackgroundDistribution", type)] = new TProfile(Form("%s_BackgroundDistribution", type), Form("%s_BackgroundDistribution", type), 20, 0, 21);
    m_plots[Form("%s_ITBackgroundDistribution", type)] = new TProfile(Form("%s_ITBackgroundDistribution", type), Form("%s_ITBackgroundDistribution", type), 20, 0, 21);
    m_plots[Form("%s_OOTBackgroundDistribution", type)] = new TProfile(Form("%s_OOTBackgroundDistribution", type), Form("%s_OOTBackgroundDistribution", type), 20, 0, 21);
    m_plots[Form("%s_NoBestTPBackgroundDistribution", type)] = new TProfile(Form("%s_NoBestTPBackgroundDistribution", type), Form("%s_NoBestTPBackgroundDistribution", type), 20, 0, 21);


    m_2Dplots[Form("%s_Q_OoTGhosts", type)] =
        new TH2D(Form("%s_Q_OoTGhosts", type), Form("%s_Q_OoTGhosts;High Quality;Out of time Ghost Quality", type), 10,
                 0, 10, 10, 0, 10);
    m_2Dplots[Form("%s_Q_ITGhosts", type)] =
        new TH2D(Form("%s_Q_ITGhosts", type), Form("%s_Q_ITGhosts;High Quality;In time Ghost Quality", type), 10, 0, 10,
                 10, 0, 10);

    m_counters[Form("%s_nClustersGhosts", type)] = 0;
    m_counters[Form("%s_ooTHQCount", type)] = 0;
    m_counters[Form("%s_nClusters", type)] = 0;
    m_counters[Form("%s_nTPs", type)] = 0;

    for (const auto st : geom.STATIONS) {
      m_plots[Form("%s_x_LowBestQ_st%d", type, st)] = new TH1D(
          Form("%s_x_LowBestQ_st%d", type, st), Form("%s_x_LowBestQ_st%d; xLoc; Entries", type, st), 100, -220, 220);
      m_plots[Form("%s_xLoc_ITGhost_st%d", type, st)] = new TH1D(
          Form("%s_xLoc_ITGhost_st%d", type, st), Form("%s_xLoc_ITGhost_st%d; xLoc; Entries", type, st), 100, -220, 220);
      m_plots[Form("%s_xLoc_OoTGhost_st%d", type, st)] = new TH1D(
          Form("%s_xLoc_OoTGhost_st%d", type, st), Form("%s_xLoc_OoTGhost_st%d; xLoc; Entries", type, st), 100, -220, 220);


      m_2Dplots[Form("%s_N_Cluster_st%d", type, st)] = new TH2I(
          Form("%s_N_Cluster_st%d", type, st), Form("%s_N_Cluster_st%d", type, st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("%s_N_Digi_st%d", type, st)] =
          new TH2I(Form("%s_N_Digi_st%d", type, st), Form("%s_N_Digi_st%d", type, st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st)] =
          new TH2D(Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st),
                   Form("AllTPs_LocalDirectionvsPosition_st%d", st), 200, -200, 200, 11, -1, 1);


      for (auto q : QUAL_PLOT){
        m_effs[Form("%s_ClusterEfficiencyVSpos_st%d_minqual%i", type, st, q)] =
            new TEfficiency(Form("%s_ClusterEfficiencyVSpos_st%d_minqual%d", type, st, q),
                            Form("%s_ClusterEfficiencyVSpos_st%d_minqual%d; sector; wheel", type, st, q), 14, -0.5, 13.5, 7, -3.5, 3.5);
        m_effs[Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d", type, st, q), Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d; pT (Gev); Efficiency", type, st, q), 50, MU::MIN_PT, 100);
        m_effs[Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d", type, st, q), Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d; phi (rad); Efficiency", type, st, q), 40, -TMath::Pi(), TMath::Pi());
        m_effs[Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d", type, st, q), Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d; eta; Efficiency", type, st, q), 50, -1, 1);
        m_effs[Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d", type, st, q), Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d; bestSeg_xLoc (cm); Efficiency", type, st, q), 100, -220, 220);
        m_effs[Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d", type, st, q), Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d; bestSeg_dirLoc; Efficieny", type, st, q), 40, -TMath::Pi(), TMath::Pi()); 
        
        m_effs[Form("%s_GhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_GhostFractionVSpT_st%d_minqual%d", type, st, q), Form("%s_GhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with ghosts", type, st, q), 50, MU::MIN_PT, 100);
        m_effs[Form("%s_ITGhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_ITGhostFractionVSpT_st%d_minqual%d", type, st, q), Form("%s_ITGhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with in-time ghosts", type, st, q), 50, MU::MIN_PT, 100);
        m_effs[Form("%s_OOTGhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(Form("%s_OOTGhostFractionVSpT_st%d_minqual%d", type, st, q), Form("%s_OOTGhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with out-of-time ghosts", type, st, q), 50, MU::MIN_PT, 100);

        m_plots[Form("%s_GhostDistributionVSpT_st%d_minqual%d", type, st, q)] = new TProfile(Form("%s_GhostDistributionVSpT_st%d_minqual%d", type, st, q), Form("%s_GhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # ghosts", type, st, q), 50, MU::MIN_PT, 100);
        m_plots[Form("%s_ITGhostDistributionVSpT_st%d_minqual%d", type, st, q)] = new TProfile(Form("%s_ITGhostDistributionVSpT_st%d_minqual%d", type, st, q), Form("%s_ITGhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # in-time ghosts", type, st, q), 50, MU::MIN_PT, 100);
        m_plots[Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d", type, st, q)] = new TProfile(Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d", type, st, q), Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # out-of-time ghosts", type, st, q), 50, MU::MIN_PT, 100);

        m_plots[Form("%s_ITBackground_st%d_minqual%d", type, st, q)] = new TH1F(Form("%s_ITBackground_st%d_minqual%d", type, st, q), Form("%s_ITBackground_st%d_minqual%d; xLoc (cm); Entries", type, st, q), 100, -220, 220);      
        m_plots[Form("%s_OOTBackground_st%d_minqual%d", type, st, q)] = new TH1F(Form("%s_OOTBackground_st%d_minqual%d", type, st, q), Form("%s_OOTBackground_st%d_minqual%d; xLoc (cm); Entries", type, st, q), 100, -220, 220);      

      }

      m_effs[Form("%s_Eff_SegMatch_st%d", type, st)] =
          new TEfficiency(Form("%s_Eff_SegMatch_st%d", type, st), Form("%s_Eff_SegMatch_st%d; sector; wheel", type, st),
                          14, -0.5, 13.5, 7, -3.5, 3.5);
      m_effs[Form("%s_Eff_DigiMatch_st%d", type, st)] =
          new TEfficiency(Form("%s_Eff_DigiMatch_st%d", type, st),
                          Form("%s_Eff_DigiMatch_st%d; sector; wheel", type, st), 14, -0.5, 13.5, 7, -3.5, 3.5);

      for (const auto wh : geom.WHEELS) {
        m_plots[Form("%s_Digi_residual_st%d_wh%d", type, st, wh)] = new TH1D(
            Form("%s_Digi_residual_st%d_wh%d", type, st, wh),
            Form("%s_Digi_residual_st%d_wh%d; cluster.bestSeg.xLoc - digi.xLoc; entries", type, st, wh), 51, -10, 10);
      }
    }
  }

  for (const auto st : geom.STATIONS) {
    m_plots[Form("Res_MuMatched_st%d", st)] = new TH1D(
        Form("Res_MuMatched_st%d", st), Form("Res_MuMatched; cluster.bestTP().xLoc - mu.xLoc; Entries"), 101, -10, 10);

    m_2Dplots[Form("N_MuMatch_st%d", st)] =
        new TH2I(Form("N_MuMatch_st%d", st), Form("N_MuMatch_st%d; sector ; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    
    m_2Dplots[Form("PosGoodMu_st%d", st)] = new TH2D(Form("PosGoodMu_st%d", st), Form("PosGoodMu_st%d; eta; phi (rad)", st), 50, -1, 1, 50, -TMath::Pi(), TMath::Pi() ); 
    m_2Dprofiles[Form("NMatchGoodMu_st%d", st)] = new TProfile2D(Form("NMatchGoodMu_st%d", st), Form("NMatchGoodMu_st%d; eta; phi (rad)", st), 50, -1, 1, 50, -TMath::Pi(), TMath::Pi(), 0, 25);

    m_effs[Form("Eff_MuMatch_st%d", st)] = new TEfficiency(
        Form("Eff_MuMatch_st%d", st), Form("Eff_MuMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);

    }
}

void Analyser::ClusterAnalisis(const std::vector<Cluster> &clusters, const std::string &typeStr,
                               const std::vector<Segment> &segments) {
  for (auto const &cluster : clusters) {
    auto wh{cluster.wheel};
    auto sec{cluster.sector};
    assert(sec < 13);  // if sector == 13 o 14 the program crashes
    auto st{cluster.station};
    const auto type{typeStr.c_str()};

    TriggerPrimitive bestTP = cluster.bestTP();
    Segment bestSeg = cluster.bestSeg();

    m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st)]->Fill(bestTP.xLoc,
                                                                            bestTP.psi);
    m_plots[Form("%s_ClusterSize", type)]->Fill(cluster.tpClusterSize());
    // ########## Study TP ghost distribution #############
    if (cluster.foundTP) ++m_counters[Form("%s_nClusters", type)];
    m_counters[Form("%s_ooTHQCount", type)] +=
        cluster.ootCountIf([=](TriggerPrimitive const &tp) { return tp.quality > HIGH_QUAL_CUT; });

    int bestQ{cluster.bestTPQuality()};
    int ootSize{cluster.ootSize()};
    int itSize{cluster.itSize()};
    
    m_plots[Form("%s_Q_Best", type)]->Fill(bestQ);
    m_plots[Form("%s_ITGhosts", type)]->Fill(itSize);
    m_plots[Form("%s_OoTGhosts", type)]->Fill(ootSize);
    m_plots[Form("%s_N_Ghost", type)]->Fill(ootSize + itSize);
    m_counters[Form("%s_nTPs", type)] += cluster.tpClusterSize();

    if (bestQ == 1) m_plots[Form("%s_x_LowBestQ_st%d", type, st)]->Fill(bestTP.xLoc);

    m_2Dplots[Form("%s_N_Cluster_st%d", type, st)]->Fill(sec, wh);

    if (cluster.hasGhosts()) {
      ++m_counters[Form("%s_nClustersGhosts", type)];

      for (const auto &ghost : cluster.ootGhosts()) {
        m_plots[Form("%s_BX_OoTGhosts", type)]->Fill(ghost.BX);
        m_plots[Form("%s_Res_OoTGhosts", type)]->Fill(ghost.xLoc - bestTP.xLoc);
        m_2Dplots[Form("%s_Q_OoTGhosts", type)]->Fill(bestQ, ghost.quality);
        m_plots[Form("%s_Q_Ghost", type)]->Fill(ghost.quality);
        m_plots[Form("%s_xLoc_OoTGhost_st%d", type, st)]->Fill(ghost.xLoc);

        if (ghost.quality == 1) {
          m_plots[Form("%s_t0_Selected", type)]->Fill(ghost.t0);
          m_plots[Form("%s_LowQ_matched", type)]->Fill(ghost.t0);
          m_plots[Form("%s_BX_LowQ_matched", type)]->Fill(ghost.BX);
        }
      }

      for (const auto &ghost : cluster.itGhosts()) {
        m_plots[Form("%s_BX_ITGhosts", type)]->Fill(ghost.BX);
        m_plots[Form("%s_Res_ITGhosts", type)]->Fill(ghost.xLoc - bestTP.xLoc);
        m_2Dplots[Form("%s_Q_ITGhosts", type)]->Fill(bestQ, ghost.quality);
        m_plots[Form("%s_Q_Ghost", type)]->Fill(ghost.quality);
        m_plots[Form("%s_xLoc_ITGhost_st%d", type, st)]->Fill(ghost.xLoc);

        if (ghost.quality == 1) {
          m_plots[Form("%s_t0_Selected", type)]->Fill(ghost.t0);
          m_plots[Form("%s_LowQ_matched", type)]->Fill(ghost.t0);
          m_plots[Form("%s_BX_LowQ_matched", type)]->Fill(ghost.BX);
        }
      }
    
    }

    if (cluster.muMatched && bestSeg.nPhiHits >= 4) {      
      int muon = cluster.MuIndex();
      double bestSegIndex = cluster.bestSegIndex();
      double dirLoc = std::atan( seg_dirLoc_x->at(bestSegIndex) / seg_dirLoc_z->at(bestSegIndex));
        for (auto qual : QUAL_PLOT){
            FillEfficiency2D(type, "pos", st, qual, sec, wh, cluster);
            FillEfficiency(type, "pT", st, qual, mu_pt->at(muon), cluster);
            FillEfficiency(type, "phi", st, qual, mu_phi->at(muon), cluster);
            FillEfficiency(type, "eta", st, qual, mu_eta->at(muon), cluster);
            FillEfficiency(type, "segxLoc", st, qual, bestSeg.xLoc, cluster);
            FillEfficiency(type, "segDirLoc", st, qual, dirLoc, cluster);
            FillGhostRatio(type, "pT", st, qual, mu_pt->at(muon), cluster);
            FillGhostProfile(type, "pT", st, qual, mu_pt->at(muon), cluster);
        }
    }
    
    // ########## Study background #############
    
    if (!cluster.muMatched){
      double xLoc = cluster.bestTP().xLoc;
      for (auto qual : QUAL_PLOT){
        FillBackground(type, st, qual, xLoc, cluster);
      }
    }

    // ########## Study segment matching #############
    m_effs[Form("%s_Eff_SegMatch_st%d", type, st)]->Fill(cluster.foundSeg, sec, wh);
    for (const auto s : segments)
      if (cluster.foundTP && s.wheel == wh && s.sector == sec && s.station == st) {
        m_plots[Form("%s_Res_SegMatched", type)]->Fill(bestTP.xLoc - s.xLoc);
      }

    // ########## Study digi clusters #############
    m_effs[Form("%s_Eff_DigiMatch_st%d", type, st)]->Fill(cluster.foundDigi, sec, wh);

    if (cluster.foundDigi) {
      for (auto &digi : cluster.matchedDigi()) {
        m_plots[Form("%s_Digi_residual_st%d_wh%d", type, st, wh)]->Fill(bestSeg.xLoc - digi.xLoc);
      }
    }

    m_plots[Form("%s_N_DigiPerCluster", type)]->Fill(cluster.nDigi());
    m_plots["DigiSL"]->Fill(cluster.digiSL());

    if (cluster.foundDigi) m_2Dplots[Form("%s_N_Digi_st%d", type, st)]->Fill(sec, wh);
  }
}

std::vector<Cluster> MissingClusters(std::vector<Cluster> FirstCLVector, std::vector<Cluster> SecondCLvector,
                                     std::string typeStr) {
  // FirstCL is the smaller set , SecondCL is the bigger one
  std::vector<Cluster> missingClusters;
  const auto type{typeStr.c_str()};

  for (const auto &FirstCluster : FirstCLVector) {
    if (!FirstCluster.foundTP) continue;
    auto wh{FirstCluster.wheel};
    auto sec{FirstCluster.sector};
    auto st{FirstCluster.station};

    int n_cluster_uguali =
        std::count_if(SecondCLvector.begin(), SecondCLvector.end(), [&](const auto &cl) { return cl == FirstCluster; });
    if (n_cluster_uguali == 0) {
      missingClusters.push_back(FirstCluster);
      std::cout << "Missing cluster for " << type << "\n" << FirstCluster << std::endl;

      std::cout << "In the first cluster there are also " << std::endl;
      for (const auto &otherFirstclusters : FirstCLVector) {
        if (otherFirstclusters.wheel != wh || otherFirstclusters.station != st || otherFirstclusters.sector != sec)
          continue;
        if (otherFirstclusters == FirstCluster) continue;
        std::cout << otherFirstclusters << std::endl;
      }

      std::cout << "In the original cluster there were:" << std::endl;
      for (const auto &SecondCluster : SecondCLvector) {
        if (SecondCluster.wheel != wh || SecondCluster.station != st || SecondCluster.sector != sec) continue;
        std::cout << SecondCluster << std::endl;
      }
      std::cout << " ################################################################################# " << std::endl;
    }
  }
  return missingClusters;
}

void Analyser::Loop() {
  Geometry geom{};

  TFile outputFile("results/outputFile_DoubleMuon_FlatPt-1To100_noPU_noRPC.root", "RECREATE");
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

  for (Long64_t jentry = 0; jentry < n_entries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    // ########## LOOP ON EVENTS #############
    if (jentry % 100 == 0) std::cout << "Processing event: " << jentry << '\r' << std::flush;
    if (gen_pdgId->size() < 1) continue;

    // ########## CREATE TPs std::vector #############
    std::vector<TriggerPrimitive> tps;

    for (std::size_t j = 0; j < ph2TpgPhiEmuAm_nTrigs; ++j) {
      tps.emplace_back(TriggerPrimitive{j, ph2TpgPhiEmuAm_wheel->at(j), ph2TpgPhiEmuAm_sector->at(j),
                                        ph2TpgPhiEmuAm_station->at(j), ph2TpgPhiEmuAm_quality->at(j),
                                        ph2TpgPhiEmuAm_phi->at(j), ph2TpgPhiEmuAm_phiB->at(j), ph2TpgPhiEmuAm_BX->at(j),
                                        ph2TpgPhiEmuAm_t0->at(j), ph2TpgPhiEmuAm_posLoc_x->at(j)});
    }

    // ########## CREATE Digis std::vector #############
    std::vector<Digi> digis;
    for (std::size_t i = 0; i < digi_nDigis; ++i) {
      digis.emplace_back(Digi(geom, i, digi_wheel->at(i), digi_sector->at(i), digi_station->at(i),
                              digi_superLayer->at(i), digi_layer->at(i), digi_wire->at(i), digi_time->at(i)));
    }

    // ########## BUILD segments std::vector #############
    std::vector<Segment> segments;
    for (std::size_t j = 0; j < seg_nSegments; ++j) {
      segments.emplace_back(Segment(j, seg_station->at(j), seg_wheel->at(j), seg_sector->at(j), seg_phi_nHits->at(j),
                                    seg_posLoc_x->at(j)));
    }

    std::vector<Cluster> matchFromLQ_clusters;
    std::vector<Cluster> matchFromHQ_clusters;

    // ########## BUILD clusters std::vector #############
    auto clusters = buildClusters(geom, tps, segments, digis, X_CUT, DIGI_CUT);
    int nClusters = clusters.size();
    m_plots["ClusterPerEvent"]->Fill(nClusters);
    // ########## BUILD cluster with phi matchin information #############
    // tag TPs that match using the extrapolation on a straight line

    std::vector<TriggerPrimitive> MatchFromLQ_tps = tps;
    std::vector<TriggerPrimitive> MatchFromHQ_tps = tps;

    // extrapolate from LQ TP
    for (TriggerPrimitive &tp : MatchFromLQ_tps) {
      if (!tp.hasMatched) {
        for (TriggerPrimitive &other_tp : MatchFromLQ_tps) {
          if (tp.index != other_tp.index) {
            tp.MatchFromLQ(other_tp, PHI_CUT, T0_CUT);
          }
        }
      }
    }

    MatchFromLQ_tps.erase(std::remove_if(MatchFromLQ_tps.begin(), MatchFromLQ_tps.end(),
                                         [](auto &tp) { return !tp.hasMatched && tp.quality == 1; }),
                          MatchFromLQ_tps.end());
    // repeat clustering after filtering
    matchFromLQ_clusters = buildClusters(geom, MatchFromLQ_tps, segments, digis, X_CUT, DIGI_CUT);

    // extrapolate from HQ TP
    for (TriggerPrimitive &tp : MatchFromHQ_tps) {
      if (!tp.hasMatched) {
        for (TriggerPrimitive &other_tp : MatchFromHQ_tps) {
          if (tp.index != other_tp.index) {
            tp.MatchFromHQ(other_tp, PHI_CUT, T0_CUT);
          }
        }
      }
    }

    MatchFromHQ_tps.erase(std::remove_if(MatchFromHQ_tps.begin(), MatchFromHQ_tps.end(),
                                         [](auto &tp) { return !tp.hasMatched && tp.quality == 1; }),
                          MatchFromHQ_tps.end());
    // repeat clustering after filtering
    matchFromHQ_clusters = buildClusters(geom, MatchFromHQ_tps, segments, digis, X_CUT, DIGI_CUT);

    std::vector<TriggerPrimitive> TPRemovedfromLQ = tps;
    for (auto tp : MatchFromLQ_tps)
      TPRemovedfromLQ.erase(std::remove(TPRemovedfromLQ.begin(), TPRemovedfromLQ.end(), tp), TPRemovedfromLQ.end());
    std::vector<TriggerPrimitive> TPRemovedfromHQ = tps;
    for (auto tp : MatchFromHQ_tps)
      TPRemovedfromHQ.erase(std::remove(TPRemovedfromHQ.begin(), TPRemovedfromHQ.end(), tp), TPRemovedfromHQ.end());

    // ########## Muon ID and GEN matching #############
    std::vector<UInt_t> goodMuons;
    for (UInt_t iMu = 0; iMu < mu_nMuons; ++iMu) {
      if (!mu_isMedium->at(iMu)) continue;
      TLorentzVector mu;
      double muEta = mu_eta->at(iMu);
      double muPhi = mu_phi->at(iMu);
      double muPT = mu_pt->at(iMu);
      mu.SetPtEtaPhiM(muPT, muEta, muPhi, MU::MASS);
      // check if the muon can be matched to one of the particles generated in this event
      for (UInt_t igen = 0; igen < gen_nGenParts; ++igen) {
        if (std::abs(gen_pdgId->at(igen)) != 13 || std::abs(gen_eta->at(igen)) > MU::MAX_ETA ||
            gen_pt->at(igen) < MU::MIN_PT)
          continue;
        double genEta = gen_eta->at(igen);
        double genPhi = gen_phi->at(igen);
        double genPT = gen_pt->at(igen);
        TLorentzVector genMu;
        genMu.SetPtEtaPhiM(genPT, genEta, genPhi, MU::MASS);
        double dR = genMu.DeltaR(mu);
        m_plots["DeltaR"]->Fill(dR);
        if (std::abs(dR) < MU::MAX_GEN_DR) {
          goodMuons.push_back(iMu);
          break;
        }
      }
    }
    
    int goodMuSize = goodMuons.size();
    m_plots["NGoodMu"]->Fill(goodMuSize);
    if  (goodMuSize > 0) m_plots["GoodMuClusterPerEvent"]->Fill(nClusters); 

    // ########## ATTEMPT cluster - muon extrapolation matching #############
    for (auto iMu : goodMuons) {
      double muPhi = mu_phi->at(iMu);
      double muEta = mu_eta->at(iMu);
      double muPT = mu_pt->at(iMu);

      m_plots["PtGoodMu"]->Fill(muPT);
      m_plots["PhiGoodMu"]->Fill(muPhi);
      m_plots["EtaGoodMu"]->Fill(muEta);

      Int_t nMatchedClusters{};
      std::array<int, 4> stationMatch{0};

      for (UInt_t i = 0; i < mu_nMatches->at(iMu); ++i) {
        int muWheel = getXY<float>(mu_matches_wheel, iMu, i);
        int muStation = getXY<float>(mu_matches_station, iMu, i);
        int muSector = getXY<float>(mu_matches_sector, iMu, i);
        if (muSector == 13) muSector = 4;
        if (muSector == 14) muSector = 10;
        double muX = getXY<float>(mu_matches_x, iMu, i);
        double edgeX = getXY<float>(mu_matches_edgeX, iMu, i);
        double edgeY = getXY<float>(mu_matches_edgeY, iMu, i);

        for (auto &cluster : clusters) {
          auto wh{cluster.wheel};
          auto sec{cluster.sector};
          if (sec == 13) sec = 4;
          if (sec == 14) sec = 10;
          auto st{cluster.station};
          if (cluster.matchMu(muWheel, muStation, muSector, edgeX, edgeY, muX, i, iMu, MU::X_CUT)) {
            ++nMatchedClusters;
            ++stationMatch[muStation-1];
            m_2Dplots[Form("PosGoodMu_st%d", muStation)]->Fill(muEta, muPhi);
            if (!cluster.foundTP) std::cout << " Matchato un cluster senza best TP " << std::endl;
          };

          m_effs[Form("Eff_MuMatch_st%d", st)]->Fill(cluster.muMatched, sec, wh);
          if (cluster.muMatched) {
            m_plots[Form("Res_MuMatched_st%d", st)]->Fill(
                cluster.bestTP().xLoc -
                getXY<float>(mu_matches_x, cluster.muMatchedIndex[0], cluster.muMatchedIndex[1]));
            m_2Dplots[Form("N_MuMatch_st%d", st)]->Fill(sec, wh);
          }
        }

        for (auto &PMcluster : matchFromLQ_clusters) {
          PMcluster.matchMu(muWheel, muStation, muSector, edgeX, edgeY, muX, i, iMu, MU::X_CUT);
        }

        for (auto &PMcluster : matchFromHQ_clusters) {
          PMcluster.matchMu(muWheel, muStation, muSector, edgeX, edgeY, muX, i, iMu, MU::X_CUT);
        }
      }
      
      m_plots["NClusterMu"]->Fill(nMatchedClusters);
      
      for ( int st = 1; st < 5; ++st ) {
        m_2Dprofiles[Form("NMatchGoodMu_st%d", st)]->Fill(muEta, muPhi, stationMatch[st-1], 1);
      }
    }

    // ########## RUN SOME ANALYSIS #############
    // ########## CLUSTER ANALYSIS #############

    ClusterAnalisis(clusters, tags[0], segments);
    ClusterAnalisis(matchFromLQ_clusters, tags[1], segments);
    ClusterAnalisis(matchFromHQ_clusters, tags[2], segments);

    std::vector<Cluster> ClusterCut_LQFilter = MissingClusters(matchFromLQ_clusters, clusters, "LQFilter");
    removedLQFilter += ClusterCut_LQFilter.size();

    std::vector<Cluster> ClusterCut_HQFilter = MissingClusters(matchFromHQ_clusters, clusters, "HQFilter");
    removedHQFilter += ClusterCut_HQFilter.size();

    // ########## PHI MATCHING TPs ANALYSIS #############
    for (TriggerPrimitive tp : tps) {
      if (tp.quality == 1) {
        m_plots["t0_LowQuality"]->Fill(tps.back().t0);
        m_plots["BX_LowQuality"]->Fill(tps.back().BX);
      }
    }
  }

  for (auto tag : tags) {
    const auto type{tag.c_str()};
    double ghostFraction = m_counters[Form("%s_nClustersGhosts", type)] / m_counters[Form("%s_nClusters", type)];
    std::cout << " \nFor " << tag << std::endl;
    std::cout << " Ratio LQ/selected with phi=  " << m_plots[Form("%s_LowQ_matched", type)]->GetEntries() << "/ "
              << m_plots["t0_LowQuality"]->GetEntries() << " = "
              << m_plots[Form("%s_LowQ_matched", type)]->GetEntries() / m_plots["t0_LowQuality"]->GetEntries()
              << std::endl;

    std::cout << " Fraction of clusters with ghost (" << m_counters[Form("%s_nClustersGhosts", type)] << ") on total ("
              << m_counters[Form("%s_nClusters", type)] << ") = " << ghostFraction << " made with "
              << m_counters[Form("%s_nTPs", type)] << " TPs" << std::endl;
    std::cout << " HQ out of time clusters: " << m_counters[Form("%s_ooTHQCount", type)] << std::endl;
  }

  std::cout << "LQFilter found " << removedLQFilter << " less clusters, HQFilter found " << removedHQFilter
            << " less clusters " << std::endl;

  outputFile.Write();
  outputFile.Close();
}

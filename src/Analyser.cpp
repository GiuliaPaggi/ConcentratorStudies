#define AnalyserBase_cxx
#include "include/Analyser.h"

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>

#include <cassert>
#include <iostream>
#include <filesystem> 

#include "include/Geometry.h"

// CB separate with comments what is about ghost cleaning,
//    clustering and analysis ...

const double PHI_CUT{0.02};
const int HIGH_QUAL_CUT{5};
const double T0_CUT{12.5};
const double X_CUT{20.0};
const double DIGI_CUT{25.0};
const int CORRECT_BX{20};

double BX_MIN{CORRECT_BX - 12};
double BX_MAX{CORRECT_BX + 12};

double T_MIN{BX_MIN * 25};
double T_MAX{BX_MAX * 25};

namespace MU {
const double MAX_ETA{1.4};
const double MIN_PT{1.5};
const double MAX_GEN_DR{0.25};
const double MASS{105.7};
const double X_CUT{20.0};  // cm
}  // namespace MU

namespace EFF {
const double MAX_ETA{0.8};
const double MIN_PT{5.0};
const double XEDGE_CUT{5.0};
const double YEDGE_CUT{5.0};
}  // namespace EFF

const std::array<int, 3> QUAL_PLOT{0, 3, 6};

void Analyser::FillEfficiency(const std::string &typeStr, const std::string &varStr, int st, int qualCut,
                              double valueToFill, const Cluster &cluster) {
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  const int BX = cluster.bestTP().BX;
  const int qual = cluster.bestTPQuality();
  const bool eff = (std::abs(BX - CORRECT_BX) < 1) && qual >= qualCut;

  m_effs[Form("%s_ClusterEfficiencyVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(eff, valueToFill);
}

void Analyser::FillEfficiency2D(const std::string &typeStr, const std::string &varStr, int st, int qualCut,
                                double valueToFillx, double valueToFilly, const Cluster &cluster) {
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  const int BX = cluster.bestTP().BX;
  const int qual = cluster.bestTPQuality();
  const bool eff = (std::abs(BX - CORRECT_BX) < 1) && qual >= qualCut;

  m_effs[Form("%s_ClusterEfficiencyVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(eff, valueToFillx,
                                                                                        valueToFilly);
}

void Analyser::FillGhostRatio(const std::string &typeStr, const std::string &varStr, int st, int qualCut,
                              double valueToFill, const Cluster &cluster) {
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  const int qual = cluster.bestTPQuality();
  if (qual >= qualCut) {
    m_effs[Form("%s_GhostFractionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(cluster.hasGhosts(), valueToFill);
    m_effs[Form("%s_ITGhostFractionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(cluster.itSize(), valueToFill);
    m_effs[Form("%s_OOTGhostFractionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(cluster.ootSize(),
                                                                                         valueToFill);
  }
}

void Analyser::FillGhostProfile(const std::string &typeStr, const std::string &varStr, int st, int qualCut,
                                double valueToFill, const Cluster &cluster) {
  const auto type{typeStr.c_str()};
  const auto var{varStr.c_str()};

  const double ITghost = cluster.itSize();
  const double OOTghost = cluster.ootSize();

  const int qual = cluster.bestTPQuality();
  if (cluster.hasGhosts() && qual >= qualCut) {
    m_plots[Form("%s_GhostDistributionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(valueToFill,
                                                                                           ITghost + OOTghost);
    m_plots[Form("%s_ITGhostDistributionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(valueToFill, ITghost);
    m_plots[Form("%s_OOTGhostDistributionVS%s_st%d_minqual%d", type, var, st, qualCut)]->Fill(valueToFill, OOTghost);
  }
}

void Analyser::FillBackground(const std::string &typeStr, int st, int qualCut, double valueToFill,
                              const Cluster &cluster) {
  const auto type{typeStr.c_str()};

  const int ITtps{cluster.itSize() + (cluster.foundTP ? 1 : 0)};
  const int OOTtps{cluster.ootSize()};

  if (cluster.bestTPQuality() >= qualCut) {
    m_plots[Form("%s_ITBackground_st%d_minqual%d", type, st, qualCut)]->Fill(valueToFill, ITtps);
    m_plots[Form("%s_OOTBackground_st%d_minqual%d", type, st, qualCut)]->Fill(valueToFill, OOTtps);
  }
}

void Analyser::DefinePlot() {

  // ############# independent from filter #############
  m_counters["%noGoodMuTPS"] = 0;

  m_plots["t0_LowQuality"] = new TH1D("t0_LowQuality", "t0_LowQuality; t0 (ns); Entries", 100, T_MIN, T_MAX);
  m_plots["BX_LowQuality"] = new TH1D("BX_LowQuality", "BX_LowQuality; BX; Entries", 24, BX_MIN, BX_MAX);

  m_plots["ClusterPerEvent"] = new TH1I("ClusterPerEvent", "ClusterPerEvent; # clusters; Entries", 30, 0, 30);
  m_plots["GoodMuClusterPerEvent"] =
      new TH1I("GoodMuClusterPerEvent", "GoodMuClusterPerEvent; # clusters; Entries", 30, 0, 30);

  // Muon matching control plots
  m_plots["DeltaR"] = new TH1D("DeltaR", "DeltaR; #DeltaR(reco,gen); Entries", 100, -0.0, 10.0);
  m_plots["NGoodMu"] = new TH1D("NGoodMu", "NGoodMu; # of \"good\" muons; Entries", 11, -.5, 10.5);
  m_plots["NClusterMu"] = new TH1D("NClusterMu", "NClusterMu; # of cluster matchers per muon; Entries", 21, -.5, 20.5);
  m_plots["PtGoodMu"] = new TH1D("PtGoodMu", "PtGoodMu; pT (GeV); Entries", 200, .5, 100.5);
  m_plots["PhiGoodMu"] = new TH1D("PhiGoodMu", "PhiGoodMu; phi(rad); Entries", 50, -TMath::Pi(), TMath::Pi());
  m_plots["EtaGoodMu"] = new TH1D("EtaGoodMu", "EtaGoodMu; eta; Entries", 50, -1, 1);

  //prefiring
  m_plots["FirstTpBX"] = new TH1I("FirstTpBX", "FirstTpBX", 10, CORRECT_BX-5, CORRECT_BX+5);
  m_plots["FirstTpt0"] = new TH1I("FirstTpt0", "FirstTpt0", 100, (CORRECT_BX-5)*25-.5, (CORRECT_BX+5)*25-.5);

  // background
  m_2Dplots["Pos_unmatchedGenMu"] = new TH2D("Pos_unmatchedGenMu", "Pos_unmatchedGenMu; eta; phi (rad)", 50, -1, 1, 50, -TMath::Pi(),
                 TMath::Pi());
  m_plots["Pt_unmatchedGenMu"] = new TH1D("Pt_unmatchedGenMu", "Pt_unmatchedGenMu; pT (GeV); Entries", 200, .5, 100.5);
  

  // ############# Plots per station #############
  for (const auto st : GEOM.STATIONS) {
      m_2Dplots[Form("N_MuMatch_st%d", st)] =
          new TH2I(Form("N_MuMatch_st%d", st), Form("N_MuMatch_st%d; sector ; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("PosGoodMu_st%d", st)] =
          new TH2D(Form("PosGoodMu_st%d", st), Form("PosGoodMu_st%d; eta; phi (rad)", st), 50, -1, 1, 50, -TMath::Pi(),
                 TMath::Pi());
      m_2Dprofiles[Form("NMatchGoodMu_st%d", st)] =
        new TProfile2D(Form("NMatchGoodMu_st%d", st), Form("NMatchGoodMu_st%d; eta; phi (rad)", st), 50, -1, 1, 50,
                       -TMath::Pi(), TMath::Pi(), 0, 25);
      m_effs[Form("Eff_MuMatch_st%d", st)] = new TEfficiency(
          Form("Eff_MuMatch_st%d", st), Form("Eff_MuMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);

      //residuals
      m_plots[Form("Res_MuMatched_st%d", st)] =
        new TH1D(Form("Res_MuMatched_st%d", st), Form("Res_MuMatched_st%d; cluster.bestTP().xLoc - mu.xLoc; Entries", st), 200,
                 -220, 220);

      m_plots[Form("Res_Mu_st%d", st)] =
        new TH1D(Form("Res_Mu_st%d", st), Form("Res_Mu; cluster.bestTP().xLoc - mu.xLoc; Entries"), 200, -220, 220);
  }

  // ############# Plots to compare filters #############
  for (auto tag : tags) {
    const auto type{tag.c_str()};
    
    // Cluster properties
    m_plots[Form("%s_N_DigiPerCluster", type)] = new TH1D(
        Form("%s_N_DigiPerCluster", type), Form("%s_N_DigiPerCluster; # digi in cluster; Entries", type), 50, 0, 50);
    m_plots[Form("%s_N_Ghosts", type)] = new TH1I(Form("%s_N_Ghosts", type), Form("%s_N_Ghosts", type), 20, 0, 20);
    m_plots[Form("%s_Q_Best", type)] = new TH1I(Form("%s_Q_Best", type), Form("%s_Q_Best", type), 10, 0, 10);
    m_plots[Form("%s_Q_Ghost", type)] = new TH1I(Form("%s_Q_Ghost", type), Form("%s_Q_Ghost", type), 10, 0, 10);
    m_plots[Form("%s_OOTGhosts", type)] = new TH1I(Form("%s_OOTGhosts", type), Form("%s_OOTGhosts", type), 20, 0, 20);
    m_plots[Form("%s_ITGhosts", type)] = new TH1I(Form("%s_ITGhosts", type), Form("%s_ITGhosts", type), 20, 0, 20);
    m_plots[Form("%s_ClusterSize", type)] =
        new TH1I(Form("%s_ClusterSize", type), Form("%s_ClusterSize; # TPs in cluster; Entries", type), 21, -.5, 20.5);
    
    m_counters[Form("%s_nClustersGhosts", type)] = 0;
    m_counters[Form("%s_ooTHQCount", type)] = 0;
    m_counters[Form("%s_nClusters", type)] = 0;
    m_counters[Form("%s_nTPs", type)] = 0;

    // residuals 
    m_plots[Form("%s_Res_ITGhosts", type)] =
        new TH1D(Form("%s_Res_ITGhosts", type), Form("%s_Res_ITGhosts", type), 110, -10.5, 10.5);
    m_plots[Form("%s_Res_OoTGhosts", type)] =
        new TH1D(Form("%s_Res_OoTGhosts", type), Form("%s_Res_OoTGhosts", type), 110, -10.5, 10.5);
    m_plots[Form("%s_Res_SegMatched", type)] =
        new TH1D(Form("%s_Res_SegMatched", type),
                 Form("%s_Res_SegMatched; cluster.bestTP().xLoc - seg.xLoc; Entries", type), 110, -10.5, 10.5);

    //background

    m_plots[Form("%s_Q_Bkg", type)] = new TH1I(Form("%s_Q_Bkg", type), Form("%s_Q_Bkg", type), 10, 0, 10);
    m_plots[Form("%s_N_bkg", type)] = new TH1I(Form("%s_N_bkg", type), Form("%s_N_bkg", type), 20, 0, 20);
    m_2Dplots[Form("%s_Q_OoTGhosts", type)] =
        new TH2D(Form("%s_Q_OoTGhosts", type), Form("%s_Q_OoTGhosts;High Quality;Out of time Ghost Quality", type), 10,
                 0, 10, 10, 0, 10);
    m_2Dplots[Form("%s_Q_ITGhosts", type)] =
        new TH2D(Form("%s_Q_ITGhosts", type), Form("%s_Q_ITGhosts;High Quality;In time Ghost Quality", type), 10, 0, 10,
                 10, 0, 10);

    m_counters[Form("%s_nBackgroundClusters", type)] = 0;
    m_counters[Form("%s_nBackgroundClusters2", type)] = 0;


    // ############# Plots per station #############
    for (const auto st : GEOM.STATIONS) {
      
      //Cluster properties
      m_plots[Form("%s_x_LowBestQ_st%d", type, st)] = new TH1D(
          Form("%s_x_LowBestQ_st%d", type, st), Form("%s_x_LowBestQ_st%d; xLoc; Entries", type, st), 100, -220, 220);
      m_plots[Form("%s_xLoc_ITGhost_st%d", type, st)] =
          new TH1D(Form("%s_xLoc_ITGhost_st%d", type, st), Form("%s_xLoc_ITGhost_st%d; xLoc; Entries", type, st), 100,
                   -220, 220);
      m_plots[Form("%s_xLoc_OoTGhost_st%d", type, st)] =
          new TH1D(Form("%s_xLoc_OoTGhost_st%d", type, st), Form("%s_xLoc_OoTGhost_st%d; xLoc; Entries", type, st), 100,
                   -220, 220);

      m_2Dplots[Form("%s_N_Cluster_st%d", type, st)] = new TH2I(
          Form("%s_N_Cluster_st%d", type, st), Form("%s_N_Cluster_st%d", type, st), 14, -0.5, 13.5, 7, -3.5, 3.5);
      m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st)] =
          new TH2D(Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st),
                   Form("AllTPs_LocalDirectionvsPosition_st%d", st), 200, -200, 200, 11, -1, 1);

      m_effs[Form("%s_Eff_SegMatch_st%d", type, st)] =
          new TEfficiency(Form("%s_Eff_SegMatch_st%d", type, st), Form("%s_Eff_SegMatch_st%d; sector; wheel", type, st),
                          14, -0.5, 13.5, 7, -3.5, 3.5);
      m_effs[Form("%s_Eff_DigiMatch_st%d", type, st)] =
          new TEfficiency(Form("%s_Eff_DigiMatch_st%d", type, st),
                          Form("%s_Eff_DigiMatch_st%d; sector; wheel", type, st), 14, -0.5, 13.5, 7, -3.5, 3.5);

      //background
      m_counters[Form("%s_bkgCluster_st%d", type, st)] = 0;
      m_counters[Form("%s_matchedCluster_st%d", type, st)] = 0;

      m_plots[Form("%s_BackgroundResidual_st%d", type, st)] =
          new TH1D(Form("%s_BackgroundResidual_st%d", type, st),
                   Form("%s_BackgroundResidual_st%d; bkg.xLoc - signal.xLoc; Entries", type, st), 100, -220, 220);
      m_plots[Form("%s_q8position_st%d", type, st)] = new TH1D(
          Form("%s_q8position_st%d", type, st), Form("%s_q8position_st%d; xLoc; Entries", type, st), 100, -220, 220);


      // ############# Plots per wheel  #############
      for (const auto wh : GEOM.WHEELS) {
        m_plots[Form("%s_Digi_residual_st%d_wh%d", type, st, wh)] = new TH1D(
            Form("%s_Digi_residual_st%d_wh%d", type, st, wh),
            Form("%s_Digi_residual_st%d_wh%d; cluster.bestSeg.xLoc - digi.xLoc; entries", type, st, wh), 51, -10, 10);
      }


      // ############# Plots per quality #############
      for (auto q : QUAL_PLOT) {

        //Cluster properties
        //Matching efficiency 
        m_effs[Form("%s_ClusterEfficiencyVSpos_st%d_minqual%i", type, st, q)] = new TEfficiency(
            Form("%s_ClusterEfficiencyVSpos_st%d_minqual%d", type, st, q),
            Form("%s_ClusterEfficiencyVSpos_st%d_minqual%d; sector; wheel", type, st, q), 14, -0.5, 13.5, 7, -3.5, 3.5);
        m_effs[Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d", type, st, q),
            Form("%s_ClusterEfficiencyVSpT_st%d_minqual%d; pT (Gev); Efficiency", type, st, q), 50, MU::MIN_PT, 100);
        m_effs[Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d", type, st, q)] =
            new TEfficiency(Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d", type, st, q),
                            Form("%s_ClusterEfficiencyVSphi_st%d_minqual%d; phi (rad); Efficiency", type, st, q), 40,
                            -TMath::Pi(), TMath::Pi());
        m_effs[Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d", type, st, q)] =
            new TEfficiency(Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d", type, st, q),
                            Form("%s_ClusterEfficiencyVSeta_st%d_minqual%d; eta; Efficiency", type, st, q), 50, -1, 1);
        m_effs[Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d", type, st, q),
            Form("%s_ClusterEfficiencyVSsegxLoc_st%d_minqual%d; bestSeg_xLoc (cm); Efficiency", type, st, q), 100, -220,
            220);
        m_effs[Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d", type, st, q),
            Form("%s_ClusterEfficiencyVSsegDirLoc_st%d_minqual%d; bestSeg_dirLoc; Efficieny", type, st, q), 40,
            -TMath::Pi(), TMath::Pi());

        //Ghost distributions
        m_effs[Form("%s_GhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_GhostFractionVSpT_st%d_minqual%d", type, st, q),
            Form("%s_GhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with ghosts", type, st, q), 50,
            MU::MIN_PT, 100);
        m_effs[Form("%s_ITGhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_ITGhostFractionVSpT_st%d_minqual%d", type, st, q),
            Form("%s_ITGhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with in-time ghosts", type, st, q),
            50, MU::MIN_PT, 100);
        m_effs[Form("%s_OOTGhostFractionVSpT_st%d_minqual%d", type, st, q)] = new TEfficiency(
            Form("%s_OOTGhostFractionVSpT_st%d_minqual%d", type, st, q),
            Form("%s_OOTGhostFractionVSpT_st%d_minqual%d; pT(GeV); Fraction of events with out-of-time ghosts", type,
                 st, q),
            50, MU::MIN_PT, 100);

        m_plots[Form("%s_GhostDistributionVSpT_st%d_minqual%d", type, st, q)] =
            new TProfile(Form("%s_GhostDistributionVSpT_st%d_minqual%d", type, st, q),
                         Form("%s_GhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # ghosts", type, st, q), 50,
                         MU::MIN_PT, 100);
        m_plots[Form("%s_ITGhostDistributionVSpT_st%d_minqual%d", type, st, q)] = new TProfile(
            Form("%s_ITGhostDistributionVSpT_st%d_minqual%d", type, st, q),
            Form("%s_ITGhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # in-time ghosts", type, st, q), 50,
            MU::MIN_PT, 100);
        m_plots[Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d", type, st, q)] = new TProfile(
            Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d", type, st, q),
            Form("%s_OOTGhostDistributionVSpT_st%d_minqual%d; pT(GeV); Average # out-of-time ghosts", type, st, q), 50,
            MU::MIN_PT, 100);
        
        //background
        m_plots[Form("%s_ITBackground_st%d_minqual%d", type, st, q)] =
            new TH1F(Form("%s_ITBackground_st%d_minqual%d", type, st, q),
                     Form("%s_ITBackground_st%d_minqual%d; xLoc (cm); Entries", type, st, q), 100, -220, 220);
        m_plots[Form("%s_OOTBackground_st%d_minqual%d", type, st, q)] =
            new TH1F(Form("%s_OOTBackground_st%d_minqual%d", type, st, q),
                     Form("%s_OOTBackground_st%d_minqual%d; xLoc (cm); Entries", type, st, q), 100, -220, 220);
      }
    }
  }
}

void Analyser::ClusterAnalysis(const std::vector<Cluster> &clusters, const std::string &typeStr,
                               const std::vector<Segment> &segments) {
  for (auto const &cluster : clusters) {
    auto wh{cluster.wheel};
    auto sec{cluster.sector};
    assert(sec < 13);  // if sector == 13 o 14 the program crashes
    auto st{cluster.station};
    const auto type{typeStr.c_str()};

    TriggerPrimitive bestTP = cluster.bestTP();
    Segment bestSeg = cluster.bestSeg();

    m_2Dplots[Form("%s_TPs_LocalDirectionvsPosition_st%d", type, st)]->Fill(bestTP.xLoc, bestTP.psi);
    m_plots[Form("%s_ClusterSize", type)]->Fill(cluster.tpClusterSize());
    
    // ############# Study TP timing #############
    m_plots["FirstTpBX"]->Fill(cluster.earliestTPBX()); 
    m_plots["FirstTpt0"]->Fill(cluster.earliestTPt0());

    // ############# Study TP ghost distribution #############
    if (cluster.foundTP) ++m_counters[Form("%s_nClusters", type)];
    m_counters[Form("%s_ooTHQCount", type)] +=
        cluster.ootCountIf([=](TriggerPrimitive const &tp) { return tp.quality > HIGH_QUAL_CUT; });

    int bestQ{cluster.bestTPQuality()};
    int ootSize{cluster.ootSize()};
    int itSize{cluster.itSize()};

    m_plots[Form("%s_Q_Best", type)]->Fill(bestQ);
    m_plots[Form("%s_ITGhosts", type)]->Fill(itSize);
    m_plots[Form("%s_OOTGhosts", type)]->Fill(ootSize);
    m_plots[Form("%s_N_Ghosts", type)]->Fill(ootSize + itSize);
    m_counters[Form("%s_nTPs", type)] += cluster.tpClusterSize();

    if (bestQ == 1) m_plots[Form("%s_x_LowBestQ_st%d", type, st)]->Fill(bestTP.xLoc);

    m_2Dplots[Form("%s_N_Cluster_st%d", type, st)]->Fill(sec, wh);

    if (cluster.hasGhosts()) {
      ++m_counters[Form("%s_nClustersGhosts", type)];

      for (const auto &ghost : cluster.ootGhosts()) {
        m_plots[Form("%s_Res_OoTGhosts", type)]->Fill(ghost.xLoc - bestTP.xLoc);
        m_2Dplots[Form("%s_Q_OoTGhosts", type)]->Fill(bestQ, ghost.quality);
        m_plots[Form("%s_Q_Ghost", type)]->Fill(ghost.quality);
        m_plots[Form("%s_xLoc_OoTGhost_st%d", type, st)]->Fill(ghost.xLoc);

      }

      for (const auto &ghost : cluster.itGhosts()) {
        m_plots[Form("%s_Res_ITGhosts", type)]->Fill(ghost.xLoc - bestTP.xLoc);
        m_2Dplots[Form("%s_Q_ITGhosts", type)]->Fill(bestQ, ghost.quality);
        m_plots[Form("%s_Q_Ghost", type)]->Fill(ghost.quality);
        m_plots[Form("%s_xLoc_ITGhost_st%d", type, st)]->Fill(ghost.xLoc);

      }
    }

    if (cluster.muMatched) ++m_counters[Form("%s_matchedCluster_st%d", type, st)];

    // ########## Study efficiency #############
    if (cluster.muMatched && bestSeg.nPhiHits >= 4 && cluster.xEdge < EFF::XEDGE_CUT && cluster.yEdge < EFF::YEDGE_CUT) {
      int muon = cluster.MuIndex();
      double muEta = mu_eta->at(muon);
      double muPt = mu_pt->at(muon);
      if (std::abs(muEta) > EFF::MAX_ETA || muPt < EFF::MIN_PT) continue;
      double bestSegIndex = cluster.bestSegIndex();
      double dirLoc = std::atan(seg_dirLoc_x->at(bestSegIndex) / seg_dirLoc_z->at(bestSegIndex));
      for (auto qual : QUAL_PLOT) {
        FillEfficiency2D(type, "pos", st, qual, sec, wh, cluster);
        FillEfficiency(type, "pT", st, qual, muPt, cluster);
        FillEfficiency(type, "phi", st, qual, mu_phi->at(muon), cluster);
        FillEfficiency(type, "eta", st, qual, muEta, cluster);
        FillEfficiency(type, "segxLoc", st, qual, bestSeg.xLoc, cluster);
        FillEfficiency(type, "segDirLoc", st, qual, dirLoc, cluster);
        FillGhostRatio(type, "pT", st, qual, muPt, cluster);
        FillGhostProfile(type, "pT", st, qual, muPt, cluster);
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
  }
}

void Analyser::BackgroundAnalysis(const std::vector<Cluster> &clusters, const std::string &typeStr) {
  for (auto const &cluster : clusters) {
    if (cluster.muMatched) continue;
    if (!cluster.foundTP) continue;

    auto wh{cluster.wheel};
    auto sec{cluster.sector};
    assert(sec < 13);  // if sector == 13 o 14 the program crashes
    auto st{cluster.station};
    const auto type{typeStr.c_str()};
    double xLoc = cluster.bestTP().xLoc;

    for (auto qual : QUAL_PLOT) {
      FillBackground(type, st, qual, xLoc, cluster);
    }

    ++m_counters[Form("%s_nBackgroundClusters2", type)];
    ++m_counters[Form("%s_bkgCluster_st%d", type, st)];

    m_plots[Form("%s_Q_Bkg", type)]->Fill(cluster.bestTPQuality());

    if (cluster.bestTPQuality()==8) m_plots[Form("%s_q8position_st%d", type, st)]->Fill(xLoc);
    m_plots[Form("%s_N_bkg", type)] ->Fill(cluster.tpClusterSize());

    for (auto const &otherCluster : clusters) {
      if (otherCluster.wheel != wh || otherCluster.station != st || otherCluster.sector != sec ||
          otherCluster.muMatched == false)
        continue;

      m_plots[Form("%s_BackgroundResidual_st%d", type, st)]->Fill(xLoc - otherCluster.bestTP().xLoc);
    }
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
  const char *out  = "results/outputFile_";
  const auto Ntuplename = std::filesystem::path{ AnalyserBase::filename }.filename().string(); 
  std::string outputFileName=out+Ntuplename;
  std::cout << "Writing output in " << outputFileName << std::endl;
  TFile outputFile(outputFileName.c_str(), "RECREATE");  // (outputFileName.c_str(), "RECREATE");
  outputFile.cd();

  DefinePlot();

  double removedHQFilter{0.0};
  double removedLQFilter{0};

  if (fChain == 0) return;

  Long64_t n_entries{fChain->GetEntriesFast()};

  for (Long64_t entry{0}; entry < n_entries; ++entry) {
    Long64_t i_entry{LoadTree(entry)};
    if (i_entry < 0) break;
    fChain->GetEntry(entry);

    // ########## LOOP ON EVENTS #############
    if (entry % 100 == 0) std::cout << "Processing event: " << entry << '\r' << std::flush;

    if (gen_pdgId->size() < 1 || gen_pt->at(0) < 10.0) continue;

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
      digis.emplace_back(Digi(i, digi_wheel->at(i), digi_sector->at(i), digi_station->at(i), digi_superLayer->at(i),
                              digi_layer->at(i), digi_wire->at(i), digi_time->at(i)));
    }

    // ########## BUILD segments std::vector #############
    std::vector<Segment> segments;
    for (std::size_t j = 0; j < seg_nSegments; ++j) {
      segments.emplace_back(Segment(j, seg_station->at(j), seg_wheel->at(j), seg_sector->at(j), seg_phi_nHits->at(j),
                                    seg_posLoc_x->at(j)));
    }

    // ########## BUILD clusters std::vector #############
    auto clusters = buildClusters(tps, segments, digis, X_CUT, DIGI_CUT);
    int nClusters = clusters.size();

    m_plots["ClusterPerEvent"]->Fill(nClusters);

    // ########## BUILD cluster with phi matching information #############
    // tag TPs that match using the extrapolation on a straight line

    auto match_tps = [&](const std::vector<int> &quals) -> std::vector<TriggerPrimitive> {
      std::vector<TriggerPrimitive> matched_tps = tps;

      for (auto &tp : matched_tps) {
        if (!tp.hasMatched) {
          for (auto &other_tp : matched_tps) {
            if (tp.index != other_tp.index) {
              tp.Match(other_tp, PHI_CUT, T0_CUT, quals);
            }
          }
        }
      }

      matched_tps.erase(std::remove_if(matched_tps.begin(), matched_tps.end(),
                                       [](auto &tp) { return !tp.hasMatched && tp.quality == 1; }),
                        matched_tps.end());

      return matched_tps;
    };

    // extrapolate from LQ TPs and repeat clustering
    auto matchFromLQ_tps = match_tps({1});
    auto matchFromLQ_clusters = buildClusters(matchFromLQ_tps, segments, digis, X_CUT, DIGI_CUT);

    // extrapolate from HQ TPs and repeat clustering
    auto matchFromHQ_tps = match_tps({3, 4, 5, 6, 7, 8, 9});
    auto matchFromHQ_clusters = buildClusters(matchFromHQ_tps, segments, digis, X_CUT, DIGI_CUT);

    // ########## Muon ID and GEN matching #############
    std::vector<UInt_t> goodMuons;
    std::vector<UInt_t> unmatchedGen;
    for (UInt_t igen = 0; igen < gen_nGenParts; ++igen) {
      bool matched{false};
      double genEta{gen_eta->at(igen)};
      double genPhi{gen_phi->at(igen)};
      double genPT{gen_pt->at(igen)};
      if (std::abs(gen_pdgId->at(igen)) != 13 || std::abs(genEta) > MU::MAX_ETA || genPT < MU::MIN_PT) continue;
      TLorentzVector genMu;
      genMu.SetPtEtaPhiM(genPT, genEta, genPhi, MU::MASS);

      for (UInt_t iMu = 0; iMu < mu_nMuons; ++iMu) {
        if (!mu_isMedium->at(iMu)) continue;
        TLorentzVector mu;
        double muEta{mu_eta->at(iMu)};
        double muPhi{mu_phi->at(iMu)};
        double muPT{mu_pt->at(iMu)};
        mu.SetPtEtaPhiM(muPT, muEta, muPhi, MU::MASS);
        // check if the muon can be matched to one of the particles generated in this event

        double dR = genMu.DeltaR(mu);
        m_plots["DeltaR"]->Fill(dR);
        if (std::abs(dR) < MU::MAX_GEN_DR) {
          goodMuons.push_back(iMu);
          matched = true;
          break;
        }
      }
      if (!matched) {
        unmatchedGen.push_back(igen);
        m_2Dplots["Pos_unmatchedGenMu"]->Fill(gen_eta->at(igen), gen_phi->at(igen));
        m_plots["Pt_unmatchedGenMu"]->Fill(gen_pt->at(igen));
      }

    }

    int goodMuSize = goodMuons.size();
    m_plots["NGoodMu"]->Fill(goodMuSize);
    if (goodMuSize > 0)
      m_plots["GoodMuClusterPerEvent"]->Fill(nClusters);
    else if (goodMuSize < 1 && nClusters > 0) {
      int nTPS = 0;
      for (const auto &cl : clusters) nTPS += cl.tpClusterSize();
      // std::cout << "event with no good muons but "<< nClusters << " clusters with " << nTPS<< std::endl;
      m_counters["noGoodMuTPS"] += nTPS;
    }
    
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
            ++stationMatch[muStation - 1];
            m_2Dplots[Form("PosGoodMu_st%d", muStation)]->Fill(muEta, muPhi);
            if (!cluster.foundTP) std::cout << " Matchato un cluster senza best TP " << std::endl;
          };

          if (wh == muWheel && sec == muSector && st == muStation)
            m_plots[Form("Res_Mu_st%d", st)]->Fill(cluster.bestTP().xLoc - muX);

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

      for (int st = 1; st < 5; ++st) {
        m_2Dprofiles[Form("NMatchGoodMu_st%d", st)]->Fill(muEta, muPhi, stationMatch[st - 1], 1);
      }
    }

    // ########## RUN SOME ANALYSIS #############
    // ########## CLUSTER ANALYSIS #############

    ClusterAnalysis(clusters, tags[0], segments);
    ClusterAnalysis(matchFromLQ_clusters, tags[1], segments);
    ClusterAnalysis(matchFromHQ_clusters, tags[2], segments);

    // ########## Study background #############
    BackgroundAnalysis(clusters, tags[0]);
    BackgroundAnalysis(matchFromLQ_clusters, tags[1]);
    BackgroundAnalysis(matchFromHQ_clusters, tags[2]);

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
    double bkgFraction = m_counters[Form("%s_nBackgroundClusters2", type)] / m_counters[Form("%s_nClusters", type)];
    std::cout << " \nFor " << tag << std::endl;

    std::cout << " Fraction of clusters with ghost (" << m_counters[Form("%s_nClustersGhosts", type)] << ") on total ("
              << m_counters[Form("%s_nClusters", type)] << ") = " << ghostFraction << " made with "
              << m_counters[Form("%s_nTPs", type)] << " TPs" << std::endl;
    std::cout << " HQ out of time clusters: " << m_counters[Form("%s_ooTHQCount", type)] << std::endl;

    std::cout << " Fraction of background cluster (" << m_counters[Form("%s_nBackgroundClusters", type)] << " / "
              << m_counters[Form("%s_nBackgroundClusters2", type)] << ") on total ("
              << m_counters[Form("%s_nClusters", type)] << ") = " << bkgFraction << std::endl;
    
    std::cout << " Fraction of cluster with prefiring TPs (" << m_plots[Form("FirstTpBX")]->GetBinContent(5) << ") with respect to the right BX ("  << m_plots[Form("FirstTpBX")]->GetBinContent(6) << ") = " 
              <<   m_plots[Form("FirstTpBX")]->GetBinContent(5)/m_plots[Form("FirstTpBX")]->GetBinContent(6) << std::endl;
    
    for (const auto st : GEOM.STATIONS) {
      double ratio = m_counters[Form("%s_bkgCluster_st%d", type, st)] /  m_counters[Form("%s_matchedCluster_st%d", type, st)];
      std::cout << "In station " << st << " the fraction of bkg clusters (" << m_counters[Form("%s_bkgCluster_st%d", type, st)] << ") with respect to matched clusters (" 
                << m_counters[Form("%s_matchedCluster_st%d", type, st)] << ") is " << ratio << std::endl;
    }

  }

  std::cout << "LQFilter found " << removedLQFilter << " less clusters, HQFilter found " << removedHQFilter
            << " less clusters " << std::endl;

  std::cout << "TPs in bkg because of no muon in event: " << m_counters["noGoodMuTPS"] << std::endl;

  outputFile.Write();
  outputFile.Close();
}

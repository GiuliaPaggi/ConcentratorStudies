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

// CB not optimal, but readable
const std::vector<int> WHEELS{-2, -1, 0, 1, 2};
const std::vector<int> SECTORS{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
const std::vector<int> STATIONS{1, 2, 3, 4};

std::vector<Cluster> buildClusters(std::vector<TriggerPrimitive> tps, std::vector<Segment> seg, double x_cut) {
  std::vector<Cluster> clusters;

  for (const auto wh : WHEELS) {
    for (const auto sec : SECTORS) {
      for (const auto st : STATIONS) {
        Cluster cluster{tps, seg, x_cut, wh, sec, st};   // aggiungi al costruttore il segmento se c'è -> match direttamente nel costruttore
        if (cluster.bestTPQuality() > -1 || cluster.bestSegPhiHits() > -1 ) {
          clusters.push_back(cluster);  // CB can be improved 
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


void Analiser::Loop() {
  const std::array<double, 4> MB{402.2, 490.5, 597.5, 700.0};
  const int LOW_QUAL_CUT{0};
  const int HIGH_QUAL_CUT{5};
  const double PSI_CUT{TMath::Pi() / 6.0};
  double PHI_CUT{0.02};
  const double PHI_CUT_2{0.01};
  const double T0_CUT{12.5};
  const double X_CUT{5.0};
  const int CORRECT_BX{380};

  double T_MIN{-9800};
  double T_MAX{-9200};
  double BX_MIN{-392};
  double BX_MAX{-368};

  /*int counter = 0;
  int counter_bad = 0;
  int counter_good = 0;
  */

  TFile *file = new TFile("./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root");   //  ./VtxSmeared/DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25_VtxSmeared.root
//./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root
  TH1D *t0_AllQuality = new TH1D("t0_AllQuality", "t0_AllQuality", 100, T_MIN, T_MAX);
  TH1D *t0_HighQuality =
      new TH1D("t0_HighQuality", "t0_HighQuality", 100, T_MIN, T_MAX);  // 3+3 4+4 3+4
  TH1D *t0_IntermediateQuality =
      new TH1D("t0_IntermediateQuality", "t0_IntermediateQuality", 100, T_MIN,
               T_MAX);  // 4+2 3+2 qualities, not there if we use slice-test configuration for emulator
  TH1D *t0_LowQuality = new TH1D("t0_LowQuality", "t0_LowQuality", 100, T_MIN, T_MAX);
  TH1D *t0_Selected = new TH1D("t0_Selected", "t0_Selected", 100, T_MIN, T_MAX);

  TH1D *LowQ_matched = new TH1D("LowQ_matched", "LowQ_matched", 100, T_MIN, T_MAX);
  TH1D *HighQ_matched = new TH1D("HighQ_matched", "HighQ_matched", 100, T_MIN, T_MAX);

  TH1D *BX_LowQuality = new TH1D("BX_LowQuality", "BX_LowQuality,BX;", 24, BX_MIN, BX_MAX);
  TH1D *BX_LowQ_matched =
      new TH1D("BX_LowQ_matched", "BX_LowQ_matched,BX;Entries", 24, BX_MIN, BX_MAX);
  TH1D *BX_LowQ_more1HQ =
      new TH1D("BX_LowQ_more1HQ", "BX_LowQ_more1HQ,BX;Entries", 24, BX_MIN, BX_MAX);

  TH1D *t0_Selected_Psi = new TH1D("t0_Selected_Psi", "t0_Selected_Psi", 100, T_MIN, T_MAX);
  TH1D *LowQ_matched_Psi = new TH1D("LowQ_matched_Psi", "LowQ_matched_Psi", 100, T_MIN, T_MAX);
  TH1D *HighQ_matched_Psi = new TH1D("HighQ_matched_Psi", "HighQ_matched_Psi", 100, T_MIN, T_MAX);

  TH1D *LowQ_more1HQ = new TH1D("LowQ_more1HQ", "LowQ_more1HQ", 100, T_MIN, T_MAX);
  TH1D *LowQ_more1HQ_Phi = new TH1D("LowQ_more1HQ_Phi", "LowQ_more1HQ_Phi", 100, T_MIN, T_MAX);

  TH2D *PhiRes_st1_2 = new TH2D("PhiRes_st1_2", "PhiRes_st1_2;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);
  TH2D *PhiRes_st1_3 = new TH2D("PhiRes_st1_3", "PhiRes_st1_3;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);
  TH2D *PhiRes_st1_4 = new TH2D("PhiRes_st1_4", "PhiRes_st1_4;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);
  TH2D *PhiRes_st2_3 = new TH2D("PhiRes_st2_3", "PhiRes_st2_3;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);
  TH2D *PhiRes_st2_4 = new TH2D("PhiRes_st2_4", "PhiRes_st2_4;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);
  TH2D *PhiRes_st3_4 = new TH2D("PhiRes_st3_4", "PhiRes_st3_4;p_{T} (GeV);Computed Phi -Phi (rad)",
                                100, 15, 105, 100, -1, 1);

  TH2D *PhiRes_vs_posizione =
      new TH2D("PhiRes_vs_posizione",
               "PhiRes_vs_posizione;Vertex distance from IP (cm);Computed Phi -Phi (rad)", 100, 0,
               60, 500, -.1, .1);
  TH2D *PhiRes_vs_dxy = new TH2D(
      "PhiRes_vs_dxy", "PhiRes_vs_dxy;Vertex distance from IP (cm);Computed Phi -Phi (rad)", 100, 0,
      35, 500, -.1, .1);
  TH2D *PhiRes_vs_dz =
      new TH2D("PhiRes_vs_dz", "PhiRes_vs_dz;Vertex distance from IP (cm);Computed Phi -Phi (rad)",
               100, -45, 45, 500, -.1, .1);

  TEfficiency *Match_vs_Eta = new TEfficiency("Match_vs_Eta", "Match_vs_Eta", 24, 0, 1);
  TEfficiency *Match_vs_Phi = new TEfficiency("Match_vs_Phi", "Match_vs_Phi", 20, -TMath::Pi(), TMath::Pi());
  TEfficiency *Match_vs_Pt = new TEfficiency("Match_vs_Pt", "Match_vs_Pt", 40, 20, 100);

  TEfficiency *MatchAndCut_vs_Eta =
      new TEfficiency("MatchAndCut_vs_Eta", "MatchPhi_vs_Eta", 24, 0, 1);
  TEfficiency *MatchAndCut_vs_Phi =
      new TEfficiency("MatchAndCut_vs_Phi", "MatchPhi_vs_Phi", 20, -TMath::Pi(), TMath::Pi());
  TEfficiency *MatchAndCut_vs_Pt =
      new TEfficiency("MatchAndCut_vs_Pt", "MatchPhi_vs_Pt", 40, 20, 100);

  // GHOST HISTOS
  // TH1D OoTGhosts_style("OoTGhosts","OoTGhosts_style", 10, 0, 10);
  // TH1D *OoTGhosts[5][4][12];

  // TH1D ITGhosts_style("ITGhosts", "ITGhosts_style", 10, 0, 10);
  // TH1D *ITGhosts[5][4][12];

  TH1D *OoTGhosts = new TH1D("OoTGhosts", "OoTGhosts", 20, 0, 20);
  TH1D *ITGhosts = new TH1D("ITGhosts", "ITGhosts", 20, 0, 20);

  TH1I *BX_ITGhosts = new TH1I("BX_ITGhosts", "BX_ITGhosts", 24, BX_MIN, BX_MAX);
  TH1I *BX_OoTGhosts = new TH1I("BX_OoTGhosts", "BX_OoTGhosts", 24, BX_MIN, BX_MAX);
  TH1D *Res_ITGhosts = new TH1D("Res_ITGhosts", "Res_ITGhosts", 110, -5.5, 5.5);
  TH1D *Res_OoTGhosts = new TH1D("Res_OoTGhosts", "Res_OoTGhosts", 110, -5.5, 5.5);

  TH2D *Q_OoTGhosts = new TH2D("Q_OoTGhosts", "Q_OoTGhosts;High Quality;Out of time Ghost Quality",
                               10, 0, 10, 10, 0, 10);
  TH2D *Q_ITGhosts = new TH2D("Q_ITGhosts", "Q_ITGhosts;High Quality;In time Ghost Quality", 10, 0, 10, 10, 0, 10);

  TH1I *N_Ghost = new TH1I("N_Ghost", "N_Ghost", 20, 0, 20);
  TH1I *Q_Best = new TH1I("Q_Best", "Q_Best", 10, 0, 10);
  TH1I *Q_Ghost = new TH1I("Q_Ghost", "Q_Ghost", 10, 0, 10);

  TH2I N_Cluster_style("N_Cluster_style", "N_Cluster_style", 14, -0.5, 13.5, 7, -3.5, 3.5);
  TH2I *N_Cluster[4];
  for (const auto st : STATIONS) {
    N_Cluster[st - 1] = new TH2I(Form("N_Cluster_st%d", st), Form("N_Cluster_st%d", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    //N_Cluster[st - 1]->SetName(Form("N_Cluster_st%d", st));
  }
  TH1D *x_LowBestQ[4];

  TH1D *Res_MuMatched = new TH1D("Res_MuMatched", "Res_MuMatched; Entries", 100, -6, 6);
  TH1D *Res_SegMatched = new TH1D("Res_SegMatched", "Res_SegMatched; Entries", 100, -6, 6);
  TH2I *N_MuMatch[4];
  TH1D *Phi_MuMatch[4];
  TEfficiency *Eff_MuMatch[4];
  TH2D *NMuMatch_vs_phi[4];

  TEfficiency *Eff_SegMatch[4];

  for (const auto st : STATIONS) {
    N_MuMatch[st - 1] = new TH2I(Form("N_MuMatch_st%d", st), Form("N_MuMatch_st%d; sector ; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    Eff_MuMatch[st-1] = new TEfficiency(Form("Eff_MuMatch_st%d", st), Form("Eff_MuMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    x_LowBestQ[st-1] = new TH1D(Form("x_LowBestQ_st%d", st), Form("x_LowBestQ_st%d; xLoc; Entries",st ), 100, -220, 220);
    NMuMatch_vs_phi[st-1] = new TH2D(Form("NMuMatch_vs_phi_st%d", st), Form("NMuMatch_vs_phi_st%d ; Muon_Phi; N_MuMatch", st), 50, -TMath::Pi(), TMath::Pi(), 6, 1, 7 );
    Phi_MuMatch[st-1] = new TH1D(Form("Phi_MuMatch_st%d", st), Form("Phi_MuMatch_st%d; Muon_phi; Entries", st), 50,-TMath::Pi(), TMath::Pi() );
    Eff_SegMatch[st-1] = new TEfficiency(Form("Eff_SegMatch_st%d", st), Form("Eff_SegMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
  }


   
  //   In a ROOT session, you can do:
  //      root> .L Analiser.C
  //      root> Analiser t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  double nClustersGhosts{};
  double ooTHQCount{};
  double nClusters{};

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
    // cout << "---------------------------------------" << endl;

    // ########## CREATE TPs std::vector #############
    std::vector<TriggerPrimitive> tps;

    for (std::size_t j = 0; j < ph2TpgPhiEmuAm_nTrigs; ++j) {
      tps.emplace_back(TriggerPrimitive{j, ph2TpgPhiEmuAm_wheel->at(j),
                                        ph2TpgPhiEmuAm_sector->at(j), ph2TpgPhiEmuAm_station->at(j),
                                        ph2TpgPhiEmuAm_quality->at(j), ph2TpgPhiEmuAm_phi->at(j),
                                        ph2TpgPhiEmuAm_phiB->at(j), ph2TpgPhiEmuAm_BX->at(j),
                                        ph2TpgPhiEmuAm_t0->at(j), ph2TpgPhiEmuAm_posLoc_x->at(j)});
    }

    //########## BUILD segments std::vector #############
    std::vector<Segment> segments;
    for (std::size_t j = 0; j < seg_nSegments; ++j) {
      segments.emplace_back(Segment(j, seg_station->at(j), seg_wheel->at(j), seg_sector->at(j), seg_phi_nHits->at(j), seg_posLoc_x->at(j)));
    }


    // ########## BUILD clusters std::vector #############
    auto clusters = buildClusters(tps, segments, X_CUT);


    // ########## ATTEMPT cluster - muon extrapolation matching #############
    for (auto &cluster: clusters){
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

          cluster.MatchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
        }
      }
    }

    // ########## RUN SOME ANALYSIS #############
    for (auto const &cluster : clusters) {
      auto wh{cluster.wheel};
      auto sec{cluster.sector};
      auto st{cluster.station};

      ++nClusters;
      ooTHQCount += cluster.ootCountIf([=](TriggerPrimitive const & tp) { return tp.quality > HIGH_QUAL_CUT; });

      int bestQ = cluster.bestTPQuality();
      int ootSize{cluster.ootSize()};
      int itSize{cluster.itSize()};

      Q_Best->Fill(bestQ);
      ITGhosts->Fill(itSize);
      OoTGhosts->Fill(ootSize);
      N_Ghost->Fill(ootSize + itSize);

      if (bestQ == 1) x_LowBestQ[st-1]->Fill(cluster.bestTP().xLoc);

      Eff_MuMatch[st-1]->Fill(cluster.muMatched, sec, wh);
        if (cluster.muMatched) {
          Res_MuMatched->Fill(cluster.bestTP().xLoc - getXY<float>(mu_matches_x, cluster.muMatchedIndex[0], cluster.muMatchedIndex[1]));
          N_MuMatch[st-1]->Fill(sec, wh);
      }
      
      N_Cluster[st - 1]->Fill(sec, wh);

      if (!cluster.hasGhosts()) continue;
      ++nClustersGhosts;

      for (const auto &ghost : cluster.ootGhosts()) {
        BX_OoTGhosts->Fill(ghost.BX);
        Res_OoTGhosts->Fill(ghost.xLoc - cluster.bestTP().xLoc);
        Q_OoTGhosts->Fill(bestQ, ghost.quality);
        Q_Ghost->Fill(ghost.quality);
      }

      for (const auto &ghost : cluster.itGhosts()) {
        BX_ITGhosts->Fill(ghost.BX);
        Res_ITGhosts->Fill(ghost.xLoc - cluster.bestTP().xLoc);
        Q_ITGhosts->Fill(bestQ, ghost.quality);
        Q_Ghost->Fill(ghost.quality);
      }

      if (cluster.foundTP && cluster.foundTP){
        Eff_SegMatch[st-1]->Fill(cluster.segMatched, sec, wh);
        Res_SegMatched->Fill( cluster.bestTP().xLoc - cluster.matchedSeg().xLoc );

        /*if (sec == 7 && wh == 2 && st == 2) {
          counter += 1;
          if (cluster.segMatched) counter_good+= 1;
          if (!cluster.segMatched) {
            counter_bad +=1;
            std::cout << "\n" <<cluster.bestTP().xLoc - cluster.matchedSeg().xLoc << " TPfound: " << cluster.foundTP << " SegFound: " << cluster.foundTP << std::endl;
          }
        }*/

      }


    }

    for (TriggerPrimitive &tp : tps) {
      if (tp.quality == 1) {
        t0_LowQuality->Fill(tps.back().t0);
        BX_LowQuality->Fill(tps.back().BX);
      }
      if (tp.quality > HIGH_QUAL_CUT && !tp.hasMatched) {
        // select HQ TPs which are not matched 
        t0_Selected->Fill(tp.t0);
        for (TriggerPrimitive &other_tp : tps) {
          if (tp.index != other_tp.index && tp.Match(other_tp, PHI_CUT, T0_CUT)) {    
            if (other_tp.quality == 1) {
              t0_Selected->Fill(other_tp.t0);
              LowQ_matched->Fill(other_tp.t0);
              BX_LowQ_matched->Fill(other_tp.BX);
            } else {
              t0_Selected->Fill(other_tp.t0);
            }
          }
        }
      }
    }

    for (TriggerPrimitive tp : tps) {
      if (tp.quality == 1) {
        MatchAndCut_vs_Eta->Fill(tp.hasMatched, std::abs(gen_eta->at(0)));
        MatchAndCut_vs_Phi->Fill(tp.hasMatched, std::abs(gen_phi->at(0)));
        MatchAndCut_vs_Pt->Fill(tp.hasMatched, std::abs(gen_pt->at(0)));
        Match_vs_Eta->Fill((tp.Matches.size() > 0), std::abs(gen_eta->at(0)));
        Match_vs_Phi->Fill((tp.Matches.size() > 0), std::abs(gen_phi->at(0)));
        Match_vs_Pt->Fill((tp.Matches.size() > 0), std::abs(gen_pt->at(0)));
      }
      if (tp.quality == 1 && tp.Matches.size() > 0) {
        LowQ_more1HQ_Phi->Fill(tp.t0);
        BX_LowQ_more1HQ->Fill(tp.BX);
      }
    }
  }
  double ghostFraction = nClustersGhosts / nClusters;

  cout << " Ratio LQ/selected with phi=  " << LowQ_matched->GetEntries() << "/ "
       << t0_LowQuality->GetEntries() << " = "
       << LowQ_matched->GetEntries() / t0_LowQuality->GetEntries() << endl;
  cout << " Ratio LQ/matched with phi=  " << LowQ_more1HQ_Phi->GetEntries() << "/ "
       << t0_LowQuality->GetEntries() << " = "
       << LowQ_more1HQ_Phi->GetEntries() / t0_LowQuality->GetEntries() << endl;
  cout << " Fraction of clusters with ghost (" << nClustersGhosts << ") on total (" << nClusters
       << ") = " << ghostFraction << endl;
  cout << " HQ out of time clusters: " << ooTHQCount << endl;

  //cout << "Low eff bin content: " << counter <<  " di cui non matchati " << counter_bad << " e matchati " << counter_good << endl;

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 500, 500, 500, 500);
  gPad->SetLogy();
  t0_LowQuality->SetLineColor(kBlue);
  t0_LowQuality->GetXaxis()->SetTitle(" t0 (ns)");
  t0_LowQuality->GetYaxis()->SetTitle(" Entries");
  t0_LowQuality->Draw();
  // LowQ_matched_Psi->SetLineColor(kRed);
  // LowQ_matched_Psi->Draw("same");
  LowQ_matched->SetLineColor(kRed);
  LowQ_matched->Draw("same");
  // LowQ_more1HQ->Draw("same");
  LowQ_more1HQ_Phi->SetLineColor(kGreen + 2);
  LowQ_more1HQ_Phi->Draw("same");

  auto leg = new TLegend(0.1, 0.7, 0.4, 0.9);
  leg->AddEntry(t0_LowQuality, "Q1");
  leg->AddEntry(LowQ_matched, "Q1 - #phi cut and time cut");
  // leg->AddEntry(LowQ_matched_Psi, "Q1 - with psi&time cut ");
  // leg->AddEntry(LowQ_more1HQ, "Q1 - associated HQ psi cut, no time cut ");
  leg->AddEntry(LowQ_more1HQ_Phi, "Q1 - #phi cut, no time cut");
  leg->Draw("same");

  TCanvas *canvas = new TCanvas("canvas", "canvas", 500, 500, 500, 500);
  canvas->Divide(2, 2);

  canvas->cd(1);
  gPad->SetLogy();

  t0_AllQuality->Draw();
  // t0_LowQuality->SetLineColor(kRed);
  t0_LowQuality->Draw("same");
  t0_HighQuality->SetLineColor(kBlack);
  t0_HighQuality->Draw("same");

  auto leg1 = new TLegend(0.1, 0.7, 0.35, 0.9);
  leg1->AddEntry(t0_LowQuality, "LQ (1-3) primitives");
  leg1->AddEntry(t0_HighQuality, "HQ (6-7-8) primitives");
  leg1->AddEntry(t0_AllQuality, "All primitives in #eta < 0.8");
  leg1->Draw("same");

  canvas->cd(2);
  gPad->SetLogy();
  // t0_LowQuality->SetLineColor(kBlack);
  t0_LowQuality->Draw();
  LowQ_matched_Psi->SetLineColor(kRed);
  LowQ_matched_Psi->Draw("same");
  // LowQ_matched->SetLineColor(kGreen);
  // LowQ_matched->Draw("same");
  LowQ_more1HQ->Draw("same");
  // LowQ_more1HQ_Phi->SetLineColor(kRed);
  LowQ_more1HQ_Phi->Draw("same");

  auto leg2 = new TLegend(0.1, 0.7, 0.4, 0.9);
  leg2->AddEntry(t0_LowQuality, "Q1");
  leg2->AddEntry(LowQ_matched, "Q1 - expected phi&time cut");
  leg2->AddEntry(LowQ_matched_Psi, "Q1 - with psi&time cut ");
  leg2->AddEntry(LowQ_more1HQ, "Q1 - associated HQ psi cut, no time cut ");
  leg2->AddEntry(LowQ_more1HQ_Phi, "Q1 - associated HQ phi cut, no time cut");
  leg2->Draw("same");

  canvas->cd(3);
  gPad->SetLogy();

  t0_Selected_Psi->Draw();
  HighQ_matched_Psi->SetLineColor(kGreen);
  HighQ_matched_Psi->Draw("same");
  LowQ_matched_Psi->SetLineColor(kRed);
  LowQ_matched_Psi->Draw("same");

  auto leg3 = new TLegend(0.1, 0.7, 0.35, 0.9);
  leg3->AddEntry(LowQ_matched_Psi, "Q1 primitives associated");
  leg3->AddEntry(HighQ_matched_Psi, "HQ primitives associated");
  leg3->AddEntry(t0_Selected_Psi, "All primitives");
  leg3->Draw("same");

  canvas->cd(4);
  gPad->SetLogy();

  t0_Selected->Draw();
  HighQ_matched->SetLineColor(kRed);
  HighQ_matched->Draw("same");
  // LowQ_matched->SetLineColor(kGreen);
  LowQ_matched->Draw("same");

  auto leg4 = new TLegend(0.1, 0.7, 0.35, 0.9);
  leg4->AddEntry(LowQ_matched, "Q1 primitives associated");
  leg4->AddEntry(HighQ_matched, "HQ primitives associated");
  leg4->AddEntry(t0_Selected, "All primitives");
  leg4->Draw("same");

  TCanvas *BXCanvas = new TCanvas("BXCanvas", "BXCanvas", 500, 500, 500, 500);
  BX_LowQuality->Draw();
  BX_LowQ_matched->SetLineColor(kRed);
  BX_LowQ_matched->Draw("same");
  BX_LowQ_more1HQ->SetLineColor(kGreen + 2);
  BX_LowQ_more1HQ->Draw("same");

  auto *legend = new TLegend(0.1, 0.7, 0.35, 0.9);
  legend->AddEntry(BX_LowQuality, "Q1 primitive");
  legend->AddEntry(BX_LowQ_matched, "Q1 selected #phi match and time cut");
  legend->AddEntry(BX_LowQ_more1HQ, "Q1 selected #phi match ");
  legend->Draw("same");

  TCanvas *EffCutCanvas = new TCanvas("EffCutCanvas", "EffCutCanvas", 500, 500, 500, 500);
  EffCutCanvas->Divide(2, 2);

  EffCutCanvas->cd(1);
  Match_vs_Phi->Draw();
  MatchAndCut_vs_Phi->SetLineColor(kRed);
  MatchAndCut_vs_Phi->Draw("same");

  EffCutCanvas->cd(2);
  Match_vs_Eta->Draw();
  MatchAndCut_vs_Eta->SetLineColor(kRed);
  MatchAndCut_vs_Eta->Draw("same");

  EffCutCanvas->cd(3);
  Match_vs_Pt->Draw();
  MatchAndCut_vs_Pt->SetLineColor(kRed);
  MatchAndCut_vs_Pt->Draw("same");

  TCanvas *ClusterProvaCanvas =
      new TCanvas("ClusterProvaCanvas", "ClusterProvaCanvas", 500, 500, 500, 500);
  ClusterProvaCanvas->Divide(2, 2);
  ClusterProvaCanvas->cd(1);
  OoTGhosts->Draw();

  ClusterProvaCanvas->cd(2);
  ITGhosts->Draw();

  ClusterProvaCanvas->cd(3);
  N_Ghost->Draw();

  ClusterProvaCanvas->cd(4);
  Q_Ghost->SetLineColor(kOrange + 7);
  Q_Ghost->Draw();
  Q_Best->Draw("same");
  auto *q_legend = new TLegend(0.1, 0.7, 0.35, 0.9);
  q_legend->AddEntry(Q_Ghost, "Ghost Quality");
  q_legend->AddEntry(Q_Best, "Best one Quality");
  q_legend->Draw("same");

  TCanvas *ClusterBXCanvas = new TCanvas("ClusterBXCanvas", "ClusterBXCanvas", 500, 500, 500, 500);
  ClusterBXCanvas->Divide(2, 2);

  ClusterBXCanvas->cd(1);
  BX_ITGhosts->SetLineColor(kRed);
  BX_ITGhosts->Draw();
  BX_OoTGhosts->Draw("same");

  ClusterBXCanvas->cd(2);

  Res_ITGhosts->SetLineColor(kRed);
  Res_ITGhosts->Draw();
  Res_OoTGhosts->Draw("same");

  ClusterBXCanvas->cd(3);
  Q_OoTGhosts->Draw("box");

  ClusterBXCanvas->cd(4);
  Q_ITGhosts->Draw("box");

  TCanvas *StatCanvas = new TCanvas("StatCanvas", "StatCanvas", 500, 500, 500, 500);
  StatCanvas->Divide(2, 2);

  for (int i = 1; i < 5; ++i) {
    StatCanvas->cd(i);
    N_Cluster[i - 1]->Draw("COLZ");
  }

  TCanvas *ResMatchCanvas = new TCanvas("ResMatchCanvas", "ResMatchCanvas", 500, 500, 500, 500);
  ResMatchCanvas->Divide(1,2);
  ResMatchCanvas->cd(1);
  Res_MuMatched->Draw();
  ResMatchCanvas->cd(2);
  Res_SegMatched->Draw();

  TCanvas *MuMatch2dCanvas = new TCanvas("MuMatch2dCanvas", "MuMatch2dCanvas", 500, 500, 500, 500);
  MuMatch2dCanvas->Divide(2,2);
  for (int i = 1; i < 5; ++i) {
    MuMatch2dCanvas->cd(i);
    N_MuMatch[i - 1]->Draw("COLZ");
  }

  TCanvas *EffMuMatchCanvas = new TCanvas("EffMuMatchCanvas", "EffMuMatchCanvas", 500, 500, 500, 500);
  EffMuMatchCanvas->Divide(2,2);
  for (int i = 1; i < 5; ++i) {
    EffMuMatchCanvas->cd(i);
    Eff_MuMatch[i - 1]->Draw("COLZ");
  }

  TCanvas *LowQBestDistributionCanvas = new TCanvas("LowQBestDistributionCanvas", "LowQBestDistributionCanvas", 500, 500, 500, 500);
  LowQBestDistributionCanvas->Divide(2,2);
  for(int i = 1; i< 5; ++i ){
    LowQBestDistributionCanvas->cd(i);
    x_LowBestQ[i-1]->Draw();
  }

  TCanvas *SegMatchCanvas = new TCanvas("SegMatchCanvas", "SegMatchCanvas", 500, 500, 500, 500);
  SegMatchCanvas->Divide(2, 2);
  for (int i = 1; i < 5; ++i) {
    SegMatchCanvas->cd(i);
    Eff_SegMatch[i - 1]->Draw("COLZ");
  }


}

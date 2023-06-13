#define Analyser_cxx
#include "include/TestAnalyser.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TClonesArray.h>

#include <iostream>
#include <vector>

void TestAnalyser::Loop() {

  Geometry geom{};

  // const std::array<double, 4> MB{402.2, 490.5, 597.5, 700.0}; // CB unused, but could go in Geometry
  // const int LOW_QUAL_CUT{0}; // CB unused
  const int HIGH_QUAL_CUT{5};
  // const double PSI_CUT{TMath::Pi() / 6.0}; // CB unused
  double PHI_CUT{0.02};
  // const double PHI_CUT_2{0.01}; // CB unused
  const double T0_CUT{12.5};
  const double X_CUT{5.0};
  const double DIGI_CUT{10.0};
  // const int CORRECT_BX{380}; // CB unused

  double T_MIN{-9800};
  double T_MAX{-9200};
  double BX_MIN{-392};
  double BX_MAX{-368};


  // CB unused
  // TFile *file = new TFile("./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root");   
  //  ./VtxSmeared/DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25_VtxSmeared.root
  //./DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root
  
  TFile outputFile("outputFile.root","RECREATE");
  
  TH1D *t0_AllQuality = new TH1D("t0_AllQuality", "t0_AllQuality", 100, T_MIN, T_MAX);
  TH1D *t0_HighQuality =
      new TH1D("t0_HighQuality", "t0_HighQuality", 100, T_MIN, T_MAX);  // 3+3 4+4 3+4
  // CB unused
  // TH1D *t0_IntermediateQuality =
  //    new TH1D("t0_IntermediateQuality", "t0_IntermediateQuality", 100, T_MIN,
  //             T_MAX);  // 4+2 3+2 qualities, not there if we use slice-test configuration for emulator
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

  TEfficiency *Match_vs_Eta = new TEfficiency("Match_vs_Eta", "Match_vs_Eta", 24, 0, 1);
  TEfficiency *Match_vs_Phi = new TEfficiency("Match_vs_Phi", "Match_vs_Phi", 20, -TMath::Pi(), TMath::Pi());
  TEfficiency *Match_vs_Pt = new TEfficiency("Match_vs_Pt", "Match_vs_Pt", 40, 20, 100);

  TEfficiency *MatchAndCut_vs_Eta =
      new TEfficiency("MatchAndCut_vs_Eta", "MatchPhi_vs_Eta", 24, 0, 1);
  TEfficiency *MatchAndCut_vs_Phi =
      new TEfficiency("MatchAndCut_vs_Phi", "MatchPhi_vs_Phi", 20, -TMath::Pi(), TMath::Pi());
  TEfficiency *MatchAndCut_vs_Pt =
      new TEfficiency("MatchAndCut_vs_Pt", "MatchPhi_vs_Pt", 40, 20, 100);

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
  TH2I *N_Digi[4];
  for (const auto st : geom.STATIONS) {
    N_Cluster[st - 1] = new TH2I(Form("N_Cluster_st%d", st), Form("N_Cluster_st%d", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    N_Digi[st - 1] = new TH2I(Form("N_Digi_st%d", st), Form("N_Digi_st%d", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
  }

  TH1D *x_LowBestQ[4];

  TH1D *Res_MuMatched = new TH1D("Res_MuMatched", "Res_MuMatched; cluster.bestTP().xLoc - mu.xLoc; Entries", 101, -10, 10);
  TH1D *Res_SegMatched = new TH1D("Res_SegMatched", "Res_SegMatched; cluster.bestTP().xLoc - seg.xLoc; Entries", 101, -10, 10);
  TH2I *N_MuMatch[4];
  TH1D *Phi_MuMatch[4];
  TEfficiency *Eff_MuMatch[4];
  TH2D *NMuMatch_vs_phi[4];

  TEfficiency *Eff_SegMatch[4];

  TEfficiency *Eff_DigiMatch[4];
  TH1D *Digi_residual[5][4];
  TH1D *N_DigiPerCluster = new TH1D("N_DigiPerCluster", "N_DigiPerCluster", 40, 0, 40); 
  TH1D *DigiSL = new TH1D("DigiSL", "DigiSL", 6, -0.5, 5.5); 

  for (const auto st : geom.STATIONS) {
    N_MuMatch[st - 1] = new TH2I(Form("N_MuMatch_st%d", st), Form("N_MuMatch_st%d; sector ; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    Eff_MuMatch[st-1] = new TEfficiency(Form("Eff_MuMatch_st%d", st), Form("Eff_MuMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    x_LowBestQ[st-1] = new TH1D(Form("x_LowBestQ_st%d", st), Form("x_LowBestQ_st%d; xLoc; Entries",st ), 100, -220, 220);
    NMuMatch_vs_phi[st-1] = new TH2D(Form("NMuMatch_vs_phi_st%d", st), Form("NMuMatch_vs_phi_st%d ; Muon_Phi; N_MuMatch", st), 50, -TMath::Pi(), TMath::Pi(), 6, 1, 7 );
    Phi_MuMatch[st-1] = new TH1D(Form("Phi_MuMatch_st%d", st), Form("Phi_MuMatch_st%d; Muon_phi; Entries", st), 50,-TMath::Pi(), TMath::Pi() );
    Eff_SegMatch[st-1] = new TEfficiency(Form("Eff_SegMatch_st%d", st), Form("Eff_SegMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    Eff_DigiMatch[st-1] = new TEfficiency(Form("Eff_DigiMatch_st%d", st), Form("Eff_DigiMatch_st%d; sector; wheel", st), 14, -0.5, 13.5, 7, -3.5, 3.5);
    //Digi_residual[st-1] = new TH1D( Form("Digi_residual_st%d", st), Form("Digi_residual_st%d; digi.xLoc-seg.wirePos; entries", st), 100, -20, 20 );
    for (const auto wh : geom.WHEELS){
      Digi_residual[wh+2][st-1] = new TH1D( Form("Digi_residual_st%d_wh%d", st, wh), Form("Digi_residual_st%d_wh%d; cluster.bestSeg.xLoc - digi.xLoc; entries", st, wh), 50, -10, 10 );
    }
   }


  double nClustersGhosts{};
  double ooTHQCount{};
  double nClusters{};

  if (fChain == 0) return;

  Long64_t n_entries = fChain->GetEntriesFast();

  // CB unused
  //Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < n_entries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //nb = 
    fChain->GetEntry(jentry);
    //nbytes += nb;

    // ########## LOOP ON EVENTS #############
    if (jentry % 100 == 0) {
      std::cout << "Processing event: " << jentry << '\r' << std::flush;
    }

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
    // ########## CREATE Digis std::vector #############
    std::vector<Digi> digis;
    for (std::size_t i = 0; i < digi_nDigis; ++i ){
      digis.emplace_back(Digi(geom, i, digi_wheel->at(i), digi_sector->at(i), 
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
    auto clusters = buildClusters(geom, tps, segments, digis, X_CUT, DIGI_CUT);

    // ########## ATTEMPT cluster - muon extrapolation matching #############
    for (auto &cluster: clusters){
      for (std::size_t iMu{}; iMu < mu_nMuons; ++iMu){

        for (std::size_t i{}; i < mu_nMatches->at(iMu); ++i){
          int muTrkWheel = getXY<float>(mu_matches_wheel, iMu, i);
          int muTrkStation = getXY<float>(mu_matches_station, iMu, i);
          int muTrkSector = getXY<float>(mu_matches_sector, iMu, i);
          if (muTrkSector == 13 ) muTrkSector = 4;
          if (muTrkSector == 14 ) muTrkSector = 10; 
          double muTrkX = getXY<float>(mu_matches_x, iMu, i);
          double edgeX = getXY<float>(mu_matches_edgeX, iMu, i);
          double edgeY = getXY<float>(mu_matches_edgeY, iMu, i);

          cluster.matchMu(muTrkWheel, muTrkStation, muTrkSector, edgeX, edgeY, muTrkX, i, iMu);
          cluster.matchDigi(digis, DIGI_CUT);
        }
      }
    }


    // ########## RUN SOME ANALYSIS #############
    for (auto const &cluster : clusters) {
      auto wh{cluster.wheel};
      auto sec{cluster.sector};
      if (sec == 13 ) sec = 4;
      if (sec == 14 ) sec = 10;
      auto st{cluster.station};

      // ########## Study TP ghost distribution #############
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

      // ########## Study segment matching #############

      Eff_SegMatch[st-1]->Fill(cluster.segMatched, sec, wh);
      for (const auto s : segments)
      if (cluster.foundTP && s.wheel == wh && s.sector == sec && s.station == st){
        Res_SegMatched->Fill( cluster.bestTP().xLoc - s.xLoc );
      }
      
      // ########## Study digi clusters #############
      Eff_DigiMatch[st-1]->Fill(cluster.digiMatched, sec, wh);
      if (cluster.digiMatched){
        for (auto &digi : cluster.matchedDigi()){
          Digi_residual[wh+2][st-1] ->Fill( cluster.bestSeg().xLoc - digi.xLoc );
        }
      }
     
      N_DigiPerCluster->Fill(cluster.nDigi());
      DigiSL->Fill(cluster.digiSL());

      if (cluster.foundDigi) N_Digi[st - 1]->Fill(sec, wh);

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
          if (tp.index != other_tp.index && tp.match(other_tp, PHI_CUT, T0_CUT)) {    
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
        Match_vs_Eta->Fill((tp.matches.size() > 0), std::abs(gen_eta->at(0)));
        Match_vs_Phi->Fill((tp.matches.size() > 0), std::abs(gen_phi->at(0)));
        Match_vs_Pt->Fill((tp.matches.size() > 0), std::abs(gen_pt->at(0)));
      }
      if (tp.quality == 1 && tp.matches.size() > 0) {
        LowQ_more1HQ_Phi->Fill(tp.t0);
        BX_LowQ_more1HQ->Fill(tp.BX);
      }
    }
  }
  double ghostFraction = nClustersGhosts / nClusters;

  std::cout << " Ratio LQ/selected with phi=  " << LowQ_matched->GetEntries() << "/ "
       << t0_LowQuality->GetEntries() << " = "
       << LowQ_matched->GetEntries() / t0_LowQuality->GetEntries() << std::endl;
  std::cout << " Ratio LQ/matched with phi=  " << LowQ_more1HQ_Phi->GetEntries() << "/ "
       << t0_LowQuality->GetEntries() << " = "
       << LowQ_more1HQ_Phi->GetEntries() / t0_LowQuality->GetEntries() << std::endl;
  std::cout << " Fraction of clusters with ghost (" << nClustersGhosts << ") on total (" << nClusters
       << ") = " << ghostFraction << std::endl;
  std::cout << " HQ out of time clusters: " << ooTHQCount << std::endl;


// ########## DRAW THE ANALYSIS HISTOS  #############

  // TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 500, 500, 500, 500);
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

  // TCanvas *BXCanvas = new TCanvas("BXCanvas", "BXCanvas", 500, 500, 500, 500);
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
  

  TCanvas *DigiMatch = new TCanvas("DigiMatch", "DigiMatch", 500, 500, 500, 500);
  DigiMatch->Divide(2, 2);
  for (int i = 1; i <5; ++i){
    DigiMatch->cd(i);
    Eff_DigiMatch[i-1]->Draw("COLZ");
  }

  TCanvas *DigiRes =  new TCanvas("DigiRes", "DigiRes", 500, 500, 500, 500);
  DigiRes->Divide(5,4);
  int cd = 1;
  for (int i = 1; i <5; ++i){
    for (int j = 0; j< 5; ++j){
      DigiRes->cd(cd);
      Digi_residual[j][i-1]->Draw("COLZ");
      cd+=1;
    }
  }

  TCanvas *DigiCluster = new TCanvas("DigiCluster", "DigiCluster", 500, 500, 500, 500);
  DigiCluster->Divide(1, 2);
  DigiCluster->cd(1);
  N_DigiPerCluster->Draw();
  DigiCluster->cd(2);
  DigiSL->Draw();
  
  TCanvas *DigiCanvas = new TCanvas("DigiCanvas", "DigiCanvas", 500, 500, 500, 500);
  DigiCanvas->Divide(2, 2);
  for (int i = 1; i <5; ++i){
    DigiCanvas->cd(i);
    N_Digi[i-1]->Draw("COLZ");
  }

  outputFile.Write();
  outputFile.Close();
}


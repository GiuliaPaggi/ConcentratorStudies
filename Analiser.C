#define Analiser_cxx
#include "Analiser.h"
#include "Cluster.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TEfficiency.h>
#include <vector>
#include <iostream>
#include "TriggerPrimitive.h"


using namespace std;

void MakeClusters(Cluster clustersArray[5][4][12], TriggerPrimitive TPS[], int sz, double xCut){
    for (int wheel = -2; wheel < 3; ++wheel) {
        for (int station = 1; station < 5; ++station) {
            for (int sector = 1; sector < 13; ++sector)
               clustersArray[wheel+2][station-1][sector-1] = Cluster(TPS, sz, xCut, station, wheel, sector);
        }
    }
};

void Analiser::Loop()
{
   TFile *file =new TFile("DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25.root");     //"VtxSmeared/DTDPGNtuple_12_4_SingleMu_20-100pT_Eta1p25_VtxSmeared.root"
   
   TH1D *t0_AllQuality = new TH1D ("t0_AllQuality", "t0_AllQuality", 100, -9800, -9200);
   TH1D *t0_HighQuality = new TH1D ("t0_HighQuality", "t0_HighQuality", 100, -9800, -9200);       //3+3 4+4 3+4
   TH1D *t0_IntermediateQuality = new TH1D ("t0_IntermediateQuality", "t0_IntermediateQuality", 100, -9800, -9200); // 3+3 4+4 3+4 4+2 3+2   PER ORA NON CI SONO
   TH1D *t0_LowQuality = new TH1D ("t0_LowQuality", "t0_LowQuality", 100 , -9800, -9200);
   TH1D *t0_Selected = new TH1D ("t0_Selected", "t0_Selected", 100, -9800, -9200);
   
   TH1D *LowQ_matched = new TH1D ("LowQ_matched", "LowQ_matched", 100, -9800, -9200);
   TH1D *HighQ_matched = new TH1D("HighQ_matched", "HighQ_matched", 100, -9800, -9200);
   
   TH1D *BX_LowQuality = new TH1D ("BX_LowQuality", "BX_LowQuality", 24 , -392, -368);
   TH1D *BX_LowQ_matched = new TH1D ("BX_LowQ_matched", "BX_LowQ_matched", 24, -392, -368);   
   TH1D *BX_LowQ_more1HQ = new TH1D("BX_LowQ_more1HQ", "BX_LowQ_more1HQ", 24, -392, -368);
   BX_LowQuality->GetXaxis()->SetTitle("BX");
   BX_LowQuality->GetYaxis()->SetTitle("Entries");

   TH1D *t0_Selected_Psi = new TH1D ("t0_Selected_Psi", "t0_Selected_Psi", 100, -9800, -9200);
   TH1D *LowQ_matched_Psi = new TH1D ("LowQ_matched_Psi", "LowQ_matched_Psi", 100, -9800, -9200);
   TH1D *HighQ_matched_Psi = new TH1D ("HighQ_matched_Psi", "HighQ_matched_Psi", 100, -9800, -9200);

   TH1D *LowQ_more1HQ= new TH1D ("LowQ_more1HQ", "LowQ_more1HQ", 100, -9800, -9200);
   TH1D *LowQ_more1HQ_Phi= new TH1D ("LowQ_more1HQ_Phi", "LowQ_more1HQ_Phi", 100, -9800, -9200);

   TH2D *PhiRes_st1_2 = new TH2D ("PhiRes_st1_2", "PhiRes_st1_2", 100, 15, 105, 100, -1, 1);
   PhiRes_st1_2->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st1_2->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_st1_3 = new TH2D ("PhiRes_st1_3", "PhiRes_st1_3", 100, 15, 105, 100, -1, 1);
   PhiRes_st1_3->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st1_3->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_st1_4 = new TH2D ("PhiRes_st1_4", "PhiRes_st1_4", 100, 15, 105, 100, -1, 1);
   PhiRes_st1_4->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st1_4->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_st2_3 = new TH2D ("PhiRes_st2_3", "PhiRes_st2_3", 100, 15, 105, 100, -1, 1);
   PhiRes_st2_3->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st2_3->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_st2_4 = new TH2D ("PhiRes_st2_4", "PhiRes_st2_4", 100, 15, 105, 100, -1, 1);
   PhiRes_st2_4->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st2_4->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_st3_4 = new TH2D ("PhiRes_st3_4", "PhiRes_st3_4", 100, 15, 105, 100, -1, 1);
   PhiRes_st3_4->GetXaxis()->SetTitle("PT (GeV)");
   PhiRes_st3_4->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");   

   TH2D *PhiRes_vs_posizione = new TH2D ("PhiRes_vs_posizione", "PhiRes_vs_posizione", 100, 0, 60, 500, -.1, .1);
   PhiRes_vs_posizione->GetXaxis()->SetTitle("Vertex distance from IP (cm)");
   PhiRes_vs_posizione->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_vs_dxy = new TH2D ("PhiRes_vs_dxy", "PhiRes_vs_dxy", 100, 0, 35, 500, -.1, .1);
   PhiRes_vs_dxy->GetXaxis()->SetTitle("Vertex distance from IP (cm)");
   PhiRes_vs_dxy->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");

   TH2D *PhiRes_vs_dz = new TH2D ("PhiRes_vs_dz", "PhiRes_vs_dz", 100, -45, 45, 500, -.1, .1);
   PhiRes_vs_dz->GetXaxis()->SetTitle("Vertex distance from IP (cm)");
   PhiRes_vs_dz->GetYaxis()->SetTitle("Computed Phi -Phi (rad)");   

   TEfficiency *Match_vs_Eta = new TEfficiency("Match_vs_Eta", "Match_vs_Eta", 24, 0, 1);
   TEfficiency *Match_vs_Phi = new TEfficiency("Match_vs_Phi", "Match_vs_Phi", 20, 0, 3.5);
   TEfficiency *Match_vs_Pt = new TEfficiency("Match_vs_Pt", "Match_vs_Pt", 40, 20, 100);
  /// TEfficiency *Match = new TEfficiency("Match", "Match", 20, -9800, -9200);

   TEfficiency *MatchAndCut_vs_Eta = new TEfficiency("MatchAndCut_vs_Eta", "MatchPhi_vs_Eta", 24, 0, 1);
   TEfficiency *MatchAndCut_vs_Phi = new TEfficiency("MatchAndCut_vs_Phi", "MatchPhi_vs_Phi", 20, 0, 3.5);
   TEfficiency *MatchAndCut_vs_Pt = new TEfficiency("MatchAndCut_vs_Pt", "MatchPhi_vs_Pt", 40, 20, 100);

   // GHOST HISTOS
   /*TH1D OoTGhosts_style("OoTGhosts","OoTGhosts_style", 10, 0, 10);
   TH1D *OoTGhosts[5][4][12];

   TH1D ITGhosts_style("ITGhosts", "ITGhosts_style", 10, 0, 10);
   TH1D *ITGhosts[5][4][12];*/
   TH1D *OoTGhosts = new TH1D("OoTGhosts", "OoTGhosts", 20, 0, 20);
   TH1D *ITGhosts = new TH1D("ITGhosts", "ITGhosts", 20, 0, 20);

   TH1I *BX_ITGhosts = new TH1I("BX_ITGhosts", "BX_ITGhosts", 24 , -392, -368);
   TH1I *BX_OoTGhosts = new TH1I("BX_OoTGhosts", "BX_OoTGhosts", 24 , -392, -368);
   TH1D *Res_ITGhosts = new TH1D("Res_ITGhosts", "Res_ITGhosts", 110, -5.5, 5.5);
   TH1D *Res_OoTGhosts = new TH1D("Res_OoTGhosts", "Res_OoTGhosts", 110, -5.5, 5.5);

   TH2D *Q_OoTGhosts = new TH2D("Q_OoTGhosts", "Q_OoTGhosts", 10, 0, 10, 10, 0, 10);
   Q_OoTGhosts->GetXaxis()->SetTitle("High Quality");
   Q_OoTGhosts->GetYaxis()->SetTitle("Out of time Quality");
   TH2D *Q_ITGhosts = new TH2D("Q_ITGhosts", "Q_ITGhosts", 10, 0, 10, 10, 0, 10);
   Q_ITGhosts->GetXaxis()->SetTitle("High Quality");
   Q_ITGhosts->GetYaxis()->SetTitle("In time Quality");

   /*for(int st=1; st < 5; st++) {
      for (int wheel = -2 ; wheel < 3; ++wheel){
         for (int sector = 1 ; sector < 13; ++sector) {
            OoTGhosts[wheel+2][st-1][sector-1]= new TH1D(OoTGhosts_style);
            OoTGhosts[wheel+2][st-1][sector-1]->SetName(Form("OoTGhosts_%d_%d_%d", wheel, st, sector));

            ITGhosts[wheel+2][st-1][sector-1]= new TH1D(ITGhosts_style);
            ITGhosts[wheel+2][st-1][sector-1]->SetName(Form("ITGhosts_%d_%d_%d", wheel, st, sector));
         }
      }
   }*/

   TH1I *N_Ghost = new TH1I("N_Ghost", "N_Ghost", 20, 0, 20);
   TH1I *Q_Best = new TH1I("Q_Best", "Q_Best", 10, 0, 10);
   TH1I *Q_Ghost = new TH1I("Q_Ghost", "Q_Ghost", 10, 0, 10);

   TEfficiency Ghost_eff_style("Ghost_eff_st", "Ghost_eff_st", 14, 0, 13, 7, -3.5, 3.5);
   TEfficiency *Ghost_eff[4];

   TH2I N_Cluster_style("N_Cluster_st","N_Cluster_st", 14, 0, 13, 7, -3.5, 3.5);    //nel secondo metto titolo;titolo assex;titolo asse y
   TH2I *N_Cluster[4];

   for (int st = 1; st < 5; ++st) {
      N_Cluster[st-1] = new TH2I(N_Cluster_style);
      N_Cluster[st-1]->SetName(Form("N_Cluster_st%d", st));
      Ghost_eff[st-1] = new TEfficiency(Form("Ghost_eff_st%d", st), Form("Ghost_eff_st%d", st), 14, 0, 13, 7, -3.5, 3.5);
      //Ghost_eff[st-1]->SetName(Form("Ghost_eff_st%d", st));

   }


//   In a ROOT session, you can do:
//      root> .L Analiser.C
//      root> Analiser t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

   const double MB[4] = {402.2, 490.5, 597.5, 700.0};
   const int LowQualityCut = 0;
   const int HighQualityCut = 5;
   const double PsiCut = TMath::Pi()/6;
   double PhiCut = 0.02;
   const double PhiCut2 = 0.01;
   const double TimeCut = 12.5;
   const double RealBX = 380; 
   const double xcut = 5.0; 
   double ClusterCount = 0;
   double OoTHQCount = 0;
   double AllEvents = 0;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      

      // ########## LOOP ON EVENTS #############
      if (jentry%1000==0) cout<<"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b Processing event "<<jentry<<flush;

      // if (Cut(ientry) < 0) continue;

            if ( abs(gen_pdgId->at(0)) != 13 || abs(gen_eta->at(0)) > 0.8  ) continue;
             //cout << "---------------------------------------" << endl;
            vector <int> sectors;
            TriggerPrimitive TPs[ph2TpgPhiEmuAm_nTrigs];

            // ########## LOOP ON TRIGGERS PRIMITIVES #############
            for (int j = 0; j < ph2TpgPhiEmuAm_nTrigs; ++j){ 
                   
               TPs[j] = TriggerPrimitive(j, ph2TpgPhiEmuAm_wheel->at(j), ph2TpgPhiEmuAm_sector->at(j), ph2TpgPhiEmuAm_station->at(j),  ph2TpgPhiEmuAm_quality ->at(j),  ph2TpgPhiEmuAm_phi ->at(j),  ph2TpgPhiEmuAm_phiB->at(j),  ph2TpgPhiEmuAm_BX->at(j), ph2TpgPhiEmuAm_t0->at(j), ph2TpgPhiEmuAm_posLoc_x->at(j));
               if (TPs[j].quality == 1) {
                  t0_LowQuality->Fill(TPs[j].t0);
                  BX_LowQuality->Fill(TPs[j].BX);
               }
            }
            
            //Find clusters
            Cluster Clusters[5][4][12];
            MakeClusters(Clusters, TPs, ph2TpgPhiEmuAm_nTrigs, xcut);

            for (int wheel = -2; wheel < 3; ++wheel) {    
               for (int st = 1; st <5; ++st){
                  for (int sec = 1; sec < 13; ++sec){
                     AllEvents +=1;
                     Ghost_eff[st-1]->Fill(!Clusters[wheel+2][st-1][sec-1].Isolated, sec, wheel);
                     if (Clusters[wheel+2][st-1][sec-1].OoTHQ == true ) OoTHQCount+=1;
                     if (Clusters[wheel+2][st-1][sec-1].Isolated == true ) continue;
                     ClusterCount +=1;
                     int OoTSize = Clusters[wheel+2][st-1][sec-1].GetOoTSize();
                     int ITSize = Clusters[wheel+2][st-1][sec-1].GetITSize();

                     OoTGhosts->Fill(OoTSize);
                     ITGhosts->Fill(ITSize);
                     N_Ghost->Fill(OoTSize+ITSize);

                     vector<double> OoTBxs = Clusters[wheel+2][st-1][sec-1].GetOoTBX();
                     vector<double> ITBxs = Clusters[wheel+2][st-1][sec-1].GetITBX();

                     vector<double> OoTRes = Clusters[wheel+2][st-1][sec-1].GetOoTResidual();
                     vector<double> ITRes = Clusters[wheel+2][st-1][sec-1].GetITResidual();

                     vector<int> OoTQual = Clusters[wheel+2][st-1][sec-1].GetOoTQualities();
                     vector<int> ITQual = Clusters[wheel+2][st-1][sec-1].GetITQualities();
                     int HQ = Clusters[wheel+2][st-1][sec-1].GetBestQuality();
                     cout << Clusters[wheel+2][st-1][sec-1]._BestQuality.BX << endl;
                     Q_Best->Fill(HQ);


                     for (int index = 0; index < OoTSize; ++index){
                        BX_OoTGhosts->Fill( OoTBxs[index] );
                        Res_OoTGhosts->Fill( OoTRes[index] );
                        Q_OoTGhosts->Fill( HQ, OoTQual[index] );
                        Q_Ghost->Fill(OoTQual[index]);
                     }
                     for (int index = 0; index < ITSize; ++index) {
                        BX_ITGhosts->Fill( ITBxs[index] );
                        Res_ITGhosts->Fill( ITRes[index] );
                        Q_ITGhosts->Fill( HQ, ITQual[index] );
                        Q_Ghost->Fill(ITQual[index]);
                     }

                     N_Cluster[st-1]->Fill(sec, wheel);
                  }
               }
            }

            for ( TriggerPrimitive &tp : TPs) {
               if (tp.quality > 5) {
                  if (tp.hasMatched) continue;
                  t0_Selected->Fill(tp.t0);
                  for ( TriggerPrimitive &trigprim : TPs) {
                     if (trigprim.index == tp.index)  continue;
                     if( tp.Match(trigprim, PhiCut, TimeCut) ){
                        if(trigprim.quality == 1) {
                                 t0_Selected->Fill(trigprim.t0); 
                                 LowQ_matched->Fill(trigprim.t0);
                                 BX_LowQ_matched->Fill(trigprim.BX);
                        }
                        else {
                           t0_Selected->Fill(trigprim.t0); 
                        }
                     } 
                  }
               }
            }

            for ( TriggerPrimitive tp : TPs) {
               //cout << tp.hasMatched  << "      ||       " << tp.Matches.size() << endl;
               if (tp.quality == 1) {
                  MatchAndCut_vs_Eta->Fill(tp.hasMatched, abs(gen_eta->at(0)));
                  MatchAndCut_vs_Phi->Fill(tp.hasMatched, abs(gen_phi->at(0)));
                  MatchAndCut_vs_Pt->Fill(tp.hasMatched, abs(gen_pt->at(0)));
                  Match_vs_Eta->Fill( (tp.Matches.size() > 0) , abs(gen_eta->at(0)));
                  Match_vs_Phi->Fill( (tp.Matches.size() > 0) , abs(gen_phi->at(0)));
                  Match_vs_Pt->Fill( (tp.Matches.size() > 0) , abs(gen_pt->at(0)) );
               }
               if (tp.quality == 1 && tp.Matches.size() > 0 ) {
                  LowQ_more1HQ_Phi->Fill(tp.t0);
                  BX_LowQ_more1HQ->Fill(tp.BX);
                  
               }
            } 
   }
   double Ghostfraction = ClusterCount/AllEvents;

   cout << " Ratio LQ/selected with phi=  " << LowQ_matched->GetEntries() << "/ " << t0_LowQuality->GetEntries()  << " = "<<  LowQ_matched->GetEntries()/t0_LowQuality->GetEntries() << endl;
   cout << " Ratio LQ/matched with phi=  " << LowQ_more1HQ_Phi->GetEntries() << "/ " << t0_LowQuality->GetEntries()  << " = "<<  LowQ_more1HQ_Phi->GetEntries()/t0_LowQuality->GetEntries() << endl;
   cout << " Fraction of events with ghost (" << ClusterCount << ") on total (" << AllEvents << ") = " << Ghostfraction << endl;
   cout << " HQ out of time events: " << OoTHQCount << endl;

   TCanvas *canvas2 = new TCanvas ("canvas2", "canvas2", 500, 500, 500, 500);
   gPad->SetLogy();
   t0_LowQuality->SetLineColor(kBlue);
   t0_LowQuality->GetXaxis()->SetTitle(" t0 (ns)");
   t0_LowQuality->GetYaxis()->SetTitle(" Entries");   
   t0_LowQuality->Draw();
   //LowQ_matched_Psi->SetLineColor(kRed);
   //LowQ_matched_Psi->Draw("same");
   LowQ_matched->SetLineColor(kRed);
   LowQ_matched->Draw("same");
   //LowQ_more1HQ->Draw("same");
   LowQ_more1HQ_Phi->SetLineColor(kGreen+2);
   LowQ_more1HQ_Phi->Draw("same");
   


   auto leg = new TLegend (0.1, 0.7, 0.4, 0.9);
   leg->AddEntry(t0_LowQuality, "Q1");
   leg->AddEntry(LowQ_matched, "Q1 - #phi cut and time cut");
   //leg->AddEntry(LowQ_matched_Psi, "Q1 - with psi&time cut ");
   //leg->AddEntry(LowQ_more1HQ, "Q1 - associated HQ psi cut, no time cut ");
   leg->AddEntry(LowQ_more1HQ_Phi, "Q1 - #phi cut, no time cut");
   leg->Draw("same");


   TCanvas *canvas = new TCanvas ("canvas", "canvas", 500, 500, 500, 500);
   canvas->Divide(2, 2);
   
   
   canvas->cd(1);
   gPad->SetLogy();

   t0_AllQuality->Draw();
   //t0_LowQuality->SetLineColor(kRed);
   t0_LowQuality->Draw("same");
   t0_HighQuality->SetLineColor(kBlack);
   t0_HighQuality->Draw("same");

   auto leg1 = new TLegend (0.1, 0.7, 0.35, 0.9);
   leg1->AddEntry(t0_LowQuality, "LQ (1-3) primitives");
   leg1->AddEntry(t0_HighQuality, "HQ (6-7-8) primitives");
   leg1->AddEntry(t0_AllQuality, "All primitives in #eta < 0.8");
   leg1->Draw("same");

   canvas->cd(2);
   gPad->SetLogy();
   //t0_LowQuality->SetLineColor(kBlack);
   t0_LowQuality->Draw();
   LowQ_matched_Psi->SetLineColor(kRed);
   LowQ_matched_Psi->Draw("same");
   //LowQ_matched->SetLineColor(kGreen);
   //LowQ_matched->Draw("same");
   LowQ_more1HQ->Draw("same");
   //LowQ_more1HQ_Phi->SetLineColor(kRed);
   LowQ_more1HQ_Phi->Draw("same");


   auto leg2 = new TLegend (0.1, 0.7, 0.4, 0.9);
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

   auto leg3 = new TLegend (0.1, 0.7, 0.35, 0.9);
   leg3->AddEntry(LowQ_matched_Psi, "Q1 primitives associated");
   leg3->AddEntry(HighQ_matched_Psi, "HQ primitives associated");
   leg3->AddEntry(t0_Selected_Psi, "All primitives");
   leg3->Draw("same");

   canvas->cd(4);
   gPad->SetLogy();
   
   t0_Selected->Draw();
   HighQ_matched->SetLineColor(kRed);
   HighQ_matched->Draw("same");
   //LowQ_matched->SetLineColor(kGreen);
   LowQ_matched->Draw("same");

   auto leg4 = new TLegend (0.1, 0.7, 0.35, 0.9);
   leg4->AddEntry(LowQ_matched, "Q1 primitives associated");
   leg4->AddEntry(HighQ_matched, "HQ primitives associated");
   leg4->AddEntry(t0_Selected, "All primitives");
   leg4->Draw("same");
   
   TCanvas *BXCanvas = new TCanvas("BXCanvas", "BXCanvas", 500, 500, 500, 500);
   BX_LowQuality->Draw();
   BX_LowQ_matched->SetLineColor(kRed);    
   BX_LowQ_matched->Draw("same");
   BX_LowQ_more1HQ->SetLineColor(kGreen+2);
   BX_LowQ_more1HQ->Draw("same");

   auto *legend = new TLegend (0.1, 0.7, 0.35, 0.9);
   legend->AddEntry(BX_LowQuality, "Q1 primitive");
   legend->AddEntry(BX_LowQ_matched, "Q1 selected #phi match and time cut");
   legend->AddEntry(BX_LowQ_more1HQ, "Q1 selected #phi match " );
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

   TCanvas *ClusterProvaCanvas = new TCanvas("ClusterProvaCanvas", "ClusterProvaCanvas", 500, 500, 500, 500);
   ClusterProvaCanvas->Divide(2, 2);
   ClusterProvaCanvas->cd(1);
   OoTGhosts->Draw();

   ClusterProvaCanvas->cd(2);
   ITGhosts->Draw();

   ClusterProvaCanvas->cd(3);
   N_Ghost->Draw();

   ClusterProvaCanvas->cd(4);
   Q_Ghost->SetLineColor(kOrange+7);
   Q_Ghost->Draw();
   Q_Best->Draw("same");
   auto *q_legend = new TLegend (0.1, 0.7, 0.35, 0.9);
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
   Q_OoTGhosts->Draw();

   ClusterBXCanvas->cd(4);
   Q_ITGhosts->Draw();

   TCanvas *StatCanvas = new TCanvas("StatCanvas", "StatCanvas", 500, 500, 500, 500);
   StatCanvas->Divide(2, 2);

   for (int i = 1; i<5; ++i){
      StatCanvas->cd(i);
      N_Cluster[i-1]->Draw("COLZ");
   }

   TCanvas *GhostEffCanvas = new TCanvas("GhostEffCanvas", "GhostEffCanvas", 500, 500, 500, 500);
   GhostEffCanvas->Divide(2, 2);
   for (int i = 1; i<5; ++i){
      GhostEffCanvas->cd(i);
      Ghost_eff[i-1]->Draw("COLZ");
   }

   //scanvas_residui->SaveAs("residui.png");

}


/*   canvas->cd(2);
   gPad->SetLogy();
   
   t0_Selected->Draw();
   HighQ_matched->SetLineColor(kGreen);
   HighQ_matched->Draw("same");
   LowQ_matched->SetLineColor(kRed);
   LowQ_matched->Draw("same");

   auto leg2 = new TLegend (0.1, 0.7, 0.35, 0.9);
   leg2->AddEntry(LowQ_matched, "LQ primitives associated (time cut on Q1) ");
   leg2->AddEntry(HighQ_matched, "HQ primitives associated");
   leg2->AddEntry(t0_Selected, "All primitives selected in #eta < 0.8");
   leg2->Draw("same");*/


               /*
               
               AM_Phi->Fill(ph2TpgPhiEmuAm_phi->at(j));
               Phi_correlation->Fill(gen_phi->at(i), ph2TpgPhiEmuAm_phi->at(j));
               Phi_vs_Pt->Fill(gen_pt->at(i), ph2TpgPhiEmuAm_phi->at(j));
               PhiB_vs_Pt->Fill(gen_pt->at(i), ph2TpgPhiEmuAm_phiB->at(j));
               PhiB_vs_Phi->Fill(gen_phi->at(i), ph2TpgPhiEmuAm_phiB->at(j));
               Am_wheel_vs_digi_wheel->Fill(digi_wheel->at(i), ph2TpgPhiEmuAm_wheel->at(j));
               AM_station_vs_digi_station->Fill(digi_station->at(i), ph2TpgPhiEmuAm_station->at(j));

               */
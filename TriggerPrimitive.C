#include "TriggerPrimitive.h"



TriggerPrimitive::TriggerPrimitive( short i, short ph2TpgPhiEmuAm_wheel, short ph2TpgPhiEmuAm_sector, short ph2TpgPhiEmuAm_station, short ph2TpgPhiEmuAm_quality, 
       int ph2TpgPhiEmuAm_phi, int ph2TpgPhiEmuAm_phiB, int ph2TpgPhiEmuAm_BX, int ph2TpgPhiEmuAm_t0, float ph2TpgPhiEmuAm_posLoc_x){
   
   index = i;
   wheel = ph2TpgPhiEmuAm_wheel;
   sector = ph2TpgPhiEmuAm_sector;
   station = ph2TpgPhiEmuAm_station;
   quality = ph2TpgPhiEmuAm_quality;
   phi = TMath::Pi()/6*ph2TpgPhiEmuAm_sector + ph2TpgPhiEmuAm_phi *.5 / 65536;    //(17 bits between 0.5 and 0.5 rad)
   if (phi > TMath::Pi()*2) phi = phi - TMath::Pi()*2 ;
   phiB = ph2TpgPhiEmuAm_phiB* 2. /4096;
   psi = phi + phiB;
   BX = ph2TpgPhiEmuAm_BX;
   t0 = ph2TpgPhiEmuAm_t0;
   xLoc = ph2TpgPhiEmuAm_posLoc_x;
};

void TriggerPrimitive::ComputeExpectedPhi() {
   const double MB[4] = {402.2, 490.5, 597.5, 700.0};
   
   for (int stat = 1; stat < 5; ++stat) {
      if (station == stat) phiExpected[stat] = phi;
      else {
         double ExpPhi = psi - TMath::ASin( TMath::Sin(phiB) * MB[station -1]/MB[stat-1] );
         phiExpected[stat] = ExpPhi;
      }
   }
   computedPhi = true;
};

bool TriggerPrimitive::Match(TriggerPrimitive &TP, double phicut, double timecut){
   if (TP.index == index) return false; 

   if (!computedPhi)  ComputeExpectedPhi();

   double DeltaPhiExp = 0;
   abs(phiExpected[TP.station] - TP.phi) < TMath::Pi() ? DeltaPhiExp = abs(phiExpected[TP.station] - TP.phi) : DeltaPhiExp = abs( 2*TMath::Pi() - abs(phiExpected[TP.station] - TP.phi));

   double Deltat0 = abs(t0 - TP.t0);

   if (DeltaPhiExp < phicut) {
      TP.Matches.push_back(index);        // NELLE QUALITà BASSE METTO L'INDICE DI QUELLA ALTA
      Matches.push_back(TP.index);        // NELLE QUALITà ALTE METTO L'INDICE DI QUELLE CHE MATCHANO
      
      if (wheel == 0 && abs(TP.wheel) < 3 ) { 

         if (TP.quality == 1 && Deltat0 < timecut) {
            TP.hasMatched = true;
            hasMatched = true;
            return true;

         } 
         else if (TP.quality > 1){
            TP.hasMatched = true;
            hasMatched = true;
            return true;
         } 
      }

      if (wheel > 0 && TP.wheel >= wheel) { 
         if (TP.quality == 1 && Deltat0 < timecut) {
            TP.hasMatched = true;
            hasMatched = true;
            return true;

            }  
         else if (TP.quality > 1){
            TP.hasMatched = true;
            hasMatched = true;
            return true;
            } 
      }

      if (wheel < 0 && TP.wheel <= wheel && TP.wheel != -5 && wheel !=5 ) { 

         if (TP.quality == 1 && Deltat0 < timecut) {
            TP.hasMatched = true;
            hasMatched = true;
            return true;
            }  
         else if (TP.quality > 1){
            TP.hasMatched = true;
            hasMatched = true;
            return true;
            } 
      }
   }
   return false;
};

//Makes cluster from the TP that calls it, checking all the TP in list, returns vector of index of the cluster

vector<int> TriggerPrimitive::MakeCluster(TriggerPrimitive listOfPrimitives[], int size,  double xcut) {
   
   vector<int> Cluster;

   for (int i = 0; i < size; ++i){
      TriggerPrimitive element = listOfPrimitives[i];   
      //double DeltaPhi = 0;
      //abs(element.phi-phi) < TMath::Pi() ? DeltaPhi = abs(element.phi-phi) : DeltaPhi = abs( 2*TMath::Pi() - abs(element.phi-phi)); 
      double DeltaxLoc = abs(element.xLoc - xLoc);

      if (element.station == station && element.wheel == wheel && DeltaxLoc < xcut ) {
         if (element.index != index) {
         inCluster = true;
         element.inCluster = true;
         Cluster.push_back(element.index);
         }
      }
   }

   return Cluster;
};

void TriggerPrimitive::FindHigherQuality (TriggerPrimitive listOfPrimitives[], vector<int> clusterIndex){
   if (!inCluster) return;
   if (clusterIndex.size() == 1 ) return;

   int maxQ = quality;
   int maxQindex = index; 
   for ( auto element : clusterIndex) {
      if (hasHigherQuality) return;
      if (listOfPrimitives[element].quality > maxQ ) { 
         maxQ = listOfPrimitives[element].quality;
         maxQindex = element;
      }
   }

   listOfPrimitives[maxQindex].hasHigherQuality = true;
};

//vector<int> TriggerPrimitive::SelectRightBX (TriggerPrimitive listOfPrimitives[], vector<int> clusterIndex){
//   vector <int> RightBXCluster;        // = {index}
//   for (auto element : clusterIndex) {
//      if (380 == listOfPrimitives[element].BX) RightBXCluster.push_back(element);
//   }
//   return RightBXCluster;
//};

void TriggerPrimitive::CheckBX() {
   if ( abs(BX+380) < 1) hasRightBX = true;
   else isGhostOutOfTime = true;
};
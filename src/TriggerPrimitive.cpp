#include "include/TriggerPrimitive.h"

#include "TMath.h"

TriggerPrimitive::TriggerPrimitive(std::size_t i, int tpg_wheel, int tpg_sector, int tpg_station,
                                   int tpg_quality, int tpg_phi, int tpg_phiB, int tpg_BX,
                                   int tpg_t0, float tpg_posLoc_x)
    : index{i},
      wheel{tpg_wheel},
      sector{tpg_sector},
      station{tpg_station},
      quality{tpg_quality},
      BX{tpg_BX},
      t0{tpg_t0},
      xLoc{tpg_posLoc_x} {
  //(17 bits between 0.5 and 0.5 rad)
  phi = TMath::Pi() / 6 * sector + tpg_phi * .5 / 65536; // CB what if phi < 0
  if (phi > TMath::Pi() * 2) {
    phi = phi - TMath::Pi() * 2;
  }

  phiB = tpg_phiB * 2. / 4096;

  psi = phi + phiB;
};

void TriggerPrimitive::computeExpectedPhi() {
  const double MB[4] = {402.2, 490.5, 597.5, 700.0}; // CB this could go in GEOM
  for (int stat = 1; stat < 5; ++stat) { // CB could use GEOM
    if (station == stat){
      phiExpected[stat-1] = phi;
    }
    else {
      double ExpPhi = psi - TMath::ASin(TMath::Sin(phiB) * MB[station - 1] / MB[stat - 1]);
      phiExpected[stat-1] = ExpPhi; 
      }
  }
  computedPhi = true;
};

bool TriggerPrimitive::match(TriggerPrimitive &tp, double phiCut, double timeCut) {
  if (tp.index == index) return false;
  if (!computedPhi) computeExpectedPhi();

  double deltaPhiExp = std::abs(phiExpected[tp.station] - tp.phi) < TMath::Pi()
                       ? std::abs(phiExpected[tp.station] - tp.phi)
                       : std::abs(2 * TMath::Pi() - std::abs(phiExpected[tp.station] - tp.phi));

  double deltat0 = std::abs(t0 - tp.t0);

  if (deltaPhiExp < phiCut) {

    tp.matches.push_back(index);  // NELLE QUALITà BASSE METTO L'INDICE DI QUELLA ALTA
    matches.push_back(tp.index);  // NELLE QUALITà ALTE METTO L'INDICE DI QUELLE CHE MATCHANO

    if (wheel == 0 && std::abs(tp.wheel) < 3) {
      if (tp.quality == 1 && deltat0 < timeCut) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;

      } else if (tp.quality > 1) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;
      }
    }

    if (wheel > 0 && tp.wheel >= wheel) {
      if (tp.quality == 1 && deltat0 < timeCut) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;

      } else if (tp.quality > 1) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;
      }
    }

    if (wheel < 0 && tp.wheel <= wheel && tp.wheel != -5 && wheel != 5) { // CB wheel?
      if (tp.quality == 1 && deltat0 < timeCut) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;
      } else if (tp.quality > 1) {
        tp.hasMatched = true;
        hasMatched = true;
        return true;
      }
    }
  }
  return false;
};

/*/ Makes cluster from the TP that calls it, checking all the TP in list, returns
// vector of index of the cluster

vector<int> TriggerPrimitive::makeCluster(TriggerPrimitive listOfPrimitives[], int size,
                                          double xcut) {
  vector<int> Cluster;

  for (int i = 0; i < size; ++i) {
    TriggerPrimitive element = listOfPrimitives[i];
    double DeltaxLoc = std::abs(element.xLoc - xLoc);

    if (element.station == station && element.wheel == wheel && DeltaxLoc < xcut) {
      if (element.index != index) {
        inCluster = true;
        element.inCluster = true;
        Cluster.push_back(element.index);
      }
    }
  }

  return Cluster;
};*/

void TriggerPrimitive::findHigherQuality(TriggerPrimitive listOfPrimitives[],
                                        const std::vector<int>& clusterIndex) {
  if (!inCluster) return;
  if (clusterIndex.size() == 1) return;

  int maxQ{quality};
  std::size_t maxQIndex{index};

  for (auto element : clusterIndex) {
    if (hasHighestQuality) return;
    if (listOfPrimitives[element].quality > maxQ) {
      maxQ = listOfPrimitives[element].quality;
      maxQIndex = element;
    }
  }

  listOfPrimitives[maxQIndex].hasHighestQuality = true;
};

void TriggerPrimitive::checkBX() {
  if (std::abs(BX + 380) < 1){ // CB should put somewhere else?
    hasRightBX = true;
  }
  else {
    isGhostOutOfTime = true;
  }
};
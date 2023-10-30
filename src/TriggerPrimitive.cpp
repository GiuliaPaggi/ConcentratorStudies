#include "include/TriggerPrimitive.h"

#include "TMath.h"

constexpr int RIGHT_BX{20};

TriggerPrimitive::TriggerPrimitive(std::size_t i, int tpg_wheel, int tpg_sector, int tpg_station, int tpg_quality,
                                   int tpg_phi, int tpg_phiB, int tpg_BX, int tpg_t0, float tpg_posLoc_x)
    : index{i},
      wheel{tpg_wheel},
      sector{tpg_sector},
      station{tpg_station},
      quality{tpg_quality},
      BX{tpg_BX},
      t0{tpg_t0},
      xLoc{tpg_posLoc_x} {
  //(17 bits between 0.5 and 0.5 rad)
  
  if (sector == 13)
    sector = 4;
  else if (sector == 14)
    sector = 10;

  phi = TMath::Pi() / 6 * sector + tpg_phi * .5 / 65536;  // CB what if phi < 0
  if (phi > TMath::Pi() * 2) {
    phi = phi - TMath::Pi() * 2;
  }

  phiB = tpg_phiB * 2. / 4096;

  psi = phi + phiB;
};

void TriggerPrimitive::computeExpectedPhi() {
  const double MB[4] = {402.2, 490.5, 597.5, 700.0};  // CB this could go in GEOM
  for (int stat = 1; stat < 5; ++stat) {              // CB could use GEOM
    if (station == stat) {
      phiExpected[stat - 1] = phi;
    } else {
      double ExpPhi = psi - TMath::ASin(TMath::Sin(phiB) * MB[station - 1] / MB[stat - 1]);
      phiExpected[stat - 1] = ExpPhi;
    }
  }
  computedPhi = true;
};

bool TriggerPrimitive::Match(TriggerPrimitive &TP, double phicut, double timecut, const std::vector<int> &quals) {
  if (std::find(quals.begin(), quals.end(), quality) == quals.end()) return false;

  if (TP.index == index) return false;
  if (!computedPhi) computeExpectedPhi();

  double DeltaPhiExp = 0;
  std::abs(phiExpected[TP.station] - TP.phi) < TMath::Pi()
      ? DeltaPhiExp = std::abs(phiExpected[TP.station] - TP.phi)
      : DeltaPhiExp = std::abs(2 * TMath::Pi() - std::abs(phiExpected[TP.station] - TP.phi));

  double Deltat0 = std::abs(t0 - TP.t0);

  if (DeltaPhiExp < phicut) {
    TP.matches.push_back(index);  // NELLE QUALITà BASSE METTO L'INDICE DI QUELLA ALTA
    matches.push_back(TP.index);  // NELLE QUALITà ALTE METTO L'INDICE DI QUELLE CHE MATCHANO

    if (wheel == 0 && std::abs(TP.wheel) < 2) {
      if (TP.quality == 1 && Deltat0 < timecut) {
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
    }

    if (wheel < 0 && TP.wheel <= wheel && TP.wheel != -5 && wheel != 5) {
      if (TP.quality == 1 && Deltat0 < timecut) {
        TP.hasMatched = true;
        hasMatched = true;
        return true;
      }
    }
  }
  return false;
};

void TriggerPrimitive::findHigherQuality(TriggerPrimitive listOfPrimitives[], const std::vector<int> &clusterIndex) {
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
  if (std::abs(BX - RIGHT_BX) < 1) {  // CB should put somewhere else?
    hasRightBX = true;
  } else {
    isGhostOutOfTime = true;
  }
};
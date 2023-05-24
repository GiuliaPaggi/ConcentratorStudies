#include "TriggerPrimitive.h"

#include "TMath.h"

TriggerPrimitive::TriggerPrimitive(std::size_t i, int tpg_wheel, int tpg_sector, int tpg_station,
                                   int tpg_quality, int tpg_phi, int tpg_phiB, int tpg_BX,
                                   int tpg_t0, float tpg_posLoc_x)
    : index(i),
      wheel(tpg_wheel),
      sector(tpg_sector),
      station(tpg_station),
      quality(tpg_quality),
      BX(tpg_BX),
      t0(tpg_t0),
      xLoc(tpg_posLoc_x) {
  //(17 bits between 0.5 and 0.5 rad)
  if (tpg_sector == 13) sector = 4;
  else if (tpg_sector == 14) sector = 10;

  phi = TMath::Pi() / 6 * sector + tpg_phi * .5 / 65536;
  if (phi > TMath::Pi() * 2) {
    phi = phi - TMath::Pi() * 2;
  }

  phiB = tpg_phiB * 2. / 4096;

  psi = tpg_phi * .5 / 65536 + phiB;
};

void TriggerPrimitive::ComputeExpectedPhi() {
  const double MB[4] = {402.2, 490.5, 597.5, 700.0};
  for (int stat = 1; stat < 5; ++stat) {
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

bool TriggerPrimitive::Match(TriggerPrimitive &TP, double phicut, double timecut) {
  if (TP.index == index) return false;
  if (!computedPhi) ComputeExpectedPhi();

  double DeltaPhiExp = 0;
  std::abs(phiExpected[TP.station] - TP.phi) < TMath::Pi()
      ? DeltaPhiExp = std::abs(phiExpected[TP.station] - TP.phi)
      : DeltaPhiExp = std::abs(2 * TMath::Pi() - std::abs(phiExpected[TP.station] - TP.phi));

  double Deltat0 = std::abs(t0 - TP.t0);

  if (DeltaPhiExp < phicut) {

    TP.Matches.push_back(index);  // NELLE QUALITà BASSE METTO L'INDICE DI QUELLA ALTA
    Matches.push_back(TP.index);  // NELLE QUALITà ALTE METTO L'INDICE DI QUELLE CHE MATCHANO

    if (wheel == 0 && std::abs(TP.wheel) < 2) {
      if (TP.quality == 1 && Deltat0 < timecut) {
        TP.hasMatched = true;
        hasMatched = true;
        return true;

      } else if (TP.quality > 1) {
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

      } else if (TP.quality > 1) {
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
      } else if (TP.quality > 1) {
        TP.hasMatched = true;
        hasMatched = true;
        return true;
      }
    }
  }
  return false;
};

void TriggerPrimitive::FindHigherQuality(TriggerPrimitive listOfPrimitives[],
                                         vector<int> clusterIndex) {
  if (!inCluster) return;
  if (clusterIndex.size() == 1) return;

  int maxQ = quality;
  int maxQindex = index;
  for (auto element : clusterIndex) {
    if (hasHighestQuality) return;
    if (listOfPrimitives[element].quality > maxQ) {
      maxQ = listOfPrimitives[element].quality;
      maxQindex = element;
    }
  }

  listOfPrimitives[maxQindex].hasHighestQuality = true;
};


void TriggerPrimitive::CheckBX() {
  if (std::abs(BX + 380) < 1)
    hasRIGHT_BX = true;
  else
    isGhostOutOfTime = true;
};
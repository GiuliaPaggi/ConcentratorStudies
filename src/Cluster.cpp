#include "include/Cluster.h"

#include <algorithm>
#include <iostream>

constexpr int RIGHT_BX{20};

Cluster::Cluster(std::vector<TriggerPrimitive> &tps, std::vector<Segment> &segs, std::vector<Digi> &digis, double xCut,
                 double digiCut, int wh, int sec, int st)
    : wheel{wh}, sector{sec}, station{st} {
  // in a given wheel-station-sector finds the in-time highest-quality TP
  // and looks for other TP in a |xCut| interval it's cluster even if no TPs but
  // here's a segment, or more than 10 digis in a SL

  if (sector == 13)
    sector = 4;
  else if (sector == 14)
    sector = 10;

  double clusterX{NAN};

  auto inChamber = [=](auto const &obj) {
    return obj.wheel == wh && obj.station == st && obj.sector == sec && obj.inCluster == false;
  };

  auto inRange = [&clusterX, &xCut](auto const &obj) { return std::abs(obj.xLoc - clusterX) < xCut; };

  auto inRangeDigi = [&clusterX, &digiCut](auto const &obj) { return std::abs(obj.xLoc - clusterX) < digiCut; };

  auto compareTPs = [](auto const &tp1, auto const &tp2) {  // needs to return true if 1 < 2
    return (tp1.BX != RIGHT_BX && tp2.BX == RIGHT_BX) ||
           (tp1.BX == RIGHT_BX && tp2.BX == RIGHT_BX && tp1.quality < tp2.quality);
  };

  auto compareSegs = [](Segment &s1, Segment &s2) { return s1.nPhiHits < s2.nPhiHits; };

  auto toggleCluster = [](auto &collection, const auto &target) {
    for (auto &obj : collection) {
      if (obj == target) {
        obj.inCluster = true;
        return;
      }
    }
  };

  // ########## TPs cluster #############
  std::vector<TriggerPrimitive> tps_in_chamber;
  // select TPs in chamber
  std::copy_if(tps.begin(), tps.end(), std::back_inserter(tps_in_chamber), inChamber);
  // find the higher quality in time
  auto bestTPIt{std::max_element(tps_in_chamber.begin(), tps_in_chamber.end(), compareTPs)};

  // the higherst quality in time is the bestTP of the cluster
  if (bestTPIt != tps_in_chamber.end() &&
      bestTPIt->BX == RIGHT_BX) {  // max_element returns last if the range is empty.
    _bestTP = *bestTPIt;
    toggleCluster(tps, _bestTP);
    tps_in_chamber.erase(bestTPIt);
    foundTP = true;
    clusterX = _bestTP.xLoc;

    if (_bestTP.BX != RIGHT_BX && tps_in_chamber.size() > 1) {
      std::cout << " la best TP di qualità " << _bestTP.quality << " non è in tempo, è al BX " << _bestTP.BX
                << " nel resto del cluster ho " << std::endl;
      for (auto const &T : tps_in_chamber) {
        if (T.index != _bestTP.index) {
          std::cout << " TP at BX " << T.BX << " and quality " << T.quality << std::endl;
        }
      }
      std::cout << "------------------------------------------------------------------------" << std::endl;
    }
  }

  if (tps_in_chamber.size() > 0) {
    // cluster close-by TPS from tps in chamber
    std::vector<TriggerPrimitive> cluster;

    if (_bestTP.index == 9999) {
      clusterX = tps_in_chamber[0].xLoc;
    }

    std::copy_if(tps_in_chamber.begin(), tps_in_chamber.end(), std::back_inserter(cluster), inRange);

    for (const auto &tp : cluster) {
      toggleCluster(tps, tp);
    }

    // sorting (actually partitioning) the vector:
    //  - first part is in-time TPs
    //  - second part is out-of-time TPs
    auto first_oot = std::partition(cluster.begin(), cluster.end(), [=](const auto &tp) { return tp.BX == RIGHT_BX; });

    // moving OOT TPs to dedicated std::vector and removing them from
    // tps_in_chamber
    std::move(first_oot, cluster.end(), std::back_inserter(_ootGhosts));
    cluster.erase(first_oot, cluster.end());

    _itGhosts = std::move(cluster);
  }

  // ########## Segments cluster #############

  std::vector<Segment> segments_in_chamber;
  std::copy_if(segs.begin(), segs.end(), std::back_inserter(segments_in_chamber), inChamber);

  // if I found TP cluster-> look for segment around it
  if (segments_in_chamber.size() > 0) {
    if (bestTPQuality() == -1) {
      auto bestSegIt = std::max_element(segments_in_chamber.begin(), segments_in_chamber.end(), compareSegs);
      _bestSeg = *bestSegIt;
      clusterX = _bestSeg.xLoc;
      foundSeg = true;
    }

    std::copy_if(segments_in_chamber.begin(), segments_in_chamber.end(), std::back_inserter(_segmentCluster), inRange);

    for (auto &seg : _segmentCluster) {
      toggleCluster(segs, seg);
    }

    if (segClusterSize() && bestTPQuality() > 0) {
      auto bestSegIt = std::max_element(_segmentCluster.begin(), _segmentCluster.end(), compareSegs);
      _bestSeg = *bestSegIt;
      foundSeg = true;
    }
  }

  // ########## Digi cluster #############
  std::vector<Digi> digis_in_chamber;
  std::copy_if(digis.begin(), digis.end(), std::back_inserter(digis_in_chamber),
               [&](const auto &digi) { return inChamber(digi) && digi.superlayer != 2; });

  auto first_sl3 = std::partition(digis_in_chamber.begin(), digis_in_chamber.end(),
                                  [=](const auto &digi) { return digi.superlayer == 1; });

  // if I found TP or segment cluster-> look for digis around it
  if (bestTPQuality() > 0 || bestSegPhiHits() > 0) {
    std::copy_if(digis_in_chamber.begin(), first_sl3, std::back_inserter(_digiCluster), inRangeDigi);

    if (_digiCluster.size() >= 10) {
      sl1Cluster = true;
    } else {
      _digiCluster.clear();
    }

    std::copy_if(first_sl3, digis_in_chamber.end(), std::back_inserter(_digiCluster), inRangeDigi);
    int digiClusterSize = _digiCluster.size();
    if ( digiClusterSize >= (10 + (10 * sl1Cluster))) {
      sl3Cluster = true;
    } else {
      _digiCluster.erase(_digiCluster.begin() + (10 * sl1Cluster), _digiCluster.end());
    }

  } else {
    std::vector<Digi> digis_in_chamber_sl1;
    std::vector<Digi> digis_in_chamber_sl3;

    std::move(digis_in_chamber.begin(), first_sl3, std::back_inserter(digis_in_chamber_sl1));
    std::move(first_sl3, digis_in_chamber.end(), std::back_inserter(digis_in_chamber_sl3));

    if (digis_in_chamber_sl1.size() > 10) {
      _digiCluster = digis_in_chamber_sl1[0].findCluster(digis_in_chamber_sl1, digiCut);
      if (_digiCluster.size() >= 10) {
        sl1Cluster = true;
      } else {
        _digiCluster.clear();
      }
    }

    if (digis_in_chamber_sl3.size() > 10) {
      auto tmpCluster = digis_in_chamber_sl3[0].findCluster(digis_in_chamber_sl3, digiCut);
      if (tmpCluster.size() >= 10) {
        sl3Cluster = true;
        std::move(tmpCluster.begin(), tmpCluster.end(), std::back_inserter(_digiCluster));
      }
    }
  }

  if (_digiCluster.size() > 0) foundDigi = true;

  for (const auto &digi : _digiCluster) {
    toggleCluster(digis, digi);
  }
}

int Cluster::ootSize() const { return _ootGhosts.size(); };
const std::vector<TriggerPrimitive> &Cluster::ootGhosts() const { return _ootGhosts; };

int Cluster::itSize() const { return _itGhosts.size(); };
const std::vector<TriggerPrimitive> &Cluster::itGhosts() const { return _itGhosts; };

int Cluster::tpClusterSize() const { return _ootGhosts.size() + _itGhosts.size() + 1; };

int Cluster::bestTPIndex() const { return _bestTP.index; };
int Cluster::bestTPQuality() const { return _bestTP.quality; };
const TriggerPrimitive &Cluster::bestTP() const { return _bestTP; };

int Cluster::itCountIf(std::function<bool(TriggerPrimitive const &)> f) const {
  return std::count_if(_itGhosts.begin(), _itGhosts.end(), f);
}

int Cluster::ootCountIf(std::function<bool(TriggerPrimitive const &)> f) const {
  return std::count_if(_ootGhosts.begin(), _ootGhosts.end(), f);
}

bool Cluster::matchMu(int muWh, int muStat, int muSec, double muXedge, double muYedge, double muX, int muIndex,
                      int iMu, double xCut) {
  if (!muMatched && muWh == wheel && muStat == station && muSec == sector && muXedge < -5 && muYedge < -5) {
    // if the extrapolated segment is within xCut from _bestTP
    if (foundSeg && std::abs(_bestSeg.xLoc - muX) < xCut) {
      muMatchedIndex[0] = iMu;
      muMatchedIndex[1] = muIndex;
      muMatched = true;
      return true;
    } else if (!foundSeg && foundTP && std::abs(_bestTP.xLoc - muX) < xCut) {
      muMatchedIndex[0] = iMu;
      muMatchedIndex[1] = muIndex;
      muMatched = true;
      return true;
    }
  }
  return false;
}

int Cluster::segClusterSize() const { return _segmentCluster.size(); };
const std::vector<Segment> &Cluster::segCluster() const { return _segmentCluster; };

int Cluster::bestSegIndex() const { return _bestSeg.index; };
int Cluster::bestSegPhiHits() const { return _bestSeg.nPhiHits; };
const Segment &Cluster::bestSeg() const { return _bestSeg; };

const std::vector<Digi> &Cluster::matchedDigi() const { return _matchedDigis; };

int Cluster::nDigi() const { return _digiCluster.size(); };

int Cluster::digiSL() const {
  if (sl1Cluster && sl3Cluster)
    return 5;
  else if (sl1Cluster)
    return 1;
  else if (sl3Cluster)
    return 3;
  else
    return 0;
};

int Cluster::MuIndex() const {
  if (muMatched) return muMatchedIndex[0];
  else return -1;
};
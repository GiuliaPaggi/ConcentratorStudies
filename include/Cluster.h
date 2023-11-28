#ifndef Cluster_h
#define Cluster_h

#include "include/Digi.h"
#include "include/Segment.h"
#include "include/TriggerPrimitive.h"

#include <iostream>
#include <algorithm>

// classe cluster-> qualità massima al bx giusto, n di ghost al bx giusto, n
// ghost al bx sbagliato, accedere a tutti (vettore indici di trigger primitive)
// parti da prompt
class Cluster {
 public:
  int wheel{-5};
  int sector{-1};
  int station{-1};  // CB can they be const?

  bool muMatched{false};
  std::array<int, 2> muMatchedIndex;  // 0 index of mu, 1 index of nmatch
  double xEdge{999.};
  double yEdge{999.};

  bool foundTP{false};  // CB what do we need out of this?
  bool foundSeg{false};
  bool foundDigi{false};

 private:
  Segment _bestSeg{};
  std::vector<Segment> _segmentCluster;

  TriggerPrimitive _bestTP{};
  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;
  TriggerPrimitive _earliestTP{};
  TriggerPrimitive _earliestTP_t0{};

  std::vector<Digi> _digiCluster;
  std::vector<Digi> _matchedDigis;

  bool sl1Cluster{false};
  bool sl3Cluster{false};

 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive>& tps, std::vector<Segment>& seg, std::vector<Digi>& digis, double xCut,
          double digiCut, int wh, int sec, int st);

  bool matchMu(int muWh, int muStat, int muSec, double muXedge, double muYedge, double muX, int muIndex, int nmu, double xCut);

  /// Size of TP clusters
  int itSize() const;
  int ootSize() const;
  int tpClusterSize() const;

  bool hasGhosts() const { return itSize() || ootSize(); };
  /// Provide list of in-time and out-of-time ghosts in cluster
  const std::vector<TriggerPrimitive>& itGhosts() const;
  const std::vector<TriggerPrimitive>& ootGhosts() const;

  /// Count number of in-time and out-of-time ghosts compatible with a given criteria
  int itCountIf(std::function<bool(TriggerPrimitive const&)> f) const;
  int ootCountIf(std::function<bool(TriggerPrimitive const&)> f) const;

  /// Information about the "best" TP in the cluster
  int bestTPIndex() const;
  int bestTPQuality() const;
  const TriggerPrimitive& bestTP() const;

  //earliest TP in cluster
  double earliestTPBX() const;
  double earliestTPt0() const;

  /// Size and vector of segment cluster
  int segClusterSize() const;
  const std::vector<Segment>& segCluster() const;

  /// Information about the "best" segment in the cluster
  int bestSegIndex() const;
  int bestSegPhiHits() const;
  const Segment& bestSeg() const;

  /// Size and SL of the digi cluster
  int nDigi() const;
  int digiSL() const;  // return 1 - 3 or 5 for both

  /// Vector of digis in the cluster
  const std::vector<Digi>& matchedDigi() const;

  //Index of matched Muon
  int MuIndex() const;
};

inline bool operator==(Cluster const& lCL, Cluster const& rCL) {
  return (lCL.wheel == rCL.wheel && lCL.sector == rCL.sector && lCL.station == rCL.station && lCL.foundTP &&
          rCL.foundTP && lCL.bestTP().quality == rCL.bestTP().quality);
  // else return(lCL.wheel == rCL.wheel && lCL.sector == rCL.sector && lCL.foundSeg && rCL.foundSeg && lCL.station ==
  // rCL.station && lCL.bestSegPhiHits() == rCL.bestSegPhiHits());
}

inline std::ostream& operator<<(std::ostream& os, Cluster const cluster) {
  os << "Cluster in wh " << cluster.wheel << " stat " << cluster.station << " sector " << cluster.sector
     << " has BestTP of quality " << cluster.bestTPQuality() << " in xLoc " << cluster.bestTP().xLoc << " at BX "
     << cluster.bestTP().BX << ", has " << cluster.segClusterSize() << " segments "
     << " in xLoc " << cluster.bestSeg().xLoc << " and has " << cluster.nDigi() << " digis " 
     << " is matched with muon " << cluster.muMatched << std::endl;
  return os;
}

#endif
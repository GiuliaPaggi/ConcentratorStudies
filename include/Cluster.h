#ifndef Cluster_h
#define Cluster_h

#include "include/TriggerPrimitive.h"
#include "include/Segment.h"
#include "include/Digi.h"

// classe cluster-> qualità massima al bx giusto, n di ghost al bx giusto, n
// ghost al bx sbagliato, accedere a tutti (vettore indici di trigger primitive)
// parti da prompt
class Cluster{
public:
  int wheel{-5};
  int sector{-1};
  int station{-1}; // CB can they be const?

  bool muMatched{false};
  std::array<int, 2> muMatchedIndex;

  bool foundTP{false}; // CB what do we need out of this?
  bool foundSeg{false};
  bool foundDigi{false};
  bool segMatched{false};
  bool digiMatched{false};

private:
  Segment _bestSeg{};
  std::vector<Segment> _segmentCluster;

  TriggerPrimitive _bestTP{};
  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;

  std::vector<Digi> _digiCluster;
  std::vector<Digi> _matchedDigis;

  bool sl1Cluster{false};
  bool sl3Cluster{false};


 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive> & tps, std::vector<Segment> & seg, std::vector<Digi> & digis, double xCut, double digiCut, int wh, int sec, int st);

  void matchMu( int muWh, int muStat, int muSec,  double muXedge, double muYedge, double muX, int muIndex, int nmu );
  void matchSegment(Segment segment, double xCut); // CB what is this for?
  void matchDigi(std::vector<Digi> const& digis, double xCut); // CB what is this for?

  /// Size of TP clusters
  int itSize() const;
  int ootSize() const;
  int tpClusterSize() const;

  /// ??
  bool isolated() const{ return bestTPQuality() < 0 && !hasGhosts(); }; // CB what is this?
  bool hasGhosts() const { return itSize() || ootSize(); };

  /// Provide list of in-time and out-of-time ghosts in cluster
  const std::vector<TriggerPrimitive>& itGhosts() const;
  const std::vector<TriggerPrimitive>& ootGhosts() const;

  /// Count number of in-time and out-of-time ghosts compatible with a given criteria
  int itCountIf(std::function<bool(TriggerPrimitive const&)> f) const;
  int ootCountIf(std::function<bool(TriggerPrimitive const&)>  f) const;

  /// Information about the "best" TP in the cluster
  int bestTPIndex() const;
  int bestTPQuality() const;
  const TriggerPrimitive& bestTP() const;

  /// Size and vector of segment cluster
  int segClusterSize() const;
  const std::vector<Segment> & segCluster() const;

  /// Information about the "best" segment in the cluster
  int bestSegIndex() const;
  int bestSegPhiHits() const;
  const Segment& bestSeg() const;

  /// Size and SL of the digi cluster
  int nDigi() const;
  int digiSL() const; // return 1 - 3 or 5 for both

  /// Vector of digis in the cluster
  const std::vector<Digi>& matchedDigi() const;


};

#endif
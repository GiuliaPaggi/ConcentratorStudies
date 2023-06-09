#ifndef Cluster_h
#define Cluster_h

#include "TriggerPrimitive.h"
#include "Segment.h"

class Cluster{
public:
  int wheel{-5};
  int station{-1};
  int sector{-1};

  bool muMatched = false;
  std::array<int, 2> muMatchedIndex;

  bool foundTP = false;
  bool foundSeg = false;
  bool foundDigi = false;
  bool segMatched = false;
  bool digiMatched = false;

private:
  Segment _bestSeg{};
  Segment _matchedSeg{};
  TriggerPrimitive _bestTP{};

  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;

  std::vector<Segment> _SegmentCluster;

  std::vector<Digi> _DigiCluster;
  std::vector<Digi> _matchedDigis;

  bool SL1Cluster = false;
  bool SL3Cluster = false;


 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive> & tps, std::vector<Segment> & seg, std::vector<Digi> & digis, double xCut, double digiCut, int wh, int sec, int st);

  void MatchMu( int muWh, int muStat, int muSec,  double muXedge, double muYedge, double muX, int muIndex, int nmu ) ;
  void MatchSegment(Segment segment, double xCut);

  int itSize() const;
  int ootSize() const;
  int tpClusterSize() const;

  bool isolated() const{ return bestTPQuality() < 0 && !hasGhosts(); };
  bool hasGhosts() const { return itSize() || ootSize(); };

  const std::vector<TriggerPrimitive>& itGhosts() const;
  const std::vector<TriggerPrimitive>& ootGhosts() const;

  int itCountIf(std::function<bool(TriggerPrimitive const&)> f) const;
  int ootCountIf(std::function<bool(TriggerPrimitive const&)>  f) const;

  int bestTPIndex() const;
  int bestTPQuality() const;
  const TriggerPrimitive& bestTP() const;

  int matchedSegIndex() const;
  int matchedSegPhiHits() const;
  const Segment& matchedSeg() const;
  int segClusterSize() const;
  const std::vector<Segment> segCluster() const;

  int bestSegPhiHits() const;
  const Segment& bestSeg() const;

  void MatchDigi(std::vector<Digi> const& digis, double xCut);
  const std::vector<Digi> matchedDigi() const;
  const int GetNDigi() const;
  const int WhichSL() const; // return 1 - 3 or 5 for both
  const double MeanDigiTime() const;
  const double MeanDigixLoc() const;

  const bool hasTP( TriggerPrimitive tp) const;
};

bool operator==(Cluster const &lCL, Cluster const &rCL) {
  return (lCL.wheel == rCL.wheel && lCL.sector == rCL.sector && lCL.station == rCL.station && lCL.foundTP && rCL.foundTP &&  lCL.bestTP().xLoc == rCL.bestTP().xLoc  );
}

#endif
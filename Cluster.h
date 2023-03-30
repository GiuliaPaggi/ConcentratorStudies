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

  bool segMatched = false;
  bool foundSeg = false;
  bool foundTP = false;
  bool digiMatched = false;

private:
  Segment _bestSeg{};
  Segment _matchedSeg{};
  TriggerPrimitive _bestTP{};

  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;

  std::vector<Digi> _matchedDigis;


 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive> const& tps, std::vector<Segment> const& seg, double xCut, int wh, int sec, int st);

  void MatchMu( int muWh, int muStat, int muSec,  double muXedge, double muYedge, double muX, int muIndex, int nmu ) ;
  void MatchSegment(Segment segment, double xCut);

  int itSize() const;
  int ootSize() const;

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

  int bestSegPhiHits() const;
  const Segment& bestSeg() const;

  void MatchDigi(std::vector<Digi> const& digis, double xCut);
  const std::vector<Digi> matchedDigi() const;

};

#endif
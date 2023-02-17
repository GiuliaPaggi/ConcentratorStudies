#ifndef Cluster_h
#define Cluster_h

#include "TriggerPrimitive.h"

class Cluster{
public:
  int wheel{-5};
  int station{-1};
  int sector{-1};

  bool muMatched = false;
  std::array<int, 2> muMatchedIndex;

 private:
  TriggerPrimitive _bestTP{};

  std::vector<TriggerPrimitive> _ootGhosts;
  std::vector<TriggerPrimitive> _itGhosts;



 public:
  Cluster(){};
  Cluster(std::vector<TriggerPrimitive> const& tps, double x_cut, int st, int wh, int sec);

  void MatchSegment( int muWh, int muStat, int muSec,  double muXedge, double muYedge, double muX, int muIndex, int nmu ) ;

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


};

#endif